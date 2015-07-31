#! /usr/bin/env perl
# Allocate Ox budgets to initial VOC using same Ox family as in mechanism comparison, compare to previous Ox budget
# Compare O3 and Ox budgets for each VOC separately
# Version 0: Jane Coates 30/7/2015

use strict;
use diagnostics;
use MECCA;
use KPP;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/New_Tagging";
my @runs = qw( O3 Ox );
my %category_mapping = (
    MOZART  =>  {   BIGALK  => [ '0.285 NC4H10', '0.151 IC4H10', '0.146 NC5H12', '0.340 IC5H12', '0.048 NC6H14', '0.020 NC7H16', '0.010 NC8H18' ],
                    BIGENE  => [ '0.333 BUT1ENE', '0.667 MEPROPENE' ],
                    TOLUENE => [ '0.166 BENZENE', '0.478 TOLUENE_M', '0.073 EBENZ', '0.142 MXYL', '0.069 OXYL', '0.073 PXYL' ],
                }
);
my (%families, %weights, %data);
$families{"HO2x"} = [ qw( HO2 HO2NO2 ) ];

my $mecca = MECCA->new("$base/MOZART-4_tagged_O3/boxmodel");
my $ntime = $mecca->time->nelem;
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;
$times = $times(1:$ntime-2);

foreach my $run (@runs) {
    my $dir = "$base/MOZART-4_tagged_$run";
    my $mecca = MECCA->new("$dir/boxmodel");
    my $kpp = KPP->new("$dir/gas.eqn");
    my $spc_file = "$dir/gas.spc";
    my $RO2_file = "$dir/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $RO2_file);
    if ($run eq "O3") {
        $families{"Ox"} = [ qw( O3_X ) ];
        $families{"Ox"} = get_tagged_species($spc_file, "Ox");
    } else {
        $families{"Ox"} = [ qw( O3 NO2 O O1D NO3 N2O5 HO2NO2 ), @no2_reservoirs ];
        $weights{"Ox"} = { NO3 => 2, N2O5 => 3 };
    }
    $data{$run} = get_data($mecca, $kpp, $run);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
        q` library(grid) `,
);

$R->set('Time', [ map { $_ } $times->dog]);
$R->run(q` d = data.frame() `);
my @VOC = qw( C2H6 C3H8 BIGALK C2H4 C3H6 BIGENE ISOP TOLUENE );

foreach my $run (sort keys %data) {
    $R->set('run', $run);
    $R->run(q` pre = data.frame(Time, Run = rep(run, length(Time))) `);
    foreach my $item (sort keys %{$data{$run}}) {
        next unless ($item ~~ @VOC);
        $R->set('item', $item);
        $R->set('rate', [ map { $_ } $data{$run}{$item}->dog ]);
        $R->run(q` pre[item] = rate `);
    }
    $R->run(q` pre = gather(pre, Item, Rate, -Time, -Run) `,
            q` d = rbind(d, pre) `,
    );
}
#my $p = $R->run(q` print(d) `);
#print $p, "\n";

$R->run(q` p = ggplot(d, aes(x = Time, y = Rate, colour = Run)) `,
        q` p = p + geom_line() `,
        q` p = p + facet_wrap( ~ Item, scales = "free") `,
        q` p = p + ylab("Production Rate (molecules cm-3)") `,
        q` p = p + theme_few() `,
        q` p = p + scale_colour_tableau() `,
        q` p = p + theme(panel.border = element_blank()) `,
        q` p = p + theme(axis.line = element_line(colour = "grey")) `,
        q` p = p + theme(legend.title = element_blank()) `, 
        q` p = p + theme(legend.position = c( 0.8, 0.2)) `,
        q` p = p + theme(strip.text = element_text(face = "bold")) `,
        q` p = p + scale_y_continuous(expand = c(0, 0)) `,
        q` p = p + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1), expand = c(0, 0)) `,
);

$R->run(q` CairoPDF(file = "O3_allocated_production_budget_for_VOC.pdf", width = 10, height = 7) `,
        q` print(p) `,
        q` dev.off() `,
);

$R->stop();

sub remove_common_processes {
    my ($production, $consumption) = @_;
    my %common_processes;
    $common_processes{$_} = 1 for (grep { defined $production->{$_} } keys %$consumption) ;

    foreach my $process (keys %common_processes) {
        my $net_effect = $production->{$process} + $consumption->{$process};
        #print $process, $net_effect->nelem, "\n";
        if ($net_effect->sum > 0) { #if net production remove consumption processes, retain net production
            if (which($net_effect < 0)->nelem > 0) { #which gets indices of non-0 values, nelem get nr of elements
                #print "which if $process $net_effect\n";
                $production->{$process} .= $net_effect;
                $production->{$process}->where($net_effect < 0) .= 0; #where returns output pdl corresponding to input pdl and meeting the condition
                $consumption->{$process} .= $net_effect;
                $consumption->{$process}->where($net_effect > 0) .= 0; #where returns output pdl corresponding to input pdl and meeting the condition
                next;
            } 
            $production->{$process} = $net_effect;
            delete $consumption->{$process};
        } else { #net consumption
            if (which($net_effect > 0)->nelem > 0) {
                $production->{$process} .= $net_effect;
                $production->{$process}->where($net_effect < 0) .= 0;
                $consumption->{$process} .= $net_effect;
                $consumption->{$process}->where($net_effect > 0) .= 0;
                next;
            }
            $consumption->{$process} = $net_effect;
            delete $production->{$process};
        }
    }
} 

sub get_no2_reservoirs { #get species that are produced when radical species react with NO2
    my ($kpp, $file) = @_; 
    open my $in, '<:encoding(utf-8)', $file or die $!; 
    my @ro2;
    for (<$in>) {
        push @ro2, split /\s+/, $_; 
    }
    close $in;
    my @no2_reservoirs;
    foreach my $ro2 (@ro2) {
        my ($reactions) = $kpp->reacting_with($ro2, 'NO2');
        foreach my $reaction (@$reactions) {
            my ($products) = $kpp->products($reaction);
            if (@$products == 1) {
                push @no2_reservoirs, $products->[0];
            }   
        }   
    }   
    return @no2_reservoirs;
} 

sub get_tagged_species {
    my  ($file, $family) = @_;
    open my $spc_in, '<:encoding(utf-8)', $file or die $!;
    my @lines = <$spc_in>;
    close $spc_in;
    my @tagged_species;
    if (exists $families{$family}) {
        foreach my $spc (@{$families{$family}}) {
            foreach my $line (@lines) {
                next unless ($line =~ /^${spc}_/);
                my ($tagged, $rest) = split / = /, $line;
                push @tagged_species, $tagged;
            }
        }
    } else {
        print "No family for $family\n";
    }
    return \@tagged_species;
}

sub get_VOC {
    my ($tag) = @_;
    my $VOC;
    if ($tag eq "CH4") {
        $VOC = "Methane";
    } elsif ($tag eq "C2H6") {
        $VOC = "Ethane";
    } elsif ($tag eq "C3H8") {
        $VOC = "Propane";
    } elsif ($tag eq "NC4H10") {
        $VOC = "Butane";
    } elsif ($tag eq "IC4H10") {
        $VOC = "2-Methylpropane";
    } elsif ($tag eq "NC5H12") {
        $VOC = "Pentane";
    } elsif ($tag eq "IC5H12") {
        $VOC = "2-Methylbutane";
    } elsif ($tag eq "NC6H14") {
        $VOC = "Hexane";
    } elsif ($tag eq "NC7H16") {
        $VOC = "Heptane";
    } elsif ($tag eq "NC8H18") {
        $VOC = "Octane";
    } elsif ($tag eq "C2H4") {
        $VOC = "Ethene";
    } elsif ($tag eq "C3H6") {
        $VOC = "Propene";
    } elsif ($tag eq "BUT1ENE") {
        $VOC = "Butene";
    } elsif ($tag eq "MEPROPENE") {
        $VOC = "2-Methylpropene";
    } elsif ($tag eq "ISOP") {
        $VOC = "Isoprene";
    } elsif ($tag eq "BENZENE") {
        $VOC = "Benzene";
    } elsif ($tag eq "TOLUENE_M") {
        $VOC = "Toluene";
    } elsif ($tag eq "MXYL") {
        $VOC = "m-Xylene";
    } elsif ($tag eq "OXYL") {
        $VOC = "o-Xylene";
    } elsif ($tag eq "PXYL") {
        $VOC = "p-Xylene";
    } elsif ($tag eq "EBENZ") {
        $VOC = "Ethylbenzene";
    } elsif ($tag eq "INI" or $tag eq "CO" or $tag eq "XTR" or $tag eq "Others") {
        $VOC = $tag;
    } else {
        print "No VOC for $tag\n";
    }
    return $VOC;
}

sub get_data {
    my ($mecca, $kpp, $run) = @_;

    my (%production_rates, %consumption_rates);
    my $species = "Ox";
    $kpp->family({
            name    => $species,
            members => $families{$species},
            weights => $weights{$species},
    });
    my $producers = $kpp->producing($species);
    my $producer_yields = $kpp->effect_on($species, $producers);
    my $consumers = $kpp->consuming($species);
    my $consumer_yields = $kpp->effect_on($species, $consumers);
    print "No consumers for $species\n" if (@$consumers == 0);
    print "No producers for $species\n" if (@$producers == 0);

    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0);
        my ($number, $tag, $next) = split /_/, $reaction;
        my $label;
        if (defined $next) {
            $label = $next;
        } elsif (defined $tag) {
            $label = $tag;
        } else {
            $label = $kpp->reaction_string($reaction);
        } 
        $production_rates{$run}{$label} += $rate(1:$ntime-2);
    }

    for (0..$#$consumers) {
        my $reaction = $consumers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
        next if ($rate->sum == 0);
        my ($number, $tag, $next) = split /_/, $reaction;
        my $label;
        if (defined $next) {
            $label = $next;
        } elsif (defined $tag) {
            $label = $tag;
        } else {
            $label = $kpp->reaction_string($reaction);
        } 
        $consumption_rates{$run}{$label} += $rate(1:$ntime-2);
    }

    if ($run eq "Ox") { #allocated HO2x production to VOC
        $kpp->family({
                name    => "HO2x",
                members => $families{"HO2x"},
                weights => $weights{"HO2x"},
        });
        my $producers = $kpp->producing("HO2x");
        my $producer_yields = $kpp->effect_on("HO2x", $producers);
        my $consumers = $kpp->consuming("HO2x");
        my $consumer_yields = $kpp->effect_on("HO2x", $consumers);
        print "No consumers for HO2x\n" if (@$consumers == 0);
        print "No producers for HO2x\n" if (@$producers == 0);

        for (0..$#$producers) {
            my $reaction = $producers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
            next if ($rate->sum == 0);
            my ($number, $tag, $next) = split /_/, $reaction;
            my $label;
            if (defined $next) {
                $label = $next;
            } elsif (defined $tag) {
                $label = $tag;
            } else {
                $label = $kpp->reaction_string($reaction);
            } 
            $production_rates{"HO2x"}{$label} += $rate(1:$ntime-2);
        }

        for (0..$#$consumers) {
            my $reaction = $consumers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
            next if ($rate->sum == 0);
            my ($number, $tag, $next) = split /_/, $reaction;
            my $label;
            if (defined $next) {
                $label = $next;
            } elsif (defined $tag) {
                $label = $tag;
            } else {
                $label = $kpp->reaction_string($reaction);
            } 
            $consumption_rates{"HO2x"}{$label} += $rate(1:$ntime-2);
        }
        remove_common_processes($production_rates{"HO2x"}, $consumption_rates{"HO2x"});
        my $total_ho2x_production = zeroes(PDL::float, $ntime-2);
        $total_ho2x_production += $production_rates{"HO2x"}{$_} foreach (keys %{$production_rates{"HO2x"}});
        
        foreach my $reaction (keys %{$production_rates{'HO2x'}}) {
            $production_rates{$run}{$reaction} += $production_rates{$run}{'HO2 + NO = NO2 + OH'} * $production_rates{'HO2x'}{$reaction} / $total_ho2x_production;
            $consumption_rates{$run}{$reaction} += $consumption_rates{$run}{'HO2 + O3 = OH'} * $consumption_rates{'HO2x'}{$reaction} / $total_ho2x_production;
            $consumption_rates{$run}{$reaction} += $consumption_rates{$run}{'HO2 + NO3 = NO2 + OH'} * $consumption_rates{'HO2x'}{$reaction} / $total_ho2x_production;
        }
        delete $production_rates{$run}{'HO2 + NO = NO2 + OH'};
        delete $consumption_rates{$run}{'HO2 + O3 = OH'};
        delete $consumption_rates{$run}{'HO2 + NO3 = NO2 + OH'};

        delete $production_rates{"HO2x"};
        delete $consumption_rates{"HO2x"};
    }
    remove_common_processes($production_rates{$run}, $consumption_rates{$run});

    return $production_rates{$run};
}
