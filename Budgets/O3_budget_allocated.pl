#! /usr/bin/env perl
# Allocate Ox budgets to initial VOC using same Ox family as in mechanism comparison, compare to previous Ox budget
# Version 0: Jane Coates 28/7/2015

use strict;
use diagnostics;
use MECCA;
use KPP;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/New_Tagging";
my @mechanisms = qw( MOZART-4 );
my %category_mapping = (
    MOZART  =>  {   BIGALK  => [ '0.285 NC4H10', '0.151 IC4H10', '0.146 NC5H12', '0.340 IC5H12', '0.048 NC6H14', '0.020 NC7H16', '0.010 NC8H18' ],
                    BIGENE  => [ '0.333 BUT1ENE', '0.667 MEPROPENE' ],
                    TOLUENE => [ '0.166 BENZENE', '0.478 TOLUENE_M', '0.073 EBENZ', '0.142 MXYL', '0.069 OXYL', '0.073 PXYL' ],
                }
);
my (%families, %weights, %data);

my $mecca = MECCA->new("$base/MOZART-4_VOC_tagged/boxmodel");
my $ntime = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $n_per_day = 43200 / $dt;
my $n_days = int $ntime / $n_per_day;

foreach my $mechanism (@mechanisms) {
    my $dir = "$base/${mechanism}_VOC_tagged";
    my $mecca = MECCA->new("$dir/boxmodel");
    my $kpp = KPP->new("$dir/gas.eqn");
    my $spc_file = "$dir/gas.spc";
    $families{"Ox"} = [ qw( O3_X ) ];
    ($data{$mechanism}) = get_data($mecca, $kpp, $spc_file);
}

foreach my $mechanism (keys %data) {
    my $key;
    if ($mechanism eq "MOZART-4") {
        $key = "MOZART";
    } else {
        $key = $mechanism;
    }
    if (exists $category_mapping{$key}) {
        my %allocated;
        foreach my $VOC (sort keys %{$data{$mechanism}}) {
            if (exists $category_mapping{$key}{$VOC}) {
                foreach my $item (@{$category_mapping{$key}{$VOC}}) {
                    my ($fraction, $mcm) = split / /, $item;
                    $allocated{$mcm} = $fraction * $data{$mechanism}{$VOC};
                }
            } else {
                $allocated{$VOC} = $data{$mechanism}{$VOC};
            }
        }
        $data{$mechanism} = \%allocated;
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
        q` library(grid) `,
);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` d = data.frame() `);

foreach my $mechanism (sort keys %data) { #hash from get_data subroutine
    $R->set('mechanism', $mechanism);
    $R->run(q` pre = data.frame(Time) `);
    foreach my $voc (sort keys %{$data{$mechanism}}) {
        my $VOC = get_VOC($voc);
        $R->set('tag', $VOC);
        $R->set('rate', [ map { $_ } $data{$mechanism}{$voc}->dog ]);
        $R->run(q` pre[tag] = rate `);
    }
    $R->run(q` pre = gather(pre, VOC, Rate, -Time) `,
            q` d = rbind(d, pre) `,
    );
}
$R->run(q` d$VOC = factor(d$VOC, levels = c("CO", "Methane", "Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane", "Ethene", "Isoprene", "Benzene", "Toluene", "m-Xylene", "o-Xylene", "p-Xylene", "Others")) `,
        #q` d$VOC = factor(d$VOC, levels = c("XTR", "CO", "Methane", "Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane", "Ethene", "Propene", "Butene", "2-Methylpropene", "Isoprene", "Benzene", "Toluene", "m-Xylene", "o-Xylene", "p-Xylene", "Ethylbenzene")) `,
        q` my.colours = c(  "CO" = "#6c254f", 
                            "Methane" = "#f9c500", 
                            "Ethane" = "#0e5c28", 
                            "Propane" = "#2b9eb3", 
                            "Butane" = "#ef6638", 
                            "2-Methylpropane" = "#0352cb",
                            "Pentane" = "#86b650",
                            "2-Methylbutane" = "#8c1531",
                            "Hexane" = "#77aecc",
                            "Heptane" = "#f7c56c",
                            "Octane" = "#c9a415",
                            "Ethene" = "#f36a71",
                            "Propene" = "#f7c56c",
                            "Butene" = "#6db875",
                            "2-Methylpropene" = "#f3aa7f",
                            "Isoprene" = "#8ed6d2",
                            "Benzene" = "#000000",
                            "Toluene" = "#b569b3",
                            "m-Xylene" = "#e7e85e",
                            "o-Xylene" = "#ae4901",
                            "p-Xylene" = "#9bb18d",
                            "Ethylbenzene" = "#f36a71",
                            "Others" = "#58591b") `,
);
#my $p = $R->run(q` print(d) `);
#print $p, "\n";

$R->run(q` p = ggplot(d, aes(x = Time, y = Rate, fill = VOC, order = VOC)) `,
        q` p = p + geom_bar(stat = "identity", position = "stack") `,
        q` p = p + ylab("O3 Production Rate (molecules cm-3)") `,
        q` p = p + theme_tufte() `,
        q` p = p + theme(axis.line = element_line(colour = "black")) `,
        q` p = p + theme(legend.title = element_blank()) `,
        q` p = p + theme(axis.title.x = element_blank()) `,
        q` p = p + scale_y_continuous(expand = c(0, 0)) `,
        q` p = p + scale_x_discrete(expand = c(0, 0)) `,
        q` p = p + theme(axis.title = element_text(face = "bold")) `,
        q` p = p + theme(axis.title.x = element_blank()) `,
        q` p = p + theme(axis.text.x = element_text(face = "bold")) `,
        q` p = p + scale_fill_manual(values = my.colours, limits = rev(levels(d$VOC))) `,
);

$R->run(q` CairoPDF(file = "O3_allocated_production_budget.pdf", width = 10, height = 7) `,
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
    my ($mecca, $kpp, $spc_file) = @_;

    my (%production_rates, %consumption_rates);
    my $species = "Ox";
    $families{$species} = get_tagged_species($spc_file, $species);
    my ($producers, $producer_yields, $consumers, $consumer_yields);
    $kpp->family({
            name    => $species,
            members => $families{$species},
            weights => $weights{$species},
    });
    $producers = $kpp->producing($species);
    $producer_yields = $kpp->effect_on($species, $producers);
    $consumers = $kpp->consuming($species);
    $consumer_yields = $kpp->effect_on($species, $consumers);
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
            print "no label for $reaction\n";
        } 
        $production_rates{$label} += $rate(1:$ntime-2);
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
            print "no label for consuming $reaction\n";
        } 
        $consumption_rates{$label} += $rate(1:$ntime-2);
    }
    remove_common_processes(\%production_rates, \%consumption_rates);

    my $others = 1e8;
    foreach my $item (keys %production_rates) {
        my $reshape = $production_rates{$item}->reshape($n_per_day, $n_days);
        $reshape = $reshape->sumover;
        $production_rates{$item} = $reshape(0:13:2);
        if ($production_rates{$item}->sum < $others) {
            $production_rates{"Others"} += $production_rates{$item};
            delete $production_rates{$item};
        }
    }
    return \%production_rates;
}
