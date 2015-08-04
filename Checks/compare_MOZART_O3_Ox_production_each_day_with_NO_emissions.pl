#! /usr/bin/env perl
# Compare the Ox and O3 production on each day with varying NO emissions
# Version 0: Jane Coates 4/8/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/MECCA/Max_Ox_Check/MOZART-4_tagged_";
my @fractions = qw( 0.5 0.8 0.9 1.0 1.1 1.2 1.5 );
my (%data, %families, %weights);

my $mecca = MECCA->new("${base}1.0/boxmodel");
my $ntime = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $n_per_day = 86400 / $dt;
my $n_days = int $ntime / $n_per_day;

foreach my $fraction (@fractions) {
    my $dir = "${base}$fraction";
    my $mecca = MECCA->new("$dir/boxmodel");
    my $kpp = KPP->new("$dir/gas.eqn");
    my $ro2_file = "$dir/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2_file);
    $families{"Ox"} = [ qw( O3 O O1D NO2 NO3 N2O5 HO2NO2 ), @no2_reservoirs ];
    $weights{"Ox"} = { NO3 => 2, N2O5 => 3 };
    $data{$fraction} = get_data($mecca, $kpp);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
        q` library(scales) `,
);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` d = data.frame() `);

foreach my $fraction (sort keys %data) {
    $R->set('fraction', $fraction);
    $R->run(q` pre = data.frame(Time, Fraction = rep(fraction, length(Time))) `);
    foreach my $spc (sort keys %{$data{$fraction}}) {
        $R->set('spc', $spc);
        $R->set('rate', [ map { $_ } $data{$fraction}{$spc}->dog ]);
        $R->run(q` pre[spc] = rate `);
    }
    $R->run(q` pre = gather(pre, Species, Rate, -Time, -Fraction) `,
            q` d = rbind(d, pre) `,
    );
}
#my $p = $R->run(q` print(d) `);
#print $p, "\n";
$R->run(q` p = ggplot(d, aes(x = Fraction, y = Rate, colour = Species)) `,
        q` p = p + geom_line() `,
        q` p = p + facet_wrap(~ Time, scales = "free", nrow = 1) `,
        q` p = p + scale_x_continuous(labels = percent) `,
        q` p = p + theme_tufte() `,
        q` p = p + theme(axis.title.x = element_blank()) `,
        q` p = p + theme(axis.line = element_line(colour = "black")) `,
        q` p = p + theme(strip.text = element_text(face = "bold")) `,
        q` p = p + theme(legend.title = element_blank()) `,
        q` p = p + theme(legend.position = "top") `,
        q` p = p + scale_colour_tableau() `,
);

$R->run(q` CairoPDF(file = "Max_O3_check_every_day.pdf", width = 15, height = 7) `,
        q` print(p) `,
        q` dev.off() `,
);

$R->stop();

sub get_no2_reservoirs { #get species that are produced when radical species react with NO2
    my ($kpp, $file) = @_; 
    open FILE, $file or die $!; 
    my @ro2;
    for (<FILE>) {
        push @ro2, split /\s+/, $_; 
    }
    close FILE;
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

sub remove_common_processes {
    my ($production, $consumption) = @_;
    my %common_processes;
    $common_processes{$_} = 1 for (grep { defined $production->{$_} } keys %$consumption) ;

    foreach my $process (keys %common_processes) {
        my $net_effect = $production->{$process} + $consumption->{$process};
        if ($net_effect->sum > 0) { #if net production remove consumption processes, retain net production
            if (which($net_effect < 0)->nelem > 0) { #which gets indices of non-0 values, nelem get nr of elements
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

sub get_data {
    my ($mecca, $kpp) = @_;
    my (%production_rates, %consumption_rates);
    foreach my $species (qw( O3 Ox )) {
        my ($producers, $producer_yields, $consumers, $consumer_yields);
        if (exists $families{$species}) {
            $kpp->family({
                            name    => $species,
                            members => $families{$species},
                            weights => $weights{$species},
            });
            $producers = $kpp->producing($species);
            $producer_yields = $kpp->effect_on($species, $producers);
            $consumers = $kpp->consuming($species);
            $consumer_yields = $kpp->effect_on($species, $consumers);
        } else {
            $producers = $kpp->producing($species);
            $producer_yields = $kpp->effect_on($species, $producers);
            $consumers = $kpp->consuming($species);
            $consumer_yields = $kpp->effect_on($species, $consumers);
        }
        print "No producers for $species\n" if (@$producers == 0);
        print "No consumers for $species\n" if (@$consumers == 0);

        for (0..$#$producers) {
            my $reaction = $producers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
            next if ($rate->sum == 0);
            $production_rates{$species} += $rate(1:$ntime-2);
        }

        for (0..$#$consumers) {
            my $reaction = $consumers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
            next if ($rate->sum == 0);
            $consumption_rates{$species} += $rate(1:$ntime-2);
        }
    }
    remove_common_processes(\%production_rates, \%consumption_rates);

    foreach my $spc (keys %production_rates) {
        my $reshape = $production_rates{$spc}->reshape($n_per_day, $n_days);
        $production_rates{$spc} = $reshape->sumover;
    }
    return \%production_rates;
}
