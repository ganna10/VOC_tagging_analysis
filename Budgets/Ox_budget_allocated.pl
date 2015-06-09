#! /usr/bin/env perl
# Allocate Ox budgets to initial VOC using same Ox family as in mechanism comparison, compare to previous Ox budget
# Version 0: Jane Coates 8/6/2015

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
my (%families, %weights, %production, %consumption);

my $mecca = MECCA->new("$base/MOZART-4_VOC_tagged/boxmodel");
my $ntime = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $n_per_day = 43200 / $dt;
my $n_days = int $ntime / $n_per_day;

foreach my $mechanism (@mechanisms) {
    my $dir = "$base/${mechanism}_VOC_tagged";
    my $mecca = MECCA->new("$dir/boxmodel");
    my $kpp = KPP->new("$dir/gas.eqn");
    my $RO2_file = "$dir/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $RO2_file);
    my @tagged_no2_reservoirs = grep { $_ .= "_X" } @no2_reservoirs;
    my $spc_file = "$dir/gas.spc";
    $families{"Ox"} = [ qw( O3_X O_X O1D_X NO2_X NO3_X NO2NO3_X HO2_X HO2NO2_X ), @tagged_no2_reservoirs ];
    $weights{"Ox"} = { NO3_X => 2, NO2NO3_X => 3 };
    ($production{$mechanism}, $consumption{$mechanism}) = get_data($mecca, $kpp, $spc_file);
}

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
        } else {
            $label = $tag;
        } 
        $production_rates{$label} = $rate(1:$ntime-2);
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
        } else {
            $label = $tag;
        } 
        $consumption_rates{$label} = $rate(1:$ntime-2);
    }
    remove_common_processes(\%production_rates, \%consumption_rates);
    return (\%production_rates, \%consumption_rates);
}
