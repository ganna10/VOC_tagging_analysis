#! /usr/bin/env perl
# TOPP calculation
# Version 0: Jane Coates 28/7/2015

use strict;
use diagnostics;
use MECCA;
use KPP;
use PDL;
use PDL::NiceSlice;

my $base = "/local/home/coates/New_Tagging";
my @mechanisms = qw( MOZART-4 );
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
    ($data{$mechanism}) = get_data($mecca, $kpp, $spc_file, $mechanism);
}

#export daily data
foreach my $mechanism (sort keys %data) {
    my $daily_out_file = "${mechanism}_TOPP_values.csv";
    open my $daily_out, '>:encoding(utf-8)', $daily_out_file or die $!;
    foreach my $VOC (sort keys %{$data{$mechanism}}) {
        print $daily_out "$VOC,";
        my @dailies = $data{$mechanism}{$VOC}->dog;
        print $daily_out join ',', @dailies;
        print $daily_out "\n";
    }
    close $daily_out;
}

#export cumulative data
foreach my $mechanism (sort keys %data) {
    my $cumulative_out_file = "${mechanism}_TOPP_cumulative_values.csv";
    open my $cumulative_out, '>:encoding(utf-8)', $cumulative_out_file or die $!;
    foreach my $VOC (sort keys %{$data{$mechanism}}) {
        print $cumulative_out "$VOC";
        my @dailies = $data{$mechanism}{$VOC}->dog;
        my $sum = 0;
        my @cumulative = ();
        foreach my $topp (@dailies) {
            $sum += $topp;
            print $cumulative_out ",$sum";
            push @cumulative, $sum;
        }
        print $cumulative_out "\n";
    }
    close $cumulative_out;
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
    my ($mecca, $kpp, $spc_file, $mechanism) = @_;

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

    my %emissions; #get emissions
    foreach my $VOC (sort keys %production_rates) {
        next if ($VOC eq "CO" or $VOC eq "XTR" or $VOC eq "INI" or $VOC eq "CH4");
        my $lookup = $VOC . "_" . $VOC;
        my $emission_reaction = $kpp->producing_from($VOC, "UNITY");
        next if (@$emission_reaction == 0);
        my $reaction_number = $kpp->reaction_number($emission_reaction->[0]);
        my $emission_rate = $mecca->rate($reaction_number);
        $emission_rate = $emission_rate(1:$ntime-2);
        $emissions{$VOC} = $emission_rate->sum * $dt;
    }

    my %TOPP; # calculate TOPPs
    foreach my $VOC (sort keys %production_rates) {
        next if ($VOC eq "CO" or $VOC eq "XTR" or $VOC eq "INI" or $VOC eq "CH4");
        my $rate = $production_rates{$VOC}->copy->reshape($n_per_day, $n_days);
        my $production = $rate->sumover * $dt;
        $TOPP{$VOC} = $production(0:13:2) / $emissions{$VOC};
    } 

    #correct TOPPs by carbon number and allocate to MCM species: TOPP actual VOC x real C number / lumped species C number
    $TOPP{"NC4H10"} = $TOPP{"BIGALK"} * 4 / 5;
    $TOPP{"IC4H10"} = $TOPP{"BIGALK"} * 4 / 5;
    $TOPP{"NC5H12"} = $TOPP{"BIGALK"} * 5 / 5;
    $TOPP{"IC5H12"} = $TOPP{"BIGALK"} * 5 / 5;
    $TOPP{"NC6H14"} = $TOPP{"BIGALK"} * 6 / 5;
    $TOPP{"NC7H16"} = $TOPP{"BIGALK"} * 7 / 5;
    $TOPP{"NC8H18"} = $TOPP{"BIGALK"} * 8 / 5;
    delete $TOPP{"BIGALK"};

    $TOPP{"BUT1ENE"} = $TOPP{"BIGENE"};
    $TOPP{"MEPROPENE"} = $TOPP{"BIGENE"};
    delete $TOPP{"BIGENE"};

    $TOPP{"BENZENE"} = $TOPP{"TOLUENE"} * 6 / 5;
    $TOPP{"MXYL"} = $TOPP{"TOLUENE"} * 8 / 7;
    $TOPP{"OXYL"} = $TOPP{"TOLUENE"} * 8 / 7;
    $TOPP{"PXYL"} = $TOPP{"TOLUENE"} * 8 / 7;
    $TOPP{"EBENZ"} = $TOPP{"TOLUENE"} * 8 / 7;

    return \%TOPP;
}
