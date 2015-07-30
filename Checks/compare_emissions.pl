#! /usr/bin/env perl
# Compare emissions of VOC between O3 and Ox tagging runs
# Version 0: Jane Coates 29/7/2015

use strict;
use diagnostics;
use MECCA;
use KPP;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/New_Tagging";
my @runs = qw( O3 Ox );
my @VOC = qw( C2H6 C3H8 BIGALK C2H4 C3H6 BIGENE TOLUENE );
my %data;

foreach my $run (@runs) {
    my $dir = "$base/MOZART-4_tagged_$run";
    my $mecca = MECCA->new("$dir/boxmodel");
    my $kpp = KPP->new("$dir/gas.eqn");
    foreach my $voc (@VOC) {
        my $lookup;
        if ($run eq "O3") {
            $lookup = $voc . "_" . $voc;
        } else {
            $lookup = $voc;
        }
        my $reactions = $kpp->producing_from($lookup, "UNITY");
        foreach my $reaction (@$reactions) {
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number);
            $data{$run}{$voc} += $rate;
        }
    }
}

my $mecca = MECCA->new("$base/MOZART-4_tagged_O3/boxmodel");
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
        q` library(tidyr) `,
);

$R->set('Time', [ map { $_ } $times->dog ]);
$R->run(q` d = data.frame() `);
foreach my $run (sort keys %data) {
    $R->set('tagging', $run);
    $R->run(q` pre = data.frame(Time, Tagging = rep(tagging, length(Time))) `);
    foreach my $VOC (sort keys %{$data{$run}}) {
        $R->set('voc', $VOC);
        $R->set('rate', [ map { $_ } $data{$run}{$VOC}->dog ]);
        $R->run(q` pre[voc] = rate `);
    }
    $R->run(q` pre = pre[1:19,] `,
            q` pre = gather(pre, VOC, Rate, -Time, -Tagging) `,
            q` d = rbind(d, pre) `,
    );
}
#my $p = $R->run(q` print(d) `);
#print $p, "\n";

$R->run(q` p = ggplot(d, aes(x = Time, y = Rate, colour = Tagging)) `,
        q` p = p + geom_line() `,
        q` p = p + facet_wrap( ~ VOC, scales = "free") `,
        q` p = p + theme_tufte() `,
        q` p = p + theme(axis.line = element_line(colour = "black")) `,
        q` p = p + theme(legend.title = element_blank()) `,
        q` p = p + theme(legend.position = c(1, 0)) `,
        q` p = p + theme(legend.justification = c(1, 0)) `,
);

$R->run(q` CairoPDF(file = "Emission_rates_comparisons.pdf", width = 15, height = 10) `,
        q` print(p) `,
        q` dev.off() `,
);

$R->stop();
