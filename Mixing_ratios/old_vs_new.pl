#! /usr/bin/env perl
# Compare mixing ratios of species included in @ARGV between previous tagged runs and new tagged runs: 'real' chemistry only
# Version 0: Jane Coates 5/6/2015

use strict;
use diagnostics;
use MECCA;
use Statistics::R;

die "Specifiy species\n" if (@ARGV == 0);

my $base = "/local/home/coates/New_Tagging";
my @mechanisms = qw( MOZART-4 );
my @runs = qw( VOC old );
my %data;

my $mecca = MECCA->new("$base/MOZART-4_old_tagged/boxmodel");
my $time = $mecca->time;
$time -= $time->at(0);
$time /= 86400;

foreach my $mechanism (@mechanisms) {
    foreach my $run (@runs) {
        my $dir = "$base/${mechanism}_${run}_tagged/boxmodel";
        my $mecca = MECCA->new($dir);
        foreach my $spc (@ARGV) {
            if ($run =~ /old/) { 
                my @tagged_species = get_tagged_species ($spc);
                foreach my $new_spc (@tagged_species) {
                    my $tracer = $mecca->tracer($new_spc);
                    next unless (defined $tracer);
                    $data{$mechanism}{$run}{$spc} += $tracer;
                }
            } else {
                $data{$mechanism}{$run}{$spc} += $mecca->tracer($spc);
            }
        }
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(ggthemes) `,
        q` library(Cairo) `,
        q` library(grid) `,
);
$R->set('Time', [ map { $_ } $time->dog ]);
$R->run(q` d = data.frame() `);

foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    foreach my $run (sort keys %{$data{$mechanism}}) {
        $R->set('run', $run);
        $R->run(q` pre = data.frame(Time, Run = rep(run, length(Time))) `);
        foreach my $spc (sort keys %{$data{$mechanism}{$run}}) {
            $R->set('spc', $spc);
            $R->set('mixing.ratio', [ map { $_ } $data{$mechanism}{$run}{$spc}->dog ]);
            $R->run(q` pre[spc] = mixing.ratio `);
        }
        $R->run(q` pre = gather(pre, Species, Mixing.Ratio, -Time, -Run) `,
                q` d = rbind(d, pre) `,
        );
    }
}

$R->run(q` my.colours = c("VOC" = "#0352cb", "old" = "#ef6638") `);
$R->run(q` p = ggplot(d, aes(x = Time, y = Mixing.Ratio, colour = Run)) `,
        q` p = p + geom_line() `,
        q` p = p + facet_wrap( ~ Species, scales = "free_y") `,
        q` p = p + theme_hc() `,
        q` p = p + theme(legend.title = element_blank()) `,
        q` p = p + theme(legend.position = "top") `,
        q` p = p + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1), expand = c(0, 0)) `,
        q` p = p + theme(strip.text = element_text(face = "bold")) `,
        q` p = p + ylab("Mixing Ratio") `,
        q` p = p + xlab("Time (Days)") `,
        q` p = p + theme(panel.margin = unit(5, "mm")) `,
        q` p = p + scale_colour_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "old_vs_new_tagging_mixing_ratios.pdf", width = 10, height = 7) `,
        q` print(p) `,
        q` dev.off() `,
);

$R->stop();

sub get_tagged_species {
    my ($spc) = @_;

    my $spc_file = "$base/MOZART-4_old_tagged/gas.spc";
    open my $spc_in, '<:encoding(utf-8)', $spc_file or die $!;
    my @lines = <$spc_in>;
    close $spc_in;
    my @species;
    foreach my $line (@lines) {
        next unless ($line =~ /IGNORE/);
        chomp $line;
        my ($spc, $rest) = split / = /, $line;
        push @species, $spc;
    }
    my @tagged = ($spc);
    foreach my $item (@species) {
        push @tagged, $item if ($item =~ /^${spc}_/);
    }
    return @tagged;
}
