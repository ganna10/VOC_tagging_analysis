#! /usr/bin/env perl
# compare mixing ratio time series between tagged and non-tagged Ox species
# Version 0: Jane Coates 5/6/2015

use strict;
use diagnostics;
use MECCA;
use PDL;
use Statistics::R;

my $base = "/local/home/coates/New_Tagging";
my @mechanisms = qw( MOZART-4 );
my %data;

foreach my $mechanism (@mechanisms) {
    my $mecca = MECCA->new("$base/${mechanism}_VOC_tagged/boxmodel");
    my @OX = get_tagged_Ox($mechanism);
    
    foreach my $ox (@OX) {
        my $mixing_ratio = $mecca->tracer($ox);
        next unless (defined $mixing_ratio);
        if ($ox =~ /_X_/) {
            my ($base, $tag) = split /_X_/, $ox;
            if ($base eq "NO2NO3" or $base eq "NO3NO2") {
                $base = "N2O5";
            } elsif ($base eq "NO2HO2") {
                $base = "HO2NO2";
            }
            $data{$mechanism}{"Tagged"}{$base} += $mixing_ratio;
        } else {
            next if ($ox =~ /_/);
            #print $ox, "\n";
            $data{$mechanism}{"Non-Tagged"}{$ox} += $mixing_ratio;
        }
    }
}

my $mecca = MECCA->new("$base/MOZART-4_VOC_tagged/boxmodel");
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(ggthemes) `,
        q` library(Cairo) `,
        q` library(grid) `,
);
$R->set('Time', [ map { $_ } $times->dog ]);
$R->run(q` d = data.frame() `);

foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    foreach my $item (sort keys %{$data{$mechanism}}) {
        $R->set('item', $item);
        $R->run(q` pre = data.frame(Time, Item = rep(item, length(Time))) `);
        foreach my $ox (sort keys %{$data{$mechanism}{$item}}) {
            $R->set('ox', $ox);
            $R->set('mixing.ratio', [ map { $_ } $data{$mechanism}{$item}{$ox}->dog ]);
            $R->run(q` pre[ox] = mixing.ratio `);
        }
        $R->run(q` pre = gather(pre, Ox, Mixing.Ratio, -Time, -Item) `,
                q` d = rbind(d, pre) `,
        );
    }
}

$R->run(q` my.colours = c("Tagged" = "#0352cb", "Non-Tagged" = "#ef6638") `);
$R->run(q` p = ggplot(d, aes(x = Time, y = Mixing.Ratio, colour = Item)) `,
        q` p = p + geom_line() `,
        q` p = p + facet_wrap( ~ Ox, scales = "free_y") `,
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

$R->run(q` CairoPDF(file = "Ox_tagged_non_tagged_mixing_ratios.pdf", width = 10, height = 7) `,
        q` print(p) `,
        q` dev.off() `,
);

$R->stop();

sub get_tagged_Ox {
    my ($mechanism) = @_;
    my @Ox = qw( O3 O O1D NO2 NO3 N2O5 NO3NO2 HNO3 HO2 HO2NO2 PAN MPAN ONIT ONITR ISOPNO3 ); ##only need one of the tagged N2O5 or HO2NO2 species
    open my $spc_in, '<:encoding(utf-8)', "$base/${mechanism}_VOC_tagged/gas.spc" or die $!;
    my @lines = <$spc_in>;
    close $spc_in;
    foreach my $line (@lines) {
        next unless ($line =~ /IGNORE/);
        foreach my $Ox (@Ox) {
            next unless ($line =~ /^${Ox}_X_/);
            my ($tagged_spc, $rest) = split / = /, $line;
            push @Ox, $tagged_spc;
        }
    }
    return @Ox;
}
