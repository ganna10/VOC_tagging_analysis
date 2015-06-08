#! /usr/bin/env perl
# Time series of all tagged HO2 mixing ratios stacked to give total HO2
# Version 0: Jane Coates 8/6/2015

use strict;
use diagnostics;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/work/users/jco/New_Tagging";
#my $base = "/local/home/coates/New_Tagging";
my @mechanisms = qw( MOZART-4 );
my $species = "HO2";
my %data;

my %category_mapping = (
    MOZART  =>  {   BIGALK  => [ '0.285 NC4H10', '0.151 IC4H10', '0.146 NC5H12', '0.340 IC5H12', '0.048 NC6H14', '0.020 NC7H16', '0.010 NC8H18' ],
                    BIGENE  => [ '0.333 BUT1ENE', '0.667 MEPROPENE' ],
                    TOLUENE => [ '0.166 BENZENE', '0.478 TOLUENE_M', '0.073 EBENZ', '0.142 MXYL', '0.069 OXYL', '0.073 PXYL' ],
                }
);

my $mecca = MECCA->new("$base/MOZART-4_VOC_tagged/boxmodel");
my $ntime = $mecca->time->nelem;
my $time = $mecca->time;
$time -= $time->at(0);
$time /= 86400;
$time = $time(1:$ntime-2);

foreach my $mechanism (@mechanisms) {
    my $dir = "$base/${mechanism}_VOC_tagged";
    my $mecca = MECCA->new("$dir/boxmodel");
    my @tagged_species = get_tagged_species($dir, $species);
    foreach my $spc (@tagged_species) {
        my ($base, $tag) = split /_X_/, $spc;
        my $mixing_ratio = $mecca->tracer($spc);
        $mixing_ratio = $mixing_ratio(1:$ntime-2);
        $data{$mechanism}{$tag} = get_contributions($tag, $mechanism, $mixing_ratio);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `,
        q` library(ggthemes) `,
);

$R->set('Time', [ map { $_ } $time->dog ]);
$R->run(q` d = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    $R->run(q` pre = data.frame(Time) `);
    foreach my $tag (sort keys %{$data{$mechanism}}) {
        foreach my $voc (sort keys %{$data{$mechanism}{$tag}}) {
            $R->set('voc', $voc);
            $R->set('mixing.ratio', [ map { $_ } $data{$mechanism}{$tag}{$voc}->dog ]);
            $R->run(q` pre[voc] = mixing.ratio `);
        }
    }
    $R->run(q` pre = gather(pre, VOC, Mixing.Ratio, -Time) `,
            q` d = rbind(d, pre) `,
    );
}

$R->run(q` d$VOC = factor(d$VOC, levels = c("INI", "XTR", "CO", "Methane", "Ethane", "Propane", "Butane", "2-Methylpropane", "Pentane", "2-Methylbutane", "Hexane", "Heptane", "Octane", "Ethene", "Propene", "Butene", "2-Methylpropene", "Isoprene", "Benzene", "Toluene", "m-Xylene", "o-Xylene", "p-Xylene", "Ethylbenzene")) `,
        q` my.colours = c( "Ethane" = "#696537", "Propane" = "#f9c500", "Butane" = "#76afca", "2-Methylpropane" = "#dc3522", "Pentane" = "#8c6238", "2-Methylbutane" = "#9bb08f", "Hexane" = "#8b1537", "Heptane" = "#ba8b01", "Octane" = "#0352cb", "Ethene" = "#86b650", "Propene" = "#6c254f", "Butene" = "#ee6738", "2-Methylpropene" = "#58691b", "Isoprene" = "#8ed6d5", "Benzene" = "#f3aa7f", "Toluene" = "#c65d6c", "m-Xylene" = "#888a87", "o-Xylene" = "#0e5c28", "p-Xylene" = "#b569b3", "Ethylbenzene" = "#2c9def", "Methane" = "#000000", "INI" = "#c9a415", "XTR" = "#1c3e3d", "CO" = "#d94a80" ) `,
);

$R->run(q` p = ggplot(d, aes(x = Time, y = Mixing.Ratio, fill = VOC, order = VOC)) `,
        q` p = p + geom_area(position = "stack") `,
        q` p = p + geom_area(position = "stack", colour = "black", show_guide = FALSE) `,
        q` p = p + theme_tufte() `,
        q` p = p + theme(axis.line = element_line(colour = "black")) `,
        q` p = p + theme(legend.title = element_blank()) `,
        q` p = p + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1), expand = c(0, 0)) `,
        q` p = p + scale_y_continuous(expand = c(0, 1e-19)) `,
        q` p = p + theme(strip.text = element_text(face = "bold")) `,
        q` p = p + ylab("Mixing Ratio") `,
        q` p = p + xlab("Time (Days)") `,
        q` p = p + scale_fill_manual(values = my.colours, limits = rev(levels(d$VOC))) `,
);

$R->run(q` CairoPDF(file = "HO2_mixing_ratio_components.pdf", width = 10, height = 7) `,
        q` print(p) `,
        q` dev.off() `,
);

$R->stop();

sub get_tagged_species {
    my ($dir, $species) = @_;
    my $spc_file = "$dir/gas.spc";
    my @tagged;
    open my $in, '<:encoding(utf-8)', $spc_file or die $!;
    while (<$in>) {
        next unless ($_ =~ /^${species}_/);
        my ($spc, $rest) = split / = /, $_;
        push @tagged, $spc;
    }
    close $in;
    return @tagged;
}

sub get_contributions {
    my ($spc, $mechanism, $mixing_ratio) = @_;
    $mechanism = "MOZART" if ($mechanism =~ /MOZ/);

    my %data;
    if (exists $category_mapping{$mechanism}) { 
        if (exists $category_mapping{$mechanism}{$spc}) {
            my $list = $category_mapping{$mechanism}{$spc}; 
            foreach (@$list) {
                my ($fraction, $VOC) = split /\s/, $_;
                $VOC = get_voc_name($VOC);
                $data{$VOC} = $fraction * $mixing_ratio;
            }
        } else {
            my $VOC = get_voc_name($spc);
            $data{$VOC} = $mixing_ratio;
        }
    }
    return \%data;
}

sub get_voc_name {
    my ($spc) = @_;
    my $VOC;
    if ($spc eq "C2H6") {
        $VOC = "Ethane";
    } elsif ($spc eq "C3H8") {
        $VOC = "Propane";
    } elsif ($spc eq "NC4H10") {
        $VOC = "Butane";
    } elsif ($spc eq "IC4H10") {
        $VOC = "2-Methylpropane";
    } elsif ($spc eq "NC5H12") {
        $VOC = "Pentane";
    } elsif ($spc eq "IC5H12") {
        $VOC = "2-Methylbutane";
    } elsif ($spc eq "NC6H14") {
        $VOC = "Hexane";
    } elsif ($spc eq "NC7H16") {
        $VOC = "Heptane";
    } elsif ($spc eq "NC8H18") {
        $VOC = "Octane";
    } elsif ($spc eq "C2H4") {
        $VOC = "Ethene";
    } elsif ($spc eq "C3H6") {
        $VOC = "Propene";
    } elsif ($spc eq "BUT1ENE") {
        $VOC = "Butene";
    } elsif ($spc eq "MEPROPENE") {
        $VOC = "2-Methylpropene";
    } elsif ($spc eq "ISOP") {
        $VOC = "Isoprene";
    } elsif ($spc eq "BENZENE") {
        $VOC = "Benzene";
    } elsif ($spc eq "TOLUENE" or $spc eq "TOLUENE_M") {
        $VOC = "Toluene";
    } elsif ($spc eq "MXYL") {
        $VOC = "m-Xylene";
    } elsif ($spc eq "OXYL") {
        $VOC = "o-Xylene";
    } elsif ($spc eq "PXYL") {
        $VOC = "p-Xylene";
    } elsif ($spc eq "EBENZ") {
        $VOC = "Ethylbenzene";
    } elsif ($spc eq "CH4") {
        $VOC = "Methane";
    } elsif ($spc eq "XTR" or $spc eq "INI" or $spc eq "CO") {
        $VOC = $spc;
    } else {
        print "No mapping for $spc\n";
    }
    return $VOC;
}
