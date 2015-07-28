#! /usr/bin/env perl
# plot TOPP values in different ways
# Version 0: Jane Coates 28/7/2015

use strict;
use diagnostics;
use Statistics::R; 

my $R = Statistics::R->new();

$R->run(q` library(tidyr) `,
		q` library(dplyr) `,
		q` library(ggplot2) `,
		q` library(Cairo)  `,
        q` library(ggthemes) `,

		q` daily.data = read.table(file = "MOZART-4_TOPP_values.csv", sep = ",") `,
		q` colnames(daily.data) = c("VOC", "Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7") `,
		q` daily.data = daily.data %>% gather(Time, TOPP, -VOC) `,
		q` daily.data$Type = rep("Tagged O3", length(daily.data$TOPP))  `,

		q` daily.plot = ggplot(daily.data, aes(x = Time, y = TOPP, colour = VOC, group = VOC)) `,
		q` daily.plot = daily.plot + geom_line() `,
        q` daily.plot = daily.plot + geom_point()  `,
        q` daily.plot = daily.plot + scale_x_discrete(expand = c(0, 0.2)) `,
        q` daily.plot = daily.plot + scale_y_continuous(expand = c(0, 0.2)) `,
        q` daily.plot = daily.plot + theme_tufte() `,
        q` daily.plot = daily.plot + theme(axis.line = element_line(colour = "black")) `,
        q` daily.plot = daily.plot + theme(legend.title = element_blank()) `,
        q` daily.plot = daily.plot + theme(axis.title = element_text(face = "bold")) `,
        q` daily.plot = daily.plot + theme(axis.title.x = element_blank()) `,
        q` daily.plot = daily.plot + theme(axis.text.x = element_text(face = "bold")) `,

        q` CairoPDF(file = "Daily_TOPP_values.pdf", width = 10, height = 7) `,
        q` print(daily.plot) `,
        q` dev.off() `,

		q` cumulative.data = read.table(file = "MOZART-4_TOPP_cumulative_values.csv", sep = ",") `,
		q` colnames(cumulative.data) = c("VOC", "Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7") `,
		q` cumulative.data = cumulative.data %>% gather(Time, TOPP, -VOC) `,
		q` cumulative.data$Type = rep("Tagged O3", length(cumulative.data$TOPP))  `,

		q` cumulative.plot = ggplot(cumulative.data, aes(x = Time, y = TOPP, colour = VOC, group = VOC)) `,
		q` cumulative.plot = cumulative.plot + geom_line() `,
        q` cumulative.plot = cumulative.plot + geom_point()  `,
        q` cumulative.plot = cumulative.plot + scale_x_discrete(expand = c(0, 0.2)) `,
        q` cumulative.plot = cumulative.plot + scale_y_continuous(expand = c(0, 0.2)) `,
        q` cumulative.plot = cumulative.plot + theme_tufte() `,
        q` cumulative.plot = cumulative.plot + theme(axis.line = element_line(colour = "black")) `,
        q` cumulative.plot = cumulative.plot + theme(legend.title = element_blank()) `,
        q` cumulative.plot = cumulative.plot + theme(axis.title = element_text(face = "bold")) `,
        q` cumulative.plot = cumulative.plot + theme(axis.title.x = element_blank()) `,
        q` cumulative.plot = cumulative.plot + theme(axis.text.x = element_text(face = "bold")) `,

        q` CairoPDF(file = "Cumulative_TOPP_values.pdf", width = 10, height = 7) `,
        q` print(cumulative.plot) `,
        q` dev.off() `,

		q` original.daily = read.table(file = "MOZART-4_original_values.csv", sep = ",") `,
		q` colnames(original.daily) = c("VOC", "Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7") `,
		q` original.daily = original.daily %>% gather(Time, TOPP, -VOC) `,
		q` original.daily$Type = rep("Tagged Ox", length(original.daily$TOPP))  `,

		q` daily.comparison = rbind(daily.data, original.daily) `,

		q` daily.comparison = ggplot(daily.comparison, aes(x = Time, y = TOPP, colour = Type, group = Type)) `,
		q` daily.comparison = daily.comparison + geom_line() `,
        q` daily.comparison = daily.comparison + geom_point()  `,
        q` daily.comparison = daily.comparison + facet_wrap( ~ VOC, scales = "free") `,
        q` daily.comparison = daily.comparison + scale_x_discrete(expand = c(0, 0.2)) `,
        q` daily.comparison = daily.comparison + scale_y_continuous(expand = c(0, 0.2)) `,
        q` daily.comparison = daily.comparison + theme_tufte() `,
        q` daily.comparison = daily.comparison + theme(axis.line = element_line(colour = "black")) `,
        q` daily.comparison = daily.comparison + theme(legend.title = element_blank()) `,
        q` daily.comparison = daily.comparison + theme(axis.title = element_text(face = "bold")) `,
        q` daily.comparison = daily.comparison + theme(axis.title.x = element_blank()) `,
        q` daily.comparison = daily.comparison + theme(axis.text.x = element_text(face = "bold")) `,
        q` daily.comparison = daily.comparison + theme(strip.text = element_text(face = "bold")) `,
        q` daily.comparison = daily.comparison + theme(legend.position = "top") `,

        q` CairoPDF(file = "Daily_TOPP_values_comparison.pdf", width = 14, height = 10) `,
        q` print(daily.comparison) `,
        q` dev.off() `,

		q` original.cumul = read.table(file = "MOZART-4_original_cumulative_values.csv", sep = ",") `,
		q` colnames(original.cumul) = c("VOC", "Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7") `,
		q` original.cumul = original.cumul %>% gather(Time, TOPP, -VOC) `,
		q` original.cumul$Type = rep("Tagged Ox", length(original.cumul$TOPP))  `,

		q` cumulative.comparison = rbind(cumulative.data, original.cumul) `,

		q` cumulative.comparison = ggplot(cumulative.comparison, aes(x = Time, y = TOPP, colour = Type, group = Type)) `,
		q` cumulative.comparison = cumulative.comparison + geom_line() `,
        q` cumulative.comparison = cumulative.comparison + geom_point()  `,
        q` cumulative.comparison = cumulative.comparison + facet_wrap( ~ VOC, scales = "free") `,
        q` cumulative.comparison = cumulative.comparison + scale_x_discrete(expand = c(0, 0.2)) `,
        q` cumulative.comparison = cumulative.comparison + scale_y_continuous(expand = c(0, 0.2)) `,
        q` cumulative.comparison = cumulative.comparison + theme_tufte() `,
        q` cumulative.comparison = cumulative.comparison + theme(axis.line = element_line(colour = "black")) `,
        q` cumulative.comparison = cumulative.comparison + theme(legend.title = element_blank()) `,
        q` cumulative.comparison = cumulative.comparison + theme(axis.title = element_text(face = "bold")) `,
        q` cumulative.comparison = cumulative.comparison + theme(axis.title.x = element_blank()) `,
        q` cumulative.comparison = cumulative.comparison + theme(axis.text.x = element_text(face = "bold")) `,
        q` cumulative.comparison = cumulative.comparison + theme(strip.text = element_text(face = "bold")) `,
        q` cumulative.comparison = cumulative.comparison + theme(legend.position = "top") `,

        q` CairoPDF(file = "Cumulative_TOPP_values_comparison.pdf", width = 14, height = 10) `,
        q` print(cumulative.comparison) `,
        q` dev.off() `,
);
