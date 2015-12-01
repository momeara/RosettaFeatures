# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "structure_component_scores",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("StructureScoreFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
library(plyr)
library(ggplot2)

library(reshape2)
sele <-"
SELECT
	CAST( s.struct_id AS TEXT) AS struct_id,
	t.score_type_name,
	s.score_value,
	p.total_residue
FROM
	score_types AS t,
	score_types AS total_score_type,
	structure_scores AS total_score,
	structure_scores AS s,
	pose_conformations AS p
WHERE
	s.score_type_id = t.score_type_id AND
	p.struct_id = s.struct_id AND
	total_score_type.score_type_name = 'total_score' AND
	total_score.struct_id = s.struct_id AND
	total_score.score_type_id = total_score_type.score_type_id AND
	total_score.score_value < 20;"

f <- query_sample_sources(sample_sources, sele)

print(summary(f))

dens <- ddply(f, .(score_type_name), function(df){
	estimate_density_1d(
		data = df,
		ids = c("sample_source"),
		variable = "score_value")
})

plot_id <- "structure_component_scores"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source), size=1.4) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	ggtitle("Rosetta Structure Component Scores") +
	facet_wrap(~score_type_name, scales="free") +
	labs(x="Rosetta Energy Units") +
	scale_y_continuous("FeatureDensity")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

ddply(dens, .(score_type_name), function(sub_dens) {
	score_type_name <- as.character(sub_dens$score_type_name[1])
	plot_id <- paste("structure_component_scores", score_type_name, sep="_")

	p <- ggplot(data=sub_dens) + theme_bw() +
		geom_line(aes(x=x, y=y, colour=sample_source), size=1.4) +
		geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
		ggtitle(paste("Rosetta Structure Component Scores ScoreType: ", score_type_name, sep="")) +
		labs(x="Rosetta Energy Units") +
		scale_y_continuous("FeatureDensity")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})

f$score_per_residue <- f$score_value / f$total_residue

dens <- ddply(f, .(score_type_name), function(df) {
	estimate_density_1d(
		data = df,
		ids = c("sample_source"),
		variable = "score_per_residue")
})

summary_stats <- ddply(f, .(sample_source, score_type_name), function(df) {
	data.frame(
		mean = round(mean(df$score_per_residue),2),
		sd = round(sd(df$score_per_residue),2))
})

plot_id <- "structure_component_scores_per_residue"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source), size=1.4) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	geom_indicator(
		data=summary_stats,
		aes(indicator=mean, colour=sample_source, group=sample_source),
		xpos="left") +
	ggtitle("Rosetta Structure Per Residue Component Scores") +
	facet_wrap(~score_type_name, scales="free") +
	labs(x="Rosetta Energy Units") +
	scale_y_continuous("FeatureDensity")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

ddply(dens, .(score_type_name), function(sub_dens) {
	score_type_name <- as.character(sub_dens$score_type_name[1])
	print(score_type_name)
	plot_id <- paste("structure_component_scores_per_residue", score_type_name, sep="_")
	sub_summary_stats <- summary_stats[summary_stats$score_type_name == score_type_name,]

	p <- ggplot(data=sub_dens) + theme_bw() +
		geom_line(aes(x=x, y=y, colour=sample_source), size=1.4) +
		geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
		geom_indicator(
			data=sub_summary_stats,
			aes(indicator=mean, colour=sample_source, group=sample_source),
			xpos="left") +
		ggtitle(paste("Rosetta Structure Component Scores Per Residue ScoreType: ", score_type_name, sep="")) +
		labs(x="Rosetta Energy Units") +
		scale_y_continuous("FeatureDensity")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})


table_id <- "structure_component_scores_per_residue"
table_title <- "Rosetta Structure Per Residue Component Scores"
save_tables(self,
	summary_stats, table_id, sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")


##################################


f_wide <- cast(f,
	sample_source + total_residue + struct_id ~ score_type_name,
	value="score_value", fill=0L)
f_combined <- melt(
	with(
		f_wide,
		data.frame(
			sample_source = sample_source,
			struct_id = struct_id,
			LJ = (fa_rep + fa_atr) / total_residue,
			LK = fa_sol / total_residue,
			HBond = (hbond_bb_sc + hbond_lr_bb + hbond_sc + hbond_sr_bb) / total_residue,
			Coulomb = (fa_elec + fa_elec) / total_residue,
			Rotamer = fa_dun / total_residue,
			Backbone = (rama + p_aa_pp + omega) / total_residue)),
	id.vars=c("sample_source", "struct_id"), variable_name="score_type_name")
names(f_combined)[4] <- "score_per_residue"

dens <- ddply(f_combined, .(score_type_name), function(df) {
	estimate_density_1d(
		data = df,
		ids = c("sample_source"),
		variable = "score_per_residue")
})

summary_stats <- ddply(f_combined, .(sample_source, score_type_name), function(df) {
	data.frame(
		mean = round(mean(df$score_per_residue),2),
		sd = round(sd(df$score_per_residue),2))
})


plot_id <- "structure_combined_component_scores_per_residue"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source), size=1.4) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	geom_indicator(
		data=summary_stats,
		aes(indicator=mean, colour=sample_source, group=sample_source),
		xpos="left") +
	ggtitle("Rosetta Structure Per Residue Combined Component Scores") +
	facet_wrap(~score_type_name, scales="free", ncol=2) +
	labs(x="Rosetta Energy Units") +
	scale_y_continuous("FeatureDensity")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

ddply(dens, .(score_type_name), function(sub_dens) {
	score_type_name <- as.character(sub_dens$score_type_name[1])
	print(score_type_name)
	plot_id <- paste("structure_combined_component_scores_per_residue", score_type_name, sep="_")
	sub_summary_stats <- summary_stats[summary_stats$score_type_name == score_type_name,]

	p <- ggplot(data=sub_dens) + theme_bw() +
		geom_line(aes(x=x, y=y, colour=sample_source), size=1.4) +
		geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
		geom_indicator(
			data=sub_summary_stats,
			aes(indicator=mean, colour=sample_source, group=sample_source),
			xpos="left") +
		ggtitle(paste("Rosetta Structure Combined Component Scores Per Residue ScoreType: ", score_type_name, sep="")) +
		labs(x="Rosetta Energy Units") +
		scale_y_continuous("FeatureDensity")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})


table_id <- "structure_combined_component_scores_per_residue"
table_title <- "Rosetta Structure Per Residue Combined Component Scores"
save_tables(self,
	summary_stats, table_id, sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")


})) # end FeaturesAnalysis
