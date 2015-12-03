# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


library(ggplot2)


library(plyr)


feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "structure_scores",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("StructureScoreFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
	s.score_value,
	p.total_residue
FROM
	score_types AS t,
	structure_scores AS s,
	pose_conformations AS p
WHERE
	t.score_type_name == 'total_score' AND
	s.score_type_id = t.score_type_id AND
	score_value < 20 AND
	p.struct_id = s.struct_id;"

f <- query_sample_sources(sample_sources, sele)

print(summary(f))

dens <- estimate_density_1d(
  data = f,
  ids = c("sample_source"),
  variable = "score_value")

plot_id <- "structure_scores"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source), size=1.4) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	ggtitle("Rosetta Structure Scores") +
	labs(x="Rosetta Energy Units") +
	scale_y_continuous("FeatureDensity", breaks=c(0, .3, .6))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


f$energy_per_residue <- f$score_value / f$total_residue

dens <- estimate_density_1d(
  data = f,
  ids = c("sample_source"),
  variable = "energy_per_residue")

summary_stats <- ddply(f, .(sample_source), function(df) {
	data.frame(
		mean = round(mean(df$energy_per_residue),2),
		sd = round(sd(df$energy_per_residue),2))
})

plot_id <- "structure_scores_per_residue"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source), size=1.4) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	geom_indicator(
		data=summary_stats,
		aes(indicator=mean, colour=sample_source, group=sample_source),
		xpos="left") +
	ggtitle("Rosetta Structure Per Residue Scores") +
	labs(x="Rosetta Energy Units") +
	scale_y_continuous("FeatureDensity")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

table_id <- "structure_scores_per_residue"
table_title <- "Rosetta Structure Per Residue Scores"
save_tables(self,
	summary_stats, table_id, sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")



})) # end FeaturesAnalysis
