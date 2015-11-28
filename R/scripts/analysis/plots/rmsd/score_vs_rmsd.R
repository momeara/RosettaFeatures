# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


check_setup()
feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "score_vs_rmsd",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("StructureScoreFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	s.score_value AS score,
	p.total_residue,
	r.all_atom AS rmsd
FROM
	score_types AS t,
	structure_scores AS s,
	pose_conformations AS p,
	protein_rmsd AS r
WHERE
	t.score_type_name == 'total_score' AND
	s.score_type_id = t.score_type_id AND
	p.struct_id = s.struct_id AND
	r.struct_id = s.struct_id;"

f <- query_sample_sources(sample_sources, sele)
f$score_per_residue <- f$score / f$total_residue

print(summary(f))

dens <- estimate_density_1d(
  data = f,
  ids = c("sample_source"),
  variable = "rmsd")

plot_id <- "rmsd_distribution"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source), size=1.4) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	ggtitle("All-Atom RMSD") +
	labs(x="Angstroms") +
	scale_y_continuous("FeatureDensity")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


f$log_rmsd <- log(f$rmsd)
dens <- estimate_density_1d(
  data = f,
  ids = c("sample_source"),
  variable = "log_rmsd")

plot_id <- "log_rmsd_distribution"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source), size=1.4) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	ggtitle("All-Atom Log(RMSD)") +
	labs(x="Log(Angstroms)") +
	scale_y_continuous("FeatureDensity")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



ddply(f, .(sample_source), function(df) {
	ss <- as.character(df$sample_source[1])
	plot_id <- paste("rmsd_vs_num_residues_", ss, sep="")
	p <- ggplot(data=df) + theme_bw() +
		geom_point(aes(x=total_residue, y=rmsd, colour=sample_source)) +
		ggtitle(paste("All-Atom RMSD vs Number of Residues SampleSource:", ss, sep="") ) +
		scale_x_continuous("Number of Residues") +
		scale_y_continuous("All Atom RMSD")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})


ddply(f, .(sample_source), function(df) {
	ss <- as.character(df$sample_source[1])
	plot_id <- paste("rmsd_vs_score_", ss, sep="")
	p <- ggplot(data=df) + theme_bw() +
		geom_point(aes(x=score, y=rmsd, colour=sample_source)) +
		ggtitle(paste("All-Atom RMSD vs Score SampleSource:", ss, sep="") ) +
		scale_x_continuous("Rosetta Score") +
		scale_y_continuous("All Atom RMSD")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})

ddply(f, .(sample_source), function(df) {
	ss <- as.character(df$sample_source[1])
	plot_id <- paste("rmsd_vs_score_per_residue_", ss, sep="")
	p <- ggplot(data=df) + theme_bw() +
		geom_point(aes(x=score_per_residue, y=rmsd, colour=sample_source)) +
		ggtitle(paste("All-Atom RMSD vs Score SampleSource:", ss, sep="") ) +
		scale_x_continuous("Rosetta Score Per Residue") +
		scale_y_continuous("All Atom RMSD")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})



summary_stats <- ddply(f, .(sample_source), function(df) {
	data.frame(
		mean = round(mean(df$rmsd),2),
		median = round(median(df$rmsd), 2),
		sd = round(sd(df$rmsd),2))
})

table_id <- "rmsd_summary_statistics"
table_title <- "All-Atom RMSD Summary Statistics"
save_tables(self,
	summary_stats, table_id, sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")

summary_stats <- ddply(f, .(sample_source), function(df) {
	data.frame(
		mean = round(mean(df$score_pre_residue),2),
		median = round(median(df$score_per_residue), 2),
		sd = round(sd(df$score_per_residue),2))
})


summary_stats <- ddply(f, .(sample_source), function(df) {
	data.frame(
		mean = round(mean(df$log_rmsd),2),
		median = round(median(df$log_rmsd), 2),
		sd = round(sd(df$log_rmsd),2))
})

table_id <- "log_rmsd_summary_statistics"
table_title <- "All-Atom Log(RMSD) Summary Statistics"
save_tables(self,
	summary_stats, table_id, sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")

summary_stats <- ddply(f, .(sample_source), function(df) {
	data.frame(
		mean = round(mean(df$score_pre_residue),2),
		median = round(median(df$score_per_residue), 2),
		sd = round(sd(df$score_per_residue),2))
})


table_id <- "total_score_per_residue_summary_statistics"
table_title <- "Total Score per Residue Summary Statistics"
save_tables(self,
	summary_stats, table_id, sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")



})) # end FeaturesAnalysis
