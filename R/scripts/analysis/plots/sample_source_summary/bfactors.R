# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "bfactors",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <- "
SELECT
	max_temperature
FROM
	residue_pdb_confidence;"
f <-  query_sample_sources(sample_sources, sele)

qs <- compute_quantiles(
	f, c("sample_source"), "max_temperature", 1000)

plot_id = "befactors"
p <- ggplot(data=qs) + theme_bw() +
	geom_line(aes(y=probs, x=quantiles, colour=sample_source)) +
	ggtitle(paste("Residue Maximum B-Factor", sep="")) +
	scale_y_continuous("Fraction", limit=c(0,1), breaks=c(0,.25, .50, .75)) +
	scale_x_continuous(expression(paste("B-Factor (", ring(A)^2, ")")), breaks=c(10, 20, 30, 50, 75, 100)) +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source), ypos=.9, yjust="bottom") +
	scale_colour_discrete("Sample Source") +
	coord_trans(x="sqrt", y="identity")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
