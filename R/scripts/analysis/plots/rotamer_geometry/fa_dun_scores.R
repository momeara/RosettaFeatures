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

feature_analyses <- c(feature_analyses, methods::new("FeaturesAnalysis",
id = "fa_dun_scores",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueScoresFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele_1b <-"
SELECT
	r.name3 AS res_type,
	rs.score_value
FROM
	residues AS r,
	residue_scores_1b AS rs,
	score_types AS st
WHERE
	rs.struct_id = r.struct_id AND
	r.resNum = rs.resNum AND
	rs.score_type_id = st.score_type_id AND
	st.score_type_name = 'fa_dun' AND
	-2 < score_value AND score_value < 5;"

scores <- query_sample_sources(sample_sources, sele_1b)
scores$res_type <- factor(scores$res_type)

dens <- estimate_density_1d(
	data = scores,
	ids = c("sample_source", "res_type"),
	variable = "score_value")

plot_id <- "fa_dun_by_res_type"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source), size=2) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	ggtitle("Rosetta Dunbrack Energy") +
	facet_wrap(~res_type);
	labs(x="Rosetta Energy Units") +
	scale_y_continuous("FeatureDensity")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)

d_ply(dens, .(res_type), function(sub_dens){
		res_type <- sub_dens[1, "res_type"]
				print(res_type)
	plot_id <- paste("fa_dun_", res_type, sep="")
	p <- ggplot(data=sub_dens) + theme_bw() +
		geom_line(aes(x=x, y=y, colour=sample_source), size=2) +
		geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
		ggtitle(paste("Rosetta Dunbrack Energy for ", res_type, sep="")) +
		labs(x="Rosetta Energy Units") +
		scale_y_continuous("FeatureDensity")
				if(nrow(sample_sources) <= 3){
		p <- p + theme(legend.position="bottom", legend.direction="horizontal")
	}

	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})

})) # end FeaturesAnalysis
