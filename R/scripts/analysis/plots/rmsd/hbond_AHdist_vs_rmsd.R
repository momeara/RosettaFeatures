# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

library(ggplot2)


feature_analyses <- c(feature_analyses, methods::new("FeaturesAnalysis",
id = "hbond_AHdist_vs_rmsd",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("StructureScoreFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
	geom.AHdist AS AHdist,
	r.all_atom AS rmsd,
	p.total_residue
FROM
	hbond_geom_coords AS geom,
	hbonds AS hb,
	hbond_sites AS don,
	hbond_sites AS acc,
	pose_conformations AS p,
	protein_rmsd AS r
WHERE
	geom.struct_id = hb.struct_id AND
	geom.hbond_id = hb.hbond_id AND
	p.struct_id = hb.struct_id AND
	r.struct_id = hb.struct_id AND
	don.struct_id = hb.struct_id AND
	don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND
	acc.site_id = hb.acc_id AND
	don.HBChemType != 'hbdon_PBA' AND
	acc.HBChemType != 'hbacc_PBA' AND
	ABS(don.resNum - acc.resNum) > 5;"

f <- query_sample_sources(sample_sources, sele)

f_sliding_windows <- sliding_windows(f, c(), "rmsd", verbose=T)

dens <- estimate_density_1d(
	data = f_sliding_windows,
	ids = c("sample_source", "windows"),
	variable = "AHdist", weight_fun = radial_3d_normalization)

plot_id <- "hbond_AHdist_by_rmsd_sliding_windows"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=windows)) +
	geom_indicator(aes(indicator=counts, colour=windows, group=windows)) +
	ggtitle("H-Bond AHdist by RMSD Sliding Windows") +
	facet_wrap(~sample_source) +
	scale_color_discrete("RMSD Windows") +
	scale_y_continuous("FeatureDensity") +
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.5))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
