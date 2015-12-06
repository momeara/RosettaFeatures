# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

library(ggplot2)


source("../hbond_geo_dim_scales.R")

feature_analyses <- c(feature_analyses, methods::new("FeaturesAnalysis",
id = "AHdist_by_resolution",
author = "Matthew O'Meara",
brief_description = "This measures the H-Bond A-H distance conditional on the resolution. Note that currently there is no features reporter for resolution so it must be included after the fact.",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){



sele <-"
SELECT
  geom.AHdist,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type,
  resolutions.resolution
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  hbond_sites AS don_site,
  hbond_sites AS acc_site,
  resolutions
WHERE
  hbond.struct_id = resolutions.struct_id AND
  hbond.struct_id = geom.struct_id AND
  hbond.hbond_id =  geom.hbond_id AND
  hbond.struct_id = don_site.struct_id AND
  hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND
  hbond.acc_id = acc_site.site_id;"
f <- query_sample_sources(sample_sources, sele)

f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type)
f <- na.omit(f, method="r")

plot_id <- "AHdist_by_resolution"
f_windows <- sliding_windows(f, c(), "resolution")
dens <- estimate_density_1d(
  f_windows, c("windows"),
  "AHdist", weight_fun = radial_3d_normalization)
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=windows)) +
	ggtitle("Hydrogen Bonds A-H Distance by Resolution Sliding Window") +
	scale_y_continuous("FeatureDensity") +
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.6))

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


plot_id <- "AHdist_by_resolution_by_chem_type"
f_windows <- sliding_windows(f, c("don_chem_type_name", "acc_chem_type_name", "sample_source"), "resolution")
dens <- estimate_density_1d(
  f_windows, c("don_chem_type_name", "acc_chem_type_name", "windows", "sample_source"),
  "AHdist", weight_fun = radial_3d_normalization)
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=windows)) +
	facet_grid(don_chem_type_name ~ acc_chem_type_name) +
	ggtitle("Hydrogen Bonds A-H Distance by Resolution Sliding Window") +
	scale_y_continuous("FeatureDensity") +
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.6))

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


plot_id <- "AHdist_by_resolution_by_chem_type_examples"
f_windows <- sliding_windows(
	f[
		f$acc_chem_type == "hbacc_CXL" &
		(f$don_chem_type == "hbdon_IME" |
		f$don_chem_type == "hbdon_HXL" |
		f$don_chem_type == "hbdon_CXA"),],
	c("don_chem_type_name", "sample_source"), "resolution")
dens <- estimate_density_1d(
  f_windows, c("don_chem_type_name", "windows", "sample_source"),
  "AHdist", weight_fun = radial_3d_normalization)
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=windows)) +
	facet_wrap( ~ don_chem_type_name, nrow=1) +
	ggtitle("Hydrogen Bonds A-H Distance by Resolution Sliding Window") +
	scale_y_continuous("FeatureDensity") +
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.6))

save_plots(self, plot_id, sample_sources, output_dir, output_formats)






})) # end FeaturesAnalysis
