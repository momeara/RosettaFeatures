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


source("../hbond_geo_dim_scales.R")

feature_analyses <- c(feature_analyses, methods::new("FeaturesAnalysis",
id = "AHdist_chem_type",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){



sele <-"
SELECT
	geom.AHdist,
	don.HBChemType AS don_chem_type, acc.HBChemType AS acc_chem_type
FROM
	hbonds AS hb,
	hbond_geom_coords AS geom,
	hbond_sites AS don, hbond_sites AS acc,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
  abs(don.resNum - acc.resNum ) > 5;"

f <- query_sample_sources(sample_sources, sele)

f$don_chem_type_name <- don_chem_type_name_wrap(f$don_chem_type)
f$acc_chem_type_name <- acc_chem_type_name_wrap(f$acc_chem_type)
f <- na.omit(f, method="r")

plot_parts <- list(
	theme_bw(),
	geom_line(aes(x=x, y=y)),
	geom_indicator(aes(indicator=counts)),
	scale_y_continuous("FeatureDensity", breaks=c(1,3,5,7)),
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), breaks=c(1.6, 1.9, 2.2, 2.5)),
	coord_trans(limx=c(1.4,2.7), limy=c(0,7.5)))

dens <- estimate_density_1d(
	f, c("sample_source"),
	"AHdist", weight_fun = radial_3d_normalization)
plot_id <- "hbond_AHdist"
p <- ggplot(data=dens, aes(colour=sample_source, group=sample_source)) + plot_parts +
	ggtitle("H-Bond A-H Distance, SeqSep > 5, B-Factor < 30") +
	scale_colour_discrete("Sample Source")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d(
	f, c("sample_source", "don_chem_type"),
	"AHdist", weight_fun = radial_3d_normalization)
dens$don_chem_type_name <- don_chem_type_name_wrap(dens$don_chem_type)
plot_id <- "hbond_AHdist_don_chem_type"
p <- ggplot(data=dens, aes(colour=sample_source, group=sample_source)) + plot_parts +
	facet_wrap(~don_chem_type_name) +
	ggtitle("H-Bond A-H Distance by Donor Chemical Type, SeqSep > 5, B-Factor < 30") +
	scale_colour_discrete("Sample Source")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

dens$don_chem_type_name <- don_chem_type_name_linear(dens$don_chem_type)

plot_id <- "hbond_AHdist_don_chem_type_ss"
p <- ggplot(data=dens, aes(colour=don_chem_type_name, group=don_chem_type_name)) + plot_parts +
	facet_wrap(~sample_source, ncol=ceiling(sqrt(nrow(sample_sources)))) +
	ggtitle("H-Bond A-H Distance by Donor Chemical Type, SeqSep > 5, B-Factor < 30") +
	scale_colour_discrete("DonChemType")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


d_ply(dens, .(sample_source), function(sub_dens) {
	ss <- as.character(sub_dens$sample_source[1])
	plot_id <- paste("hbond_AHdist_don_chem_type", ss, sep="_")
	p <- ggplot(data=sub_dens, aes(colour=don_chem_type_name, group=don_chem_type_name)) + plot_parts +
		ggtitle(paste("H-Bond A-H Distance by Donor Chemical Type, SampleSource: ", ss, "\nSeqSep > 5, B-Factor < 30",sep="")) +
		scale_colour_discrete("DonChemType")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})



dens <- estimate_density_1d(
	f, c("sample_source", "acc_chem_type"),
	"AHdist", weight_fun = radial_3d_normalization)
dens$acc_chem_type_name <- acc_chem_type_name_wrap(dens$acc_chem_type)
plot_id <- "hbond_AHdist_acc_chem_type"
p <- ggplot(data=dens, aes(colour=sample_source, group=sample_source)) + plot_parts +
	facet_wrap(~acc_chem_type_name) +
	ggtitle("H-Bond A-H Distance by Acceptor Chemical Type, SeqSep > 5, B-Factor < 30") +
	scale_colour_discrete("Sample Source")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

dens$acc_chem_type_name <- acc_chem_type_name_linear(dens$acc_chem_type)

plot_id <- "hbond_AHdist_acc_chem_type_ss"
p <- ggplot(data=dens, aes(colour=acc_chem_type_name, group=acc_chem_type_name)) + plot_parts +
	facet_wrap(~sample_source, ncol=ceiling(sqrt(nrow(sample_sources)))) +
	ggtitle("H-Bond A-H Distance by Acceptor Chemical Type, SeqSep > 5, B-Factor < 30") +
	scale_colour_discrete("AccChemType")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

d_ply(dens, .(sample_source), function(sub_dens) {
	ss <- as.character(sub_dens$sample_source[1])
	plot_id <- paste("hbond_AHdist_acc_chem_type", ss, sep="_")
	p <- ggplot(data=sub_dens, aes(colour=acc_chem_type_name, group=acc_chem_type_name)) + plot_parts +
		ggtitle(paste("H-Bond A-H Distance by Acceptor Chemical Type, SampleSource: ", ss, "\nSeqSep > 5, B-Factor < 30", sep="")) +
		scale_colour_discrete("AccChemType")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})



f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type)

dens <- estimate_density_1d(
	f, c("sample_source", "acc_chem_type_name", "don_chem_type_name"),
	"AHdist", weight_fun = radial_3d_normalization)
plot_id <- "hbond_AHdist_chem_type"
p <- ggplot(data=dens, aes(colour=sample_source, group=sample_source)) + plot_parts +
	facet_grid(don_chem_type_name ~ acc_chem_type_name) +
	ggtitle("H-Bond A-H Distance by Chemical Type, SeqSep > 5, B-Factor < 30") +
	scale_colour_discrete("Sample Source")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
