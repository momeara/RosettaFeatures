# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "AHD_cdf_chem_type",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")

sele <-"
SELECT
  geom.cosAHD,
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

f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type)
f <- na.omit(f, method="r")


f$AHD <- acos(f$cosAHD) * 180/pi




qs <- compute_quantiles(
	f, c("sample_source"), "AHD", 1000)
plot_id = "hbond_AHD_CDF"
p <- ggplot(data=qs) + theme_bw() +
	geom_line(aes(y=probs, x=quantiles, colour=sample_source)) +
	ggtitle(paste("H-Bond AHD Cumulative Distribution Function, SeqSep > 5, B-Fact < 30", sep="")) +
	scale_y_continuous("Fraction", limit=c(0,1), breaks=c(0,.25, .50, .75)) +
	scale_x_continuous("Angle Deviation from Linear for Acceptor -- Hydrogen -- Donor (degrees)") +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source), ypos=.9, yjust="bottom") +
	scale_colour_discrete("Sample Source")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



qs <- compute_quantiles(
	f, c("sample_source", "don_chem_type"), "AHD", 1000)
qs$don_chem_type_name <- don_chem_type_name_wrap(qs$don_chem_type)

plot_id = "hbond_AHD_CDF_don_chem_type"
p <- ggplot(data=qs) + theme_bw() +
	geom_line(aes(y=probs, x=quantiles, colour=sample_source)) +
	ggtitle(paste("H-Bond AHD Cumulative Distribution Function by Donor Chemical Type, SeqSep > 5, B-Fact < 30", sep="")) +
	scale_y_continuous("Fraction", limit=c(0,1), breaks=c(0,.25, .50, .75)) +
	scale_x_continuous("Angle Deviation from Linear for Acceptor -- Hydrogen -- Donor (degrees)") +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source), ypos=.7, yjust="bottom") +
	facet_wrap(~don_chem_type_name) +
	scale_colour_discrete("Sample Source")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

qs$don_chem_type_name <- don_chem_type_name_linear(qs$don_chem_type)
plot_id = "hbond_AHD_CDF_don_chem_type_ss"
p <- ggplot(data=qs) + theme_bw() +
	geom_line(aes(y=probs, x=quantiles, colour=don_chem_type_name)) +
	ggtitle(paste("H-Bond AHD Cumulative Distribution Function by Donor Chemical Type, SeqSep > 5, B-Fact < 30", sep="")) +
	scale_y_continuous("Fraction", limit=c(0,1), breaks=c(0,.25, .50, .75)) +
	scale_x_continuous("Angle Deviation from Linear for Acceptor -- Hydrogen -- Donor (degrees)") +
	geom_indicator(aes(colour=don_chem_type_name, indicator=counts, group=don_chem_type_name), ypos=.7, yjust="bottom") +
	facet_wrap(~sample_source, ncol=ceiling(sqrt(nrow(sample_sources)))) +
	scale_colour_discrete("DonChemType")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

d_ply(qs, .(sample_source), function(sub_qs) {
	ss <- as.character(sub_qs$sample_source[1])
	plot_id <- paste("hbond_AHD_CDF_don_chem_type", ss, sep="_")
	p <- ggplot(data=sub_qs) + theme_bw() +
		geom_line(aes(y=probs, x=quantiles, colour=don_chem_type_name)) +
		geom_indicator(aes(colour=don_chem_type_name, indicator=counts, group=don_chem_type_name), ypos=.7, yjust="bottom") +
		ggtitle(paste("H-Bond AHD Cumulative Distribution Function by Donor Chemical Type SampleSource: ", ss, ", SeqSep > 5, B-Fact < 30", sep="")) +
		scale_y_continuous("Fraction", limit=c(0,1), breaks=c(0,.25, .50, .75)) +
		scale_x_continuous("Angle Deviation from Linear for Acceptor -- Hydrogen -- Donor (degrees)") +
		scale_colour_discrete("DonChemType")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})


qs <- compute_quantiles(
	f, c("sample_source", "acc_chem_type"), "AHD", 1000)
qs$acc_chem_type_name <- acc_chem_type_name_wrap(qs$acc_chem_type)

plot_id = "hbond_AHD_CDF_acc_chem_type"
p <- ggplot(data=qs) + theme_bw() +
	geom_line(aes(y=probs, x=quantiles, colour=sample_source)) +
	ggtitle(paste("H-Bond AHD Cumulative Distribution Function by Acceptor Chemical Type, SeqSep > 5, B-Fact < 30", sep="")) +
	scale_y_continuous("Fraction", limit=c(0,1), breaks=c(0,.25, .50, .75)) +
	scale_x_continuous("Angle Deviation from Linear for Acceptor -- Hydrogen -- Donor (degrees)") +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source), ypos=.7, yjust="bottom") +
	facet_wrap( ~ acc_chem_type_name) +
	scale_colour_discrete("Sample Source")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

qs$acc_chem_type_name <- acc_chem_type_name_linear(qs$acc_chem_type)
plot_id = "hbond_AHD_CDF_acc_chem_type_ss"
p <- ggplot(data=qs) + theme_bw() +
	geom_line(aes(y=probs, x=quantiles, colour=acc_chem_type_name)) +
	ggtitle(paste("H-Bond AHD Cumulative Distribution Function by Acceptor Chemical Type, SeqSep > 5, B-Fact < 30", sep="")) +
	scale_y_continuous("Fraction", limit=c(0,1), breaks=c(0,.25, .50, .75)) +
	scale_x_continuous("Angle Deviation from Linear for Acceptor -- Hydrogen -- Donor (degrees)") +
	geom_indicator(aes(colour=acc_chem_type_name, indicator=counts, group=acc_chem_type_name), ypos=.7, yjust="bottom") +
	facet_wrap( ~ sample_source, ncol=ceiling(sqrt(nrow(sample_sources)))) +
	scale_colour_discrete("AccChemType")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

d_ply(qs, .(sample_source), function(sub_qs) {
	ss <- as.character(sub_qs$sample_source[1])
	plot_id <- paste("hbond_AHD_CDF_acc_chem_type", ss, sep="_")
	p <- ggplot(data=sub_qs) + theme_bw() +
		geom_line(aes(y=probs, x=quantiles, colour=acc_chem_type_name)) +
		geom_indicator(aes(colour=acc_chem_type_name, indicator=counts, group=acc_chem_type_name), ypos=.7, yjust="bottom") +
		ggtitle(paste("H-Bond AHD Cumulative Distribution Function by Acceptor Chemical Type SampleSource: ", ss, ", SeqSep > 5, B-Fact < 30", sep="")) +
		scale_y_continuous("Fraction", limit=c(0,1), breaks=c(0,.25, .50, .75)) +
		scale_x_continuous("Angle Deviation from Linear for Acceptor -- Hydrogen -- Donor (degrees)") +
		scale_colour_discrete("AccChemType")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})


qs <- compute_quantiles(
	f, c("sample_source", "don_chem_type_name", "acc_chem_type_name"), "AHD", 1000)

plot_id = "hbond_AHD_CDF_chem_type"
p <- ggplot(data=qs) + theme_bw() +
	geom_line(aes(y=probs, x=quantiles, colour=sample_source)) +
	ggtitle(paste("H-Bond AHD Cumulative Distribution Function by Chemical Type, SeqSep > 5, B-Fact < 30", sep="")) +
	scale_y_continuous("Fraction", limit=c(0,1), breaks=c(0,.25, .50, .75)) +
	scale_x_continuous("Angle Deviation from Linear for Acceptor -- Hydrogen -- Donor (degrees)") +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source), ypos=.9, yjust="bottom") +
	facet_grid(don_chem_type_name ~ acc_chem_type_name) +
	scale_colour_discrete("Sample Source")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



})) # end FeaturesAnalysis
