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
id = "angle_dependence_on_AHdist",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


source("../hbond_geo_dim_scales.R")

sele <-"
SELECT
  geom.AHdist,
	geom.cosAHD,
	geom.cosBAH,
	don.HBChemType AS don_chem_type, acc.HBChemType AS acc_chem_type
FROM
	hbonds AS hb cross join
	hbond_sites_pdb AS don_pdb cross join
	hbond_sites_pdb AS acc_pdb cross join
	hbond_sites AS don cross join
 	hbond_sites AS acc cross join
	hbond_geom_coords AS geom
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	abs(don.resNum - acc.resNum ) > 5 AND
	geom.AHdist < 4;"

f <- query_sample_sources(sample_sources, sele)


f <- transform(f,
	AHD = acos(f$cosAHD) * 180/pi,
	don_chem_type_name = don_chem_type_name_wrap(don_chem_type),
 	acc_chem_type_name = acc_chem_type_name_wrap(acc_chem_type))
f <- na.omit(f, method="r")

compute_AHD_qs <- function(f, ids, verbose=F){
	w <- sliding_windows(f, ids, "AHdist", verbose=verbose)
	compute_quantiles(w, c(ids, "windows"), "AHD", 1000)
}

plot_parts_AHD <- list(
	theme_bw(),
	geom_line(aes(y=probs, x=quantiles)),
	scale_y_continuous("Fraction", limit=c(0,1), breaks=c(0,.25, .50, .75)),
	scale_x_continuous("Angle Deviation from Linear for Acceptor -- Hydrogen -- Donor (degrees)"),
	geom_indicator(aes(indicator=counts), ypos=.9, yjust="bottom"))


make_id_AHD <- function(ids) {
	paste("hbond_AHD_dependence_on_AHdist_", paste(ids, collapse="_"), sep="")
}
make_title_AHD <- function(ids) {
	ggtitle(
		paste(
			"H-Bond AHD CDF by ",
			paste(ids, collapse=", "),
			" SeqSep > 5, B-Fact < 30, AHDist sliding windows", sep=""))
}

ids <- c("sample_source")
plot_id <- make_id_AHD(ids)
qs <- compute_AHD_qs(f, ids)
p <- ggplot(data=qs, aes(colour=windows, group=windows)) + plot_parts_AHD +
	make_title_AHD(c("")) +
	facet_wrap(~sample_source, ncol=ceiling(sqrt(nrow(sample_sources)))) +
	scale_colour_discrete("AHdist Window")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

ids <- c("sample_source")
plot_id <- paste(make_id_AHD(ids), "_restricted", sep="")
qs <- compute_AHD_qs(f[f$AHdist <= 2.3 & f$AHdist >= 1.9,], ids)
p <- ggplot(data=qs, aes(colour=windows, group=windows)) + plot_parts_AHD +
	make_title_AHD(c("")) +
	facet_wrap(~sample_source, ncol=ceiling(sqrt(nrow(sample_sources)))) +
	scale_colour_discrete("AHdist Window")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


ids <- c("sample_source", "don_chem_type_name")
plot_id <- make_id_AHD(ids)
f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
qs <- compute_AHD_qs(f, ids)
p <- ggplot(data=qs, aes(colour=windows, group=windows)) + plot_parts_AHD +
	make_title_AHD("DonChemType") +
	facet_grid(don_chem_type_name ~ sample_source) +
	scale_colour_discrete("AHdist Window")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

ids <- c("sample_source", "acc_chem_type_name")
plot_id <- make_id_AHD(ids)
f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type)
qs <- compute_AHD_qs(f, ids)
p <- ggplot(data=qs, aes(colour=windows, group=windows)) + plot_parts_AHD +
	make_title_AHD("AccChemType") +
	facet_grid(acc_chem_type_name ~ sample_source) +
	scale_colour_discrete("AHdist Window")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

ids <- c("sample_source", "acc_chem_type_name", "don_chem_type_name")
f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type)
qs <- compute_AHD_qs(f, ids)

d_ply(qs, .(sample_source), function(qsf){
	ss <- as.character(qsf$sample_source[1])
	plot_id <- make_id_AHD(c(ids, ss))
	p <- ggplot(data=qsf, aes(colour=windows, group=windows)) + plot_parts_AHD +
		make_title_AHD(paste("ChemType, SS: ", ss, sep="")) +
		facet_grid(don_chem_type_name ~ acc_chem_type_name) +
		scale_colour_discrete("AHdist Window")
	if(nrow(sample_sources) <= 3){
		p <- p + theme(legend.position="bottom", legend.direction="horizontal")
	}
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})

qs <- compute_AHD_qs(
	f[f$AHdist <= 2.3 & f$AHdist >= 1.9,], "sample_source", verbose=F)
d_ply(qs, c("sample_source"), function(sub_qs){
	ss_id <- as.character(sub_qs[1,"sample_source"])
	plot_id <- paste(make_id_AHD(ss_id), "_restricted", sep="")
	p <- ggplot(data=sub_qs, aes(colour=windows, group=windows)) + plot_parts_AHD +
		make_title_AHD(c("")) +
		scale_colour_discrete("AHdist Window")
	if(nrow(sample_sources) <= 3){
		p <- p + theme(legend.position="bottom", legend.direction="horizontal")
	}
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})


ref_sample_sources <- as.character(sample_sources[sample_sources$reference,"sample_source"])
new_sample_sources <- as.character(sample_sources[!sample_sources$reference,"sample_source"])

qs <- compute_AHD_qs(
	f[f$AHdist <= 2.3 & f$AHdist >= 1.9,], "sample_source", verbose=F)

d_ply(
	qs[qs$sample_source %in% ref_sample_sources,],
	c("sample_source"), function(qs_ref) {
	ref_ss_id <- as.character(qs_ref[1,"sample_source"])
	d_ply(
		qs[qs$sample_source %in% new_sample_sources,],
		c("sample_source"), function(qs_new){
		new_ss_id <- as.character(qs_new[1,"sample_source"])
		plot_id <- paste(make_id_AHD(ref_ss_id), new_ss_id, "_restricted", sep="")
		p <- ggplot(rbind(qs_ref, qs_new)) +
			theme_bw() +
			scale_y_continuous("Fraction", limit=c(0,1), breaks=c(0,.25, .50, .75)) +
			scale_x_continuous("Angle Deviation from Linear for Acceptor -- Hydrogen -- Donor (degrees)") +
			scale_line_type_manual("Sample Source", c(3,1)) +
#			geom_indicator(aes(indicator=counts), ypos=.9, yjust="bottom")) +
			list(make_title_AHD(c(""))) +
			scale_colour_discrete("AHdist Window")
		if(nrow(sample_sources) <= 3){
			p <- p + theme(legend.position="bottom", legend.direction="horizontal")
		}
		save_plots(self, plot_id, sample_sources, output_dir, output_formats)
	})
})





})) # end FeaturesAnalysis
