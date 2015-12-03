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

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "Sp2_summary",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){



sele <-"
SELECT
	geom.chi, geom.cosBAH,
	acc.HBChemType AS acc_chem_type,
	don.HBChemType AS don_chem_type,
	acc_ss_code.label AS acc_ss,
	abs(don.resNum - acc.resNum) AS seq_sep
FROM
	hbonds AS hb cross join
	hbond_sites AS acc cross join
	hbond_sites AS don cross join
	hbond_sites_pdb AS don_pdb cross join
	hbond_sites_pdb AS acc_pdb cross join
	hbond_geom_coords AS geom cross join
	residue_secondary_structure AS acc_ss cross join
	dssp_codes AS acc_ss_code
WHERE
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	(acc.HBChemType == 'hbacc_PBA' OR
		acc.HBChemType == 'hbacc_CXL' OR
		acc.HBChemType == 'hbacc_CXA') AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	acc_ss.struct_id = acc.struct_id AND acc_ss.resNum = acc.resNum AND
	acc_ss_code.code = acc_ss.dssp AND
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id;"

f <- query_sample_sources(sample_sources, sele)
f <- f[as.character(f$don_chem_type) != "hbdon_NONE",]

f <- transform(f,
	don_chem_type_name_linear = don_chem_type_name_linear(don_chem_type),
 	don_chem_type_name_wrap = don_chem_type_name_wrap(don_chem_type),
	acc_chem_type_name = factor(
  	acc_chem_type,
		levels = c("hbacc_PBA", "hbacc_CXA", "hbacc_CXL"),
		labels = c("aPBA: bb", "aCXA: n,q", "aCXL: d,e")),
	acc_ss_wrap = factor(
		acc_ss,
		levels = c(
			"H: a-Helix", "E: b-Sheet", "T: HB Turn",
			"G: 3/10 Helix", "B: b-Bridge", "S: Bend",
			"I: pi-Helix", "Irregular"),
		labels = c(
			"H: a-Helix", "E: b-Sheet", "T: HB Turn",
			"G: 3/10 Helix", "B: b-Bridge", "S: Bend",
			"I: pi-Helix", "Irregular")),
	acc_ss_linear = factor(
		acc_ss,
		levels = c(
			"H: a-Helix", "G: 3/10 Helix", "I: pi-Helix",
			"E: b-Sheet", "B: b-Bridge",
			"T: HB Turn", "S: Bend", "Irregular"),
		labels = c(
			"H: a-Helix", "G: 3/10 Helix", "I: pi-Helix",
			"E: b-Sheet", "B: b-Bridge",
			"T: HB Turn", "S: Bend", "Irregular")),
	chi = (((chi*180/pi) + 90) %% 360) - 90,
	quarter_chi = asin(abs(sin(f$chi))) * 180/pi)
f <- na.omit(f, method="r")




generic_plot_parts <- list(
	theme_bw(),
	geom_indicator(aes(indicator=counts)))

make_plots <- function( ss, dim, dim_name, plot_parts, compute_line ) {
	dens <- compute_line(c("acc_chem_type_name"))
	plot_id <- make_id(dim, "acc_chem_type", ss)
	p <- ggplot(
		data=dens,
		aes(color=acc_chem_type_name, group=acc_chem_type_name)) +
		plot_parts +
		make_title(dim_name, "Acceptor Chemical Type") +
		scale_colour_discrete("AccChem") +
		theme(legend.position=c(.85,.45)) +
		theme(legend.justification=c("left", "top"))
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	dens <- compute_line(
		c("acc_chem_type_name",
			"don_chem_type_name_linear",
			"don_chem_type_name_wrap"))
	plot_id <- make_id(dim, "acc_chem_type_don_chem_type", ss)
	p <- ggplot(
		data=dens,
		aes(color=don_chem_type_name_linear, group=don_chem_type_name_linear)) +
		facet_wrap( ~ acc_chem_type_name, ncol=2) +
		plot_parts +
		make_title(dim_name, c("AccChem Type", "DonChem Type")) +
		scale_colour_discrete("DonChem") +
		theme(legend.position=c(.85,.45)) +
		theme(legend.justification=c("left", "top"))
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	plot_id <- make_id(dim, "don_chem_type_acc_chem_type", ss)
	p <- ggplot(
		data=dens,
		aes(color=acc_chem_type_name, group=acc_chem_type_name)) +
		facet_wrap( ~ don_chem_type_name_wrap) +
		plot_parts +
		make_title(dim_name, c("DonChem Type", "AccChem Type")) +
		scale_colour_discrete("AccChem")  +
		theme(legend.position=c(.58,.35)) +
		theme(legend.justification=c("left", "top"))
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	dens <- compute_line(c("acc_chem_type_name", "acc_ss_wrap", "acc_ss_linear"))
	plot_id <- make_id(dim, "acc_chem_type_acc_dssp", ss)
	p <- ggplot(
		data=dens,
		aes(color=acc_ss_linear, group=acc_ss_linear)) +
		facet_wrap( ~ acc_chem_type_name, ncol=2) +
		plot_parts +
		make_title(dim_name, c("AccChem Type", "Acc DSSP")) +
		scale_colour_discrete("AccDSSP") +
		theme(legend.position=c(.85,.45)) +
		theme(legend.justification=c("left", "top"))
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	plot_id <- make_id(dim, "acc_dssp_acc_chem_type", ss)
	p <- ggplot(
		data=dens,
		aes(color=acc_chem_type_name, group=acc_chem_type_name)) +
		facet_wrap( ~ acc_ss_wrap) +
		plot_parts +
		make_title(dim_name, c("Acc DSSP", "AccChem Type")) +
		scale_colour_discrete("AccChem") +
		theme(legend.position=c(.75,.35)) +
		theme(legend.justification=c("left", "top"))
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
}


make_all_plots <- function(zf, make_title, make_id) {
	d_ply(zf, .(sample_source), function(df) {
		ss <- as.character(df$sample_source[1])

		########################### BAH ##########################
		plot_parts <- c(generic_plot_parts, list(
			geom_line(aes(x=acos(x)*180/pi, y=y)),
			scale_x_continuous(
				expression(paste("BAH (", degree, ")"))),
			scale_y_continuous("Feature Density", breaks=c(0, 2, 4)),
			coord_trans(limx=c(0, 90), limy=c(0,6))))

		compute_line <- function(ids){
			estimate_density_1d_reflect_boundary(
				data = df,
				ids = ids,
				variable = "cosBAH",
				reflect_left=T, left_boundary=0,
				reflect_right=T, right_boundary=pi, adjust=.2)
		}

		make_plots(ss, "cosBAH", "BAH", plot_parts, compute_line)

		########## BA-Chi ##########
		plot_parts <- c(generic_plot_parts, list(
			geom_line(aes(x=x, y=y*1000)),
			scale_x_continuous(
				expression(paste("BA", chi, " (", degree, ")")),
				breaks=c(0, 180)),
			scale_y_continuous('Feature Density', breaks=c(0,2,4,6,8)),
			coord_trans(limx=c(-90, 270), limy=c(0,9))))

		compute_line <- function(ids) {
			estimate_density_1d_wrap(df, ids, "chi", xlim=c(-90, 270))
		}

		make_plots(ss, "chi", "BA-Chi", plot_parts, compute_line)


		########################### quarter_chi ##########################
		plot_parts <- list(
			theme_bw(),
			geom_line(aes(x=quantiles, y=probs)),
			scale_y_continuous("Fraction", limit=c(0,1), breaks=c(0,.25, .50, .75)),
			scale_x_continuous(
				expression(paste(frac("1", "4"), "BA", chi, " (", degree, ")"))),
			geom_indicator(aes(indicator=counts), ypos=.7, yjust="bottom"))

		compute_line <- function(ids){
			compute_quantiles(df, ids, "quarter_chi")
		}

		make_plots(ss, "quarter_chi", "Quarter BA-Chi", plot_parts, compute_line)

	})
}


################# All the data #######################

make_title <- function(dim, cond=c()) {
	ggtitle(paste(
		"H-Bond ", dim, "; ",
		ifelse(
			length(cond),
			paste(paste(cond, collapse=", "), "; ", sep=""),
			""),
		"B-Factor < 30", sep=""))
}

make_id <- function(dim, id, ss) {
	paste( "hbond_Sp2", dim, id, ss, sep="_", collapse="_")
}

make_all_plots(f, make_title, make_id)


################ Seq sep < 5 #########################

make_title <- function(dim, cond=c()) {
	ggtitle(paste(
		"H-Bond ", dim, "; ",
		ifelse(
			length(cond),
			paste(paste(cond, collapse=", "), "; ", sep=""),
			""),
		"SeqSep < 5, B-Factor < 30", sep=""))
}

make_id <- function(dim, id, ss) {
	paste( "hbond_Sp2_short", dim, id, ss, sep="_", collapse="_")
}

make_all_plots(f[f$seq_sep < 5,], make_title, make_id)

################ Seq sep >= 5 #########################

make_title <- function(dim, cond=c()) {
	ggtitle(paste(
		"HBond ", dim, "; ",
		ifelse(
			length(cond),
			paste(paste(cond, collapse=", "), "; ", sep=""),
			""),
		"SeqSep >= 5, B-Factor < 30", sep=""))
}

make_id <- function(dim, id, ss) {
	paste( "hbond_Sp2_long", dim, id, ss, sep="_", collapse="_")
}

make_all_plots(f[f$seq_sep >= 5,], make_title, make_id)



})) # end FeaturesAnalysis
