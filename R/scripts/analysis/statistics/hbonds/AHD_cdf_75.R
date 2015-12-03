# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

library(ggplot2)


feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "AHD_cdf_75",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){




source("../../plots/hbonds/hbond_geo_dim_scales.R")

sele <-"
SELECT
	geom.cosAHD,
	acc.HBChemType AS acc_chem_type,
	don.HBChemType AS don_chem_type
FROM
	hbonds AS hb,
	hbond_geom_coords AS geom,
	hbond_sites AS don,
	hbond_sites AS acc,
	hbond_sites_pdb AS don_pdb,
	hbond_sites_pdb AS acc_pdb
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	ABS(don.resNum - acc.resNum) > 4;"

f <- query_sample_sources(sample_sources, sele)

f <- transform(f,
	don_chem_type_name = don_chem_type_name_linear(don_chem_type),
	acc_chem_type_name = acc_chem_type_name_linear(acc_chem_type),
	AHD=acos(cosAHD))

f <- na.omit(f, method="r")

ref_ss <- sample_sources[sample_sources$reference, "sample_source"]
if(length(ref_ss) != 1) {
	stop("ERROR: This analysis script requires a single reference sample source")
}
new_ss <- sample_sources[!sample_sources$reference,"sample_source"]



cdf.don <- compute_quantiles(f[f$acc_chem_type != "hbacc_PBA",], c("sample_source", "don_chem_type_name"), "AHD", .75)
names(cdf.don)[2] <- "chem_type_name"
cdf.acc <- compute_quantiles(f[f$don_chem_type != "hbdon_PBA",], c("sample_source", "acc_chem_type_name"), "AHD", .75)
names(cdf.acc)[2] <- "chem_type_name"

cdf.chem <- rbind(cdf.don, cdf.acc)
cdf.chem$quantiles <- cdf.chem$quantiles * 180/pi

t <- cast(cdf.chem, chem_type_name + probs ~ sample_source, value="quantiles")

table_id <- "AHD_cdf_don_or_acc_chem_type"
table_title <- "A-H-D Angle containing 75% of H-Bonds (Degrees)\nB-Factor < 30, SeqSep > 4, SC-Partner"
save_tables(self,
	t, table_id,
	sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")



##########################

cdf.ref <- cdf.chem[cdf.chem$sample_source == ref_ss,c("sample_source", "chem_type_name", "quantiles")]
names(cdf.ref)[1] <- "ref_sample_source"
names(cdf.ref)[3] <- "ref_quantile"
cdf.new <- cdf.chem[cdf.chem$sample_source %in% new_ss,c("sample_source", "chem_type_name", "quantiles")]
names(cdf.new)[1] <- "new_sample_source"
names(cdf.new)[3] <- "new_quantile"
cdf.chem_ref_new <- merge(cdf.new, cdf.ref)

plot_id <- "AHD_cdf_don_or_acc_chem_type_qq"
p <- ggplot(data=cdf.chem_ref_new) + theme_bw() +
	geom_abline(slope=1) +
	geom_point(aes(x=ref_quantile, y=new_quantile, colour=new_sample_source)) +
	coord_equal(ratio=1) +
	scale_x_continuous(
		paste("A-H-D Angle (Degrees), Ref: ", ref_ss, sep="")) +
	scale_y_continuous(paste("A-H-D Angle Candidate", sep="")) +
	scale_colour_discrete("") +
	ggtitle("A-H-D Angle Containing 75% of H-Bonds\n\nB-Factor < 30, SeqSep > 4, SC-Partner")

if(length(new_ss) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

alt_output_formats <- transform(output_formats, width=height*.65)
save_plots(self, plot_id, sample_sources, output_dir, alt_output_formats)


plot_id <- "AHD_cdf_don_or_acc_chem_type_qq_text"
p <- ggplot(data=cdf.chem_ref_new) + theme_bw() +
	geom_abline(slope=1) +
	geom_text(aes(x=ref_quantile, y=new_quantile, colour=new_sample_source, label=chem_type_name), size=2) +
#	coord_equal(ratio=1) +
	scale_x_continuous(
		paste("A-H-D Angle (Degrees), Ref: ", ref_ss, sep="")) +
	scale_y_continuous(paste("A-H-D Angle Candidate", sep="")) +
	scale_colour_discrete("") +
	ggtitle("A-H-D Angle Containing 75% of H-Bonds\n\nB-Factor < 30, SeqSep > 4, SC-Partner")

if(length(new_ss) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

#alt_output_formats <- transform(output_formats, width=height*.65)
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

#
#######################################
#
#modes.all <- estimate_primary_modes_1d(f, c("sample_source", "don_chem_type_name", "acc_chem_type_name"), "AHdist")
#
#t <- cast(modes.all, don_chem_type_name + acc_chem_type_name ~ sample_source, value="primary_mode")
#
#ref_ss <- sample_sources[sample_sources$reference, "sample_source"]
#if(length(ref_ss) != 1) {
#	stop("ERROR: This analysis script requires a single reference sample source")
#}
#new_ss <- sample_sources[!sample_sources$reference,"sample_source"]
#
#modes.all.ref <- modes.all[modes.all$sample_source == ref_ss,]
#names(modes.all.ref)[1] <- "ref_sample_source"
#names(modes.all.ref)[4] <- "ref_primary_mode"
#modes.all.new <- modes.all[modes.all$sample_source %in% new_ss,]
#names(modes.all.new)[1] <- "new_sample_source"
#names(modes.all.new)[4] <- "new_primary_mode"
#modes.all <- merge(modes.all.new, modes.all.ref)
#
#plot_id <- "AHdist_primary_mode_by_don_acc_chem_type"
#p <- ggplot(data=modes.all) + theme_bw() +
#	geom_abline(slope=1) +
#	geom_point(aes(x=ref_primary_mode, y=new_primary_mode, colour=new_sample_source)) +
#	coord_equal(ratio=1) +
#	scale_x_continuous(
#		paste("A-H Length (A), Ref: ", ref_ss, sep="")) +
#	scale_y_continuous(paste("A-H Length (A) Candidate", sep="")) +
#	scale_colour_discrete("") +
#	ggtitle("A-H Lengths: Reference vs Candidate Sample Source\nGrouped by Donor and Acceptor Chemical Types")
#
#if(length(new_ss) <= 3){
#	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
#}
#
#alt_output_formats <- transform(output_formats, width=height*.8)
#save_plots(self, plot_id, sample_sources, output_dir, alt_output_formats)
#
#
#plot_id <- "AHdist_primary_mode_by_don_acc_chem_type_text"
#p <- ggplot(data=modes.all) + theme_bw() +
#	geom_abline(slope=1) +
#		geom_text(aes(x=ref_primary_mode, y=new_primary_mode, colour=new_sample_source, label=interaction(don_chem_type_name, acc_chem_type_name, sep="\n")), size=2) +
#	coord_equal(ratio=1) +
#	scale_x_continuous(
#		paste("A-H Length (A), Ref: ", ref_ss, sep="")) +
#	scale_y_continuous(paste("A-H Length (A) Candidate", sep="")) +
#	scale_colour_discrete("") +
#	ggtitle("A-H Lengths: Reference vs Candidate Sample Source\nGrouped by Donor and Acceptor Chemical Types")
#
#if(length(new_ss) <= 3){
#	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
#}
#
#alt_output_formats <- transform(output_formats, width=height*.8)
#save_plots(self, plot_id, sample_sources, output_dir, alt_output_formats)

})) # end FeaturesAnalysis

