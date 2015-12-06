# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

library(reshape2)
library(ggplot2)
source("../../plots/hbonds/hbond_geo_dim_scales.R")

feature_analyses <- c(feature_analyses, methods::new("FeaturesAnalysis",
id = "AHdist_min_chem_type_with_params",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
	geom.AHdist,
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
	acc_chem_type_name = acc_chem_type_name_linear(acc_chem_type))

ref_ss <- sample_sources[sample_sources$reference, "sample_source"]
if(length(ref_ss) != 1) {
	stop("ERROR: This analysis script requires a single reference sample source")
}
new_ss <- sample_sources[!sample_sources$reference,"sample_source"]



modes.don <- estimate_primary_modes_1d(f[f$acc_chem_type != "hbacc_PBA",], c("sample_source", "don_chem_type_name"), "AHdist")
names(modes.don)[2] <- "chem_type_name"
modes.acc <- estimate_primary_modes_1d(f[f$don_chem_type != "hbdon_PBA",], c("sample_source", "acc_chem_type_name"), "AHdist")
names(modes.acc)[2] <- "chem_type_name"

modes <- rbind(modes.don, modes.acc)

t <- dcast(modes, chem_type_name ~ sample_source, value="primary_mode")

table_id <- "AHdist_primary_mode_by_don_or_acc_chem_type"
table_title <- "Primary Mode of H-Bond A-H Length"
save_tables(self,
	t, table_id,
	sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")

###########################

modes.ref <- modes[modes$sample_source == ref_ss,]
names(modes.ref)[1] <- "ref_sample_source"
names(modes.ref)[3] <- "ref_primary_mode"
modes.new <- modes[modes$sample_source %in% new_ss,]
names(modes.new)[1] <- "new_sample_source"
names(modes.new)[3] <- "new_primary_mode"
modes <- merge(modes.new, modes.ref)

plot_id <- "AHdist_primary_mode_by_don_or_acc_chem_type"
p <- ggplot(data=modes) + theme_bw() +
	geom_abline(slope=1) +
	geom_point(aes(x=ref_primary_mode, y=new_primary_mode, colour=new_sample_source)) +
	coord_equal(ratio=1) +
	scale_x_continuous(
		paste("A-H Length (A), Ref: ", ref_ss, sep="")) +
	scale_y_continuous(paste("A-H Length (A) Candidate", sep="")) +
	scale_colour_discrete("") +
	ggtitle("A-H Lengths: Reference vs Candidate Sample Source\nGrouped by Chemical Types")

if(length(new_ss) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

alt_output_formats <- transform(output_formats, width=height*.65)
save_plots(self, plot_id, sample_sources, output_dir, alt_output_formats)

plot_id <- "AHdist_primary_mode_by_don_or_acc_chem_type_text"
p <- ggplot(data=modes) + theme_bw() +
	geom_abline(slope=1) +
	geom_text(aes(x=ref_primary_mode, y=new_primary_mode, colour=new_sample_source, label=chem_type_name), size=2) +
#	geom_point(aes(x=ref_primary_mode, y=new_primary_mode, colour=new_sample_source)) +
	coord_equal(ratio=1) +
	scale_x_continuous(
		paste("A-H Length (A), Ref: ", ref_ss, sep="")) +
	scale_y_continuous(paste("A-H Length (A) Candidate", sep="")) +
	scale_colour_discrete("") +
	ggtitle("A-H Lengths: Reference vs Candidate Sample Source\nGrouped by Chemical Types, SeqSep > 4, SC Partner")

if(length(new_ss) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

alt_output_formats <- transform(output_formats, width=height*.65)
save_plots(self, plot_id, sample_sources, output_dir, alt_output_formats)

######################################

modes.all <- estimate_primary_modes_1d(f, c("sample_source", "don_chem_type_name", "acc_chem_type_name"), "AHdist")

t <- dcast(modes.all, don_chem_type_name + acc_chem_type_name ~ sample_source, value="primary_mode")

ref_ss <- sample_sources[sample_sources$reference, "sample_source"]
if(length(ref_ss) != 1) {
	stop("ERROR: This analysis script requires a single reference sample source")
}
new_ss <- sample_sources[!sample_sources$reference,"sample_source"]

modes.all.ref <- modes.all[modes.all$sample_source == ref_ss,]
names(modes.all.ref)[1] <- "ref_sample_source"
names(modes.all.ref)[4] <- "ref_primary_mode"
modes.all.new <- modes.all[modes.all$sample_source %in% new_ss,]
names(modes.all.new)[1] <- "new_sample_source"
names(modes.all.new)[4] <- "new_primary_mode"
modes.all <- merge(modes.all.new, modes.all.ref)

plot_id <- "AHdist_primary_mode_by_don_acc_chem_type"
p <- ggplot(data=modes.all) + theme_bw() +
	geom_abline(slope=1) +
	geom_point(aes(x=ref_primary_mode, y=new_primary_mode, colour=new_sample_source)) +
	coord_equal(ratio=1) +
	scale_x_continuous(
		paste("A-H Length (A), Ref: ", ref_ss, sep="")) +
	scale_y_continuous(paste("A-H Length (A) Candidate", sep="")) +
	scale_colour_discrete("") +
	ggtitle("A-H Lengths: Reference vs Candidate Sample Source\nGrouped by Donor and Acceptor Chemical Types")

if(length(new_ss) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

alt_output_formats <- transform(output_formats, width=height*.8)
save_plots(self, plot_id, sample_sources, output_dir, alt_output_formats)


plot_id <- "AHdist_primary_mode_by_don_acc_chem_type_text"
p <- ggplot(data=modes.all) + theme_bw() +
	geom_abline(slope=1) +
		geom_text(aes(x=ref_primary_mode, y=new_primary_mode, colour=new_sample_source, label=interaction(don_chem_type_name, acc_chem_type_name, sep="\n")), size=2) +
	coord_equal(ratio=1) +
	scale_x_continuous(
		paste("A-H Length (A), Ref: ", ref_ss, sep="")) +
	scale_y_continuous(paste("A-H Length (A) Candidate", sep="")) +
	scale_colour_discrete("") +
	ggtitle("A-H Lengths: Reference vs Candidate Sample Source\nGrouped by Donor and Acceptor Chemical Types")

if(length(new_ss) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

alt_output_formats <- transform(output_formats, width=height*.8)
save_plots(self, plot_id, sample_sources, output_dir, alt_output_formats)

})) # end FeaturesAnalysis

