# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "quarter_chi_chem_type",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
library(ggplot2)


source("../hbond_geo_dim_scales.R")

sele <-"
SELECT
	geom.chi,
	acc.HBChemType AS acc_chem_type,
	don.HBChemType AS don_chem_type,
	CASE acc.HBChemType
		WHEN 'hbacc_IMD' THEN 'ring' WHEN 'hbacc_IME' THEN 'ring'
		WHEN 'hbacc_AHX' THEN 'sp3'  WHEN 'hbacc_HXL' THEN 'sp3'
		WHEN 'hbacc_CXA' THEN 'sp2'  WHEN 'hbacc_CXL' THEN 'sp2'
		WHEN 'hbacc_PBA' THEN 'bb_sp2' END AS hybrid,
	acc_atoms.base2_x AS ab2x, acc_atoms.base2_y AS ab2y, acc_atoms.base2_z AS ab2z, -- acceptor base 2 atom
 	acc_atoms.base_x  AS abx,  acc_atoms.base_y  AS aby,  acc_atoms.base_z  AS abz,  -- acceptor base atom
	acc_atoms.atm_x   AS ax,   acc_atoms.atm_y   AS ay,   acc_atoms.atm_z   AS az,   -- acceptor atom
	don_atoms.atm_x   AS hx,   don_atoms.atm_y   AS hy,   don_atoms.atm_z   AS hz    -- hydrogen atom
FROM
	hbonds AS hb,
	hbond_geom_coords AS geom,
	hbond_sites AS don, hbond_sites AS acc,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb,
	hbond_site_atoms AS don_atoms, hbond_site_atoms AS acc_atoms
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
	acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id AND
	abs(don.resNum - acc.resNum ) > 5;"
f <- query_sample_sources(sample_sources, sele)
f <- f[as.character(f$don_chem_type) != "hbdon_NONE",]

f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type)
f <- na.omit(f, method="r")

alt_chi_dihedral_angle <- function(ab2, ab, a, h){
	alt_ab <- (ab + ab2)/2
	alt_ab2 <- vector_crossprod(ab - ab2, a - ab) - alt_ab
	vector_dihedral(alt_ab2, alt_ab, a, h)
}

f[f$hybrid %in% c("sp3", "ring"), "chi"] <-
	with(f[f$hybrid %in% c("sp3", "ring"),], alt_chi_dihedral_angle(
		cbind(ab2x, ab2y, ab2z), cbind(abx, aby, abz),
		cbind(ax, ay, az), cbind(hx, hy, hz)))

f$all_sp2 <- as.character(f$hybrid) == "sp2" | as.character(f$hybrid) == "bb_sp2"

# angle out of the plane
f$quarter_chi <- asin(abs(sin(f$chi))) * 180/pi

qs <- compute_quantiles(
	f[f$all_sp2,], c("sample_source"), "quarter_chi", 1000)

plot_id = "hbond_quarter_BAchi"
p <- ggplot(data=qs) + theme_bw() +
	geom_line(aes(y=probs, x=quantiles, colour=sample_source)) +
	ggtitle(paste("H-Bond BA-Chi Angle Out of Sp2 plane Cumulative Distribution Function, SeqSep > 5, B-Fact < 30", sep="")) +
	scale_y_continuous("Fraction", limit=c(0,1), breaks=c(0,.25, .50, .75)) +
	scale_x_continuous("BA-Chi Angle Out of the Plane (degrees)") +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source), ypos=.7, yjust="bottom")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


stats <- ddply(f[f$all_sp2,], .(sample_source), function(df){
	data.frame(
		count = length(df$quarter_chi),
		mean = mean(df$quarter_chi),
		median = median(df$quarter_chi),
		sd = sd(df$quarter_chi))
})
table_id <- "hbond_sp2_quarter_BAchi_statistics"
table_title <- "H-Bond BA-Chi Out of Sp2 plane Summary Statistics, SeqSep > 5, B-Factor < 30"
save_tables(self,
	stats, table_id,
	sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")





##############################################################
qs <- compute_quantiles(
	f[f$all_sp2,], c("sample_source", "don_chem_type_name", "acc_chem_type_name"), "quarter_chi", 1000)

plot_id = "hbond_quarter_BAchi_chem_type"
p <- ggplot(data=qs) + theme_bw() +
	geom_line(aes(y=probs, x=quantiles, colour=sample_source)) +
	ggtitle(paste("H-Bond BA-Chi Angle Out of Sp2 plane Cumulative Distribution Function, SeqSep > 5, B-Fact < 30", sep="")) +
	scale_y_continuous("Fraction", limit=c(0,1), breaks=c(0,.25, .50, .75)) +
	scale_x_continuous("BA-Chi Angle Out of the Plane (degrees)") +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source), ypos=.7, yjust="bottom") +
	facet_grid(don_chem_type_name ~ acc_chem_type_name)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


stats <- ddply(f[f$all_sp2,], .(sample_source, don_chem_type_name, acc_chem_type_name), function(df){
	data.frame(
		count = length(df$quarter_chi),
		mean = mean(df$quarter_chi),
		median = median(df$quarter_chi),
		sd = sd(df$quarter_chi))
})
table_id <- "hbond_sp2_quarter_BAchi_chem_type_statistics"
table_title <- "H-Bond BA-Chi Out of Sp2 plane Summary Statistics, SeqSep > 5, B-Factor < 30"
save_tables(self,
	stats, table_id,
	sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")


##############################################################
f$don_chem_type_name <- don_chem_type_name_wrap(f$don_chem_type)

qs <- compute_quantiles(
	f[f$all_sp2,], c("sample_source", "don_chem_type_name"), "quarter_chi", 1000)

plot_id = "hbond_quarter_BAchi_don_chem_type"
p <- ggplot(data=qs) + theme_bw() +
	geom_line(aes(y=probs, x=quantiles, colour=sample_source)) +
	ggtitle(paste("H-Bond BA-Chi Angle Out of Sp2 plane Cumulative Distribution Function, SeqSep > 5, B-Fact < 30", sep="")) +
	scale_y_continuous("Fraction", limit=c(0,1), breaks=c(0,.25, .50, .75)) +
	scale_x_continuous("BA-Chi Angle Out of the Plane (degrees)") +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source), ypos=.7, yjust="bottom") +
	facet_wrap(~don_chem_type_name)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


stats <- ddply(f[f$all_sp2,], .(sample_source, don_chem_type_name), function(df){
	data.frame(
		count = length(df$quarter_chi),
		mean = mean(df$quarter_chi),
		median = median(df$quarter_chi),
		sd = sd(df$quarter_chi))
})
table_id <- "hbond_sp2_quarter_BAchi_don_chem_type_statistics"
table_title <- "H-Bond BA-Chi Out of Sp2 plane Summary Statistics, SeqSep > 5, B-Factor < 30"
save_tables(self,
	stats, table_id,
	sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")



################################################################
f$acc_chem_type_name <- acc_chem_type_name_wrap(f$acc_chem_type)

qs <- compute_quantiles(
	f[f$all_sp2,], c("sample_source", "acc_chem_type_name"), "quarter_chi", 1000)

plot_id = "hbond_quarter_BAchi_acc_chem_type"
p <- ggplot(data=qs) + theme_bw() +
	geom_line(aes(y=probs, x=quantiles, colour=sample_source)) +
	ggtitle(paste("H-Bond BA-Chi Angle Out of Sp2 plane Cumulative Distribution Function, SeqSep > 5, B-Fact < 30", sep="")) +
	scale_y_continuous("Fraction", limit=c(0,1), breaks=c(0,.25, .50, .75)) +
	scale_x_continuous("BA-Chi Angle Out of the Plane (degrees)") +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source), ypos=.7, yjust="bottom") +
	facet_wrap(~acc_chem_type_name)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


stats <- ddply(f[f$all_sp2,], .(sample_source, acc_chem_type_name), function(df){
	data.frame(
		count = length(df$quarter_chi),
		mean = mean(df$quarter_chi),
		median = median(df$quarter_chi),
		sd = sd(df$quarter_chi))
})
table_id <- "hbond_sp2_quarter_BAchi_acc_chem_type_statistics"
table_title <- "H-Bond BA-Chi Out of Sp2 plane Summary Statistics, SeqSep > 5, B-Factor < 30"
save_tables(self,
	stats, table_id,
	sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")








})) # end FeaturesAnalysis
