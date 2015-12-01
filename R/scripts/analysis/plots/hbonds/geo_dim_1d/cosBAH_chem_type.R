# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "cosBAH_chem_type",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
library(ggplot2)


source("../hbond_geo_dim_scales.R")

sele <-"
SELECT
	acc_atoms.base_x  AS bx,  acc_atoms.base_y  AS by,  acc_atoms.base_z  AS bz,
	acc_atoms.base2_x AS b2x, acc_atoms.base2_y AS b2y, acc_atoms.base2_z AS b2z,
	acc_atoms.atm_x   AS ax,  acc_atoms.atm_y   AS ay,  acc_atoms.atm_z   AS az,
	don_atoms.atm_x   AS hx,  don_atoms.atm_y   AS hy,  don_atoms.atm_z   AS hz,
  don.HBChemType AS don_chem_type, acc.HBChemType AS acc_chem_type,
  CASE acc.HBChemType
		WHEN 'hbacc_IMD' THEN 'ring' WHEN 'hbacc_IME' THEN 'ring'
		WHEN 'hbacc_AHX' THEN 'sp3'  WHEN 'hbacc_HXL' THEN 'sp3'
		WHEN 'hbacc_CXA' THEN 'sp2'  WHEN 'hbacc_CXL' THEN 'sp2'
		WHEN 'hbacc_PBA' THEN 'bb_sp2' END AS hybrid
FROM
  hbonds AS hb,
  hbond_geom_coords AS geom,
  hbond_sites AS don, hbond_sites AS acc,
	hbond_site_atoms AS don_atoms, hbond_site_atoms AS acc_atoms,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
	acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	abs(don.resNum - acc.resNum ) > 5;"
f <- query_sample_sources(sample_sources, sele)

f <- f[f$don_chem_type != "hbdon_NONE" & f$acc_chem_type != "hbacc_NONE",]
f <- with(f, data.frame(
	sample_source = sample_source,
	cosBAH = ifelse(hybrid %in% c("sp3", "ring"),
		vector_dotprod(
		vector_normalize(cbind(ax-(bx+b2x)/2, ay-(by+b2y)/2, az-(bz+b2z)/2)),
		vector_normalize(cbind(hx-ax, hy-ay, hz-az))),
		vector_dotprod(
		vector_normalize(cbind(ax-bx, ay-by, az-bz)),
		vector_normalize(cbind(hx-ax, hy-ay, hz-az)))),
	don_chem_type = don_chem_type,
	acc_chem_type = acc_chem_type,
	don_chem_type_name = don_chem_type_name_linear(don_chem_type),
	acc_chem_type_name = acc_chem_type_name_linear(acc_chem_type),
	hybrid = factor(
		hybrid,
		levels=c("ring", "sp3", "sp2", "bb_sp2"),
		labels=c("ring", "Sp3", "Sp2", "BB Sp2"))))
f <- na.omit(f, method="r")


plot_parts <- list(
	theme_bw(),
	geom_line(aes(x=acos(x)*180/pi, y=y)),
	geom_indicator(aes(indicator=counts)),
	scale_x_continuous(paste('Base -- Acceptor -- Hydrogen (degrees)')),
	scale_y_continuous("Feature Density", breaks=c(0, 2, 4)),
	coord_trans(limx=c(0, 90), limy=c(0,6)))

compute_dens <- function(ids){
	dens <- estimate_density_1d_reflect_boundary(
		data = f,
		ids = ids,
		variable = "cosBAH",
		reflect_left=T, left_boundary=0,
		reflect_right=T, right_boundary=pi, adjust=.2)
}

make_id <- function(id) {
	paste("hbond_BAH_", id, sep="")
}

make_title <- function(cond) {
	ggtitle(paste("H-Bond BAH Angle by ", cond, " SeqSep > 5, B-Factor < 30\n(normalized for equal volume per unit distance)", sep=""))
}
################################################################################

dens <- compute_dens(c("sample_source", "hybrid"))

plot_id <- make_id("acc_hybrid")
p <- ggplot(data=dens, aes(color=sample_source, group=sample_source)) + plot_parts +
	facet_wrap(~ hybrid) +
	make_title("Acceptor Hybrid") +
	scale_colour_discrete("Sample Source")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

d_ply(dens, .(sample_source), function(sub_dens) {
	ss <- as.character(sub_dens$sample_source[1])
	plot_id <- make_id(paste("acc_hybrid", ss, sep="_"))
	p <- ggplot(data=sub_dens, aes(color=hybrid, group=hybrid)) + plot_parts +
		make_title(paste("Accceptor Hybrid; SampleSource: ", ss, sep="")) +
		scale_colour_discrete("AccHybrid")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})

################################################################################

dens <- compute_dens(c("sample_source", "hybrid", "don_chem_type_name"))

plot_id <- make_id("don_chem_type_acc_hybrid")
p <- ggplot(data=dens, aes(color=sample_source, group=sample_source)) + plot_parts +
	facet_grid(don_chem_type_name ~ hybrid) +
	make_title("Donor Chemical Type, Acceptor Hybrid") +
	scale_colour_discrete("Sample Source")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- make_id("don_chem_type_acc_hybrid_ss")
p <- ggplot(data=dens, aes(color=don_chem_type_name, group=don_chem_type_name)) + plot_parts +
	facet_grid(sample_source ~ hybrid) +
	make_title("Donor Chemical Type, Acceptor Hybrid") +
	scale_colour_discrete("DonChemType")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

d_ply(dens, .(hybrid), function(sub_dens){
	hybrid <- as.character(sub_dens$hybrid[1])
	plot_id <- make_id(paste("don_chem_type_acc_hybrid", hybrid, sep="_"))
	p <- ggplot(data=sub_dens, aes(color=sample_source, group=sample_source)) + plot_parts +
		facet_wrap(~ don_chem_type_name) +
		make_title(paste("Donor Chemical Type; Acceptor Hybrid: ", hybrid, sep="")) +
		scale_colour_discrete("Sample Source")
	if(nrow(sample_sources) <= 3){
		p <- p + theme(legend.position="bottom", legend.direction="horizontal")
	}
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})



####################################

dens <- compute_dens(c("sample_source", "acc_chem_type", "hybrid"))
dens$acc_chem_type_name <- acc_chem_type_name_wrap(dens$acc_chem_type)

plot_id <- make_id("acc_chem_type")
p <- ggplot(data=dens, aes(colour=sample_source, group=sample_source)) + plot_parts +
	facet_wrap( ~ acc_chem_type_name) +
	make_title("Acceptor Chemical Type") +
	scale_colour_discrete("Sample Source")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

dens$acc_chem_type_name <- acc_chem_type_name_linear(dens$acc_chem_type)
plot_id <- make_id("acc_chem_type_acc_hybrid_ss")
p <- ggplot(data=dens, aes(colour=acc_chem_type_name, group=acc_chem_type_name)) + plot_parts +
	facet_grid( hybrid ~ sample_source ) +
	make_title("Acceptor Chemical Type By Acceptor Hybrid") +
	scale_colour_discrete("AccChemType")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

d_ply(dens, .(sample_source), function(sub_dens) {
	ss <- as.character(sub_dens$sample_source[1])
	plot_id <- make_id(paste("acc_chem_type", ss, sep="_"))
	p <- ggplot(data=sub_dens, aes(color=acc_chem_type_name, group=acc_chem_type_name)) + plot_parts +
		geom_line(aes(x=acos(x)*180/pi, y=y, colour=acc_chem_type_name)) +
		make_title(paste("Acceptor Chemical Type; Sample Source: ", ss, sep="")) +
		scale_colour_discrete("AccChemType")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})




###########################################################

dens <- compute_dens(c("sample_source", "acc_chem_type", "don_chem_type"))
dens$don_chem_type_name <- don_chem_type_name_linear(dens$don_chem_type)
dens$acc_chem_type_name <- acc_chem_type_name_linear(dens$acc_chem_type)

plot_id <- make_id("chem_type")
p <- ggplot(data=dens, aes(colour=sample_source, group=sample_source)) + plot_parts +
	facet_grid(don_chem_type_name ~ acc_chem_type_name) +
	make_title("Chemical Type") +
	scale_colour_discrete("Sample Source")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

########################################################################

f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type)


compute_stats <- function(ids, table_id, table_title){
	stats <- estimate_primary_modes_1d(f, ids, "cosBAH")
	stats$primary_mode <- round(acos(stats$primary_mode) * 180/pi, 2)
	save_tables(
		self,
		stats, make_id(paste("primary_mode_", table_id, sep="")),
		sample_sources, output_dir, output_formats,
		caption=paste("H-Bond BAH Angle Primary Mode By", table_title, "; SeqSep > 5, B-Fact < 30", sep=""),
		caption.placement="top")

}

stats <- compute_stats(
	ids=c("sample_source", "acc_chem_type_name", "don_chem_type_name"),
	table_id = "chem_type",
	table_title = "Chem Type")

stats <- compute_stats(
	ids=c("sample_source", "acc_chem_type_name"),
	table_id = "acc_chem_type",
	table_title = "Acceptor Chem Type")

stats <- compute_stats(
	ids=c("sample_source", "don_chem_type_name"),
	table_id = "don_chem_type",
	table_title = "Donor Chem Type")

stats <- compute_stats(
	ids=c("sample_source", "hybrid"),
	table_id = "acc_hybrid",
	table_title = "Acceptor Hybrid")

stats <- compute_stats(
	ids=c("sample_source", "hybrid", "don_chem_type_name"),
	table_id = "acc_hybrid_don_chem_type",
	table_title = "Acceptor Hybrid Donor Chem Type")

})) # end FeaturesAnalysis
