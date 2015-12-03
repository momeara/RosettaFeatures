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

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "VBAH_chem_type",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){



sele <-"
SELECT
	acc_atoms.base_x  AS abx,  acc_atoms.base_y  AS aby,  acc_atoms.base_z  AS abz,
	acc_atoms.base2_x AS ab2x, acc_atoms.base2_y AS ab2y, acc_atoms.base2_z AS ab2z,
	acc_atoms.atm_x   AS ax,  acc_atoms.atm_y   AS ay,  acc_atoms.atm_z   AS az,
	don_atoms.atm_x   AS hx,  don_atoms.atm_y   AS hy,  don_atoms.atm_z   AS hz,
	acc.HBChemType AS acc_chem_type,
	don.HBChemType AS don_chem_type
FROM
	hbonds AS hb cross join
	hbond_sites AS don cross join
	hbond_sites AS acc cross join
	hbond_sites_pdb AS don_pdb cross join
	hbond_sites_pdb AS acc_pdb cross join
	hbond_site_atoms AS don_atoms cross join
	hbond_site_atoms AS acc_atoms
WHERE
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	(acc.HBChemType = 'hbacc_HXL' OR acc.HBChemType = 'hbacc_AHX') AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
	acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id;"


f <- query_sample_sources(sample_sources, sele)

f <- with(f, data.frame(
	sample_source = sample_source,
	don_chem_type = don_chem_type,
	acc_chem_type = acc_chem_type,
	don_chem_type_name = factor(don_chem_type_name_linear(f$don_chem_type)),
	acc_chem_type_name = factor(acc_chem_type,
		levels = c("hbacc_HXL", "hbacc_AHX"),
		labels = c("aHXL: s,t", "aAHX: y")),
	cosBAH = vector_dotprod(
		vector_normalize(cbind(ax-(abx+ab2x)/2, ay-(aby+ab2y)/2, az-(abz+ab2z)/2)),
		vector_normalize(cbind(hx-ax, hy-ay, hz-az)))))
f <- na.omit(f, method="r")


f <- na.omit(f, method="r")

dens <- estimate_density_1d_reflect_boundary(
  data = f,
  ids = c("sample_source", "acc_chem_type_name", "don_chem_type_name"),
  variable = "cosBAH",
  reflect_left=T, left_boundary=0,
  reflect_right=T, right_boundary=pi, adjust=.2
	)

stats <- estimate_primary_modes_1d(
	f, c("sample_source", "don_chem_type_name", "acc_chem_type_name"), "cosBAH")
stats$primary_mode <- acos(stats$primary_mode) * 180/pi

plot_id <- "hbond_vBAH_hydroxyl_acceptor_chem_type"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=acos(x)*180/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source)) +
	geom_indicator(
		data=stats,
		aes(colour=sample_source, indicator=primary_mode, group=sample_source),
		xpos="left") +
	facet_grid(don_chem_type_name ~ acc_chem_type_name) +
	ggtitle("HBond BAH Angle by Chemical Type, SeqSep > 5, B-Factor < 30\n(normalized for equal volume per unit distance)") +
	scale_x_continuous(paste('Base -- Acceptor -- Hydrogen (degrees)')) +
	scale_y_continuous("Feature Density", limits=c(0,5), breaks=c(0, 2, 4))
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

table_id <- "hbond_vBAH_hydroxyl_acceptor_primary_mode_chem_type"
table_title <- "H-Bond BAH Angle Primary Mode By Chem Type"
save_tables(self,
	stats, table_id,
	sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")



########################################################################

f <- transform(f,
	don_chem_type_name = don_chem_type_name_wrap(f$don_chem_type))

dens <- estimate_density_1d_reflect_boundary(
  data = f,
  ids = c("sample_source", "acc_chem_type_name"),
  variable = "cosBAH",
  reflect_left=T, left_boundary=0,
  reflect_right=T, right_boundary=pi, adjust=.2)

stats <- estimate_primary_modes_1d(
	f, c("sample_source", "acc_chem_type_name"), "cosBAH")
stats$BAH_primary_mode <- acos(stats$primary_mode) * 180/pi

table_id <- "hbond_vBAH_hydroxyl_acceptor_primary_mode_acc_chem_type"
table_title <- "H-Bond BAH Angle Primary Mode By Acceptor Chem Type"
save_tables(self,
	stats, table_id,
	sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")

plot_id <- "hbond_vBAH_hydroxyl_acceptor_acc_chem_type"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=acos(x)*180/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source)) +
	geom_indicator(
		data=stats,
		aes(colour=sample_source, indicator=primary_mode, group=sample_source),
		xpos="left") +
	facet_wrap( ~ acc_chem_type_name) +
	ggtitle("HBond Sp3 Acceptor Virtual-BAH Angle by Acceptor Chemical Type, SeqSep > 5, B-Factor < 30\n(normalized for equal volume per unit distance)") +
	scale_x_continuous(paste('Base -- Acceptor -- Hydrogen (degrees)')) +
	scale_y_continuous("Feature Density", limits=c(0,5), breaks=c(0, 2, 4))
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d_reflect_boundary(
  data = f,
  ids = c("sample_source", "don_chem_type_name"),
  variable = "cosBAH",
  reflect_left=T, left_boundary=0,
  reflect_right=T, right_boundary=pi, adjust=.2)

plot_id <- "hbond_vBAH_hydroxyl_acceptor_don_chem_type"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=acos(x)*180/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source)) +
	facet_wrap(~ don_chem_type_name) +
	ggtitle("HBond Sp3 Acceptor Virtual-BAH Angle by Donor Chemical Type, SeqSep > 5, B-Factor < 30\n(normalized for equal volume per unit distance)") +
	scale_x_continuous(paste('Base -- Acceptor -- Hydrogen (degrees)')) +
	scale_y_continuous("Feature Density", limits=c(0,4), breaks=c(0, 1, 2, 3))
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



###########################################################

stats <- estimate_primary_modes_1d(
	f[as.character(f$don_chem_type) != "hbdon_PBA",],
	c("sample_source", "acc_chem_type_name"), "cosBAH")
stats$primary_mode <- acos(stats$primary_mode) * 180/pi
table_id <- "hbond_vBAH_hydroxyl_acceptor_primary_mode_acc_chem_type_to_sc"
table_title <- "H-Bond BAH Angle Primary Mode By Acceptor Chem Type To Side chain Donors"
save_tables(self,
	stats, table_id,
	sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")




})) # end FeaturesAnalysis
