# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()
feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "chi_BAH_eq_polar_density",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")

sele <-"
SELECT
	geom.chi,
	acc_atoms.base_x  AS bx,  acc_atoms.base_y  AS by,  acc_atoms.base_z  AS bz,
	acc_atoms.base2_x AS b2x, acc_atoms.base2_y AS b2y, acc_atoms.base2_z AS b2z,
	acc_atoms.atm_x   AS ax,  acc_atoms.atm_y   AS ay,  acc_atoms.atm_z   AS az,
	don_atoms.atm_x   AS hx,  don_atoms.atm_y   AS hy,  don_atoms.atm_z   AS hz,
	acc_site.HBChemType AS acc_chem_type,
	don_site.HBChemType AS don_chem_type,
	CASE acc_site.HBChemType
		WHEN 'hbacc_IMD' THEN 'ring' WHEN 'hbacc_IME' THEN 'ring'
		WHEN 'hbacc_AHX' THEN 'sp3'  WHEN 'hbacc_HXL' THEN 'sp3'
		WHEN 'hbacc_CXA' THEN 'sp2'  WHEN 'hbacc_CXL' THEN 'sp2'
		WHEN 'hbacc_PBA' THEN 'bb_sp2'  END AS hybrid
FROM
	hbond_geom_coords AS geom,
	hbonds AS hb,
	hbond_sites AS don_site,
	hbond_sites AS acc_site,
	hbond_site_atoms AS don_atoms,
	hbond_site_atoms AS acc_atoms
WHERE
	hb.struct_id = geom.struct_id AND hb.hbond_id = geom.hbond_id AND
	hb.struct_id = don_site.struct_id AND hb.don_id = don_site.site_id AND
	hb.struct_id = acc_site.struct_id AND hb.acc_id = acc_site.site_id AND
	don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
	acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id AND
	ABS(don_site.resNum - acc_site.resNum) > 5;"

f <- query_sample_sources(sample_sources, sele)

f <- transform(f,
	don_chem_type_name = don_chem_type_name_linear(don_chem_type),
	acc_chem_type_name = acc_chem_type_name_linear(acc_chem_type),
	cosBAH = ifelse(hybrid %in% c("sp3", "ring"),
		vector_dotprod(
			vector_normalize(cbind(ax-(bx+b2x)/2, ay-(by+b2y)/2, az-(bz+b2z)/2)),
			vector_normalize(cbind(hx-ax, hy-ay, hz-az))),
		vector_dotprod(
			vector_normalize(cbind(ax-bx, ay-by, az-bz)),
			vector_normalize(cbind(hx-ax, hy-ay, hz-az)))),
	chi = ifelse(hybrid %in% c("sp3", "ring"),
		vector_dihedral(
			cbind(b2x, b2y, b2z),
			cbind((bx+b2x)/2, (by+b2y)/2, (bz+b2z)/2),
			cbind(ax, ay, az),
			cbind(hx, hy, hz)),
		vector_dihedral(
			cbind(b2x, b2y, b2z),
			cbind(bx, by, bz),
			cbind(ax, ay, az),
			cbind(hx, hy, hz))),
	hybrid = factor(
		hybrid,
		levels=c("ring", "sp3", "sp2", "bb_sp2"),
		labels=c("ring", "Sp3", "Sp2", "BB Sp2")))
f <- na.omit(f, method="r")


#equal area projection
f <- transform(f,
	capx = 2*sin(acos(cosBAH)/2)*cos(chi),
	capy = 2*sin(acos(cosBAH)/2)*sin(chi))

capx_limits <- c(-1.5,1.5)
capy_limits <- c(-1.5,1.5)

max_BAH_angle = 110 # in degrees
capx_limits <- c(-1.5, 1.5)
capy_limits <- c(-1.5, 1.5)

narrow_output_formats <- transform(output_formats, width=height)

plot_parts <- list(
	theme_bw(),
	theme(panel.background=element_rect(fill="#00007F", colour="#00007F")),
	geom_raster(aes(x=x, y=y, fill=z)),
	geom_indicator(aes(indicator=counts), color="white", group=1),
	polar_equal_area_grids_bw(scale=.4, label_scale=.6),
	coord_equal(ratio=1),
	scale_fill_gradientn('Density', colours=jet.colors(15)),
	scale_x_continuous('', limits=capx_limits, breaks=c()),
	scale_y_continuous('', limits=capy_limits, breaks=c()),
	theme(
		axis.text.x = element_blank(),
		axis.text.y = element_blank(),
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		axis.ticks.x = element_blank(),
		axis.ticks.y = element_blank()))

dens <- estimate_density_2d(
	f, c("sample_source", "hybrid"), "capx", "capy", n_pts=500, scaled=T)

plot_id <- "hbond_chi_BAH_eq_polar_density"
d_ply(dens, .(sample_source), function(sub_dens){
	ss_id <- sub_dens$sample_source[1]
	hybrid <- sub_dens$hybrid[1]
	sub_plot_id <- paste(plot_id, ss_id, sep="_")
	ggplot(data=sub_dens) + plot_parts +
		facet_wrap(~hybrid, nrow=2) +
		polar_equal_area_grids_bw() +
		ggtitle(paste(hybrid, " Acceptor H-Bonds SeqSep > 5: chi vs BAH; ss_id: ", ss_id, sep="")) +
	save_plots(self, sub_plot_id, sample_sources,	output_dir, narrow_output_formats)
})

plot_id <- "hbond_chi_BAH_eq_polar_density"
d_ply(dens, .(sample_source, hybrid), function(sub_dens){
	ss_id <- sub_dens$sample_source[1]
	hybrid <- sub_dens$hybrid[1]
	sub_plot_id <- paste(plot_id, ss_id, hybrid, sep="_")
	ggplot(data=sub_dens) + plot_parts +
		polar_equal_area_grids_bw() +
		ggtitle(paste(hybrid, " Acceptor H-Bonds SeqSep > 5: chi vs BAH; ss_id: ", ss_id, sep="")) +
	save_plots(self, sub_plot_id, sample_sources,	output_dir, narrow_output_formats)
})





dens <- estimate_density_2d(
	f, c("sample_source", "acc_chem_type_name", "don_chem_type_name"), "capx", "capy")

plot_id = "hbond_chi_BAH_eq_polar_density_by_chem_type"
ddply(dens, c("sample_source"), function(sub_dens){
	ss_id <- sub_dens$sample_source[1]
	sub_plot_id <- paste(plot_id, ss_id, sep="_")
	ggplot(data=sub_dens) + plot_parts +
		polar_equal_area_grids_bw(scale=.3, label_scale=.5) +
		facet_grid(acc_chem_type_name ~ don_chem_type_name) +
		ggtitle(paste("H-Bonds SeqSep > 5: chi vs BAH by Chemical Type; ss_id: ", ss_id, sep="")) +
	save_plots(self, sub_plot_id, sample_sources, output_dir, output_formats)
})


plot_id = "hbond_chi_BAH_eq_polar_density"
d_ply(dens, .(sample_source, don_chem_type_name, acc_chem_type_name), function(sub_dens){
	ss_id <- sub_dens$sample_source[1]

	sub_plot_id <- paste(
		plot_id, ss_id, sub_dens$don_chem_type_name[1], sub_dens$acc_chem_type_name[1], sep="_")

	ggplot(data=sub_dens) + plot_parts +
		polar_equal_area_grids_bw() +
		ggtitle(paste(
			"H-Bonds: chi vs BAH SeqSep > 5 ",
			"Don:", sub_dens$don_chem_type_name[1], " ",
			"Acc:", sub_dens$acc_chem_type_name[1], "\n",
			"ss_id: ", ss_id, sep=""))
	save_plots(self, sub_plot_id, sample_sources, output_dir, output_formats)
})


})) # end FeaturesAnalysis
