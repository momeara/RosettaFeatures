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
id = "out_of_plane_chem_type",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){



sele <-"
SELECT
	geom.AHdist,
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
	hbonds AS hb cross join
	hbond_sites_pdb AS don_pdb cross join
	hbond_sites_pdb AS acc_pdb cross join
	hbond_sites AS don cross join
	hbond_sites AS acc cross join
	hbond_geom_coords AS geom cross join
	hbond_site_atoms AS don_atoms cross join
	hbond_site_atoms AS acc_atoms
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

f$don_chem_type_name <- don_chem_type_name_wrap(f$don_chem_type)
f$acc_chem_type_name <- acc_chem_type_name_wrap(f$acc_chem_type)
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

f$out_of_plane <- cos(f$chi) * f$AHdist
f$out_of_plane_perp <- sin(f$chi) * f$AHdist

plot_parts <- list(
	theme_bw(),
	scale_x_continuous(
		'Out of Ideal Acceptor -- Lone Pair Plane'),
	scale_y_continuous('Feature Density', breaks=c(0,2,4,6)),
	theme(legend.position=c(.58,.35)),
	theme(legend.justification=c("left", "top")))

dens <- estimate_density_1d(
	f, c("sample_source", "hybrid", "don_chem_type_name"), "out_of_plane")

l_ply(levels(dens$hybrid), function(hybrid){
	plot_id = paste("hbond_out_of_plane_", hybrid, "_acceptor_by_don_chem_type", sep="")
	ggplot(data=dens[dens$hybrid==hybrid,]) + plot_parts +
		geom_line(aes(x=x, y=y, colour=sample_source)) +
		geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
		facet_wrap( ~ don_chem_type_name) +
		ggtitle(paste("HBonds Out of Plane for ", hybrid, " acceptors, SeqSep > 5, B-Factor < 30")) +
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})

dens <- estimate_density_1d(
	f, c("sample_source", "hybrid", "don_chem_type_name"), "out_of_plane_perp")

l_ply(levels(dens$hybrid), function(hybrid){
	plot_id = paste("hbond_out_of_plane_perp_", hybrid, "_acceptor_by_don_chem_type", sep="")
	ggplot(data=dens[dens$hybrid==hybrid,]) + plot_parts +
		geom_line(aes(x=x, y=y, colour=sample_source)) +
		geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
		facet_wrap( ~ don_chem_type_name) +
		ggtitle(paste("HBonds Out of Perpendicular Plane for ", hybrid, " acceptors, SeqSep > 5, B-Factor < 30")) +
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})


f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type)
f <- na.omit(f, method="r")

plot_id = "hbond_out_of_plane_chem_type"
dens <- estimate_density_1d(
	f, c("sample_source", "acc_chem_type_name", "don_chem_type_name"), "out_of_plane")
p <- ggplot(data=dens) + plot_parts +
	geom_line(aes(x=x, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_grid(don_chem_type_name ~ acc_chem_type_name) +
	ggtitle("H-Bond Out of Plane by Chemical Type, B-Factor < 30\n(normalized for equal volume per unit distance)")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id = "hbond_out_of_plane_perp_chem_type"
dens <- estimate_density_1d(
	f, c("sample_source", "acc_chem_type_name", "don_chem_type_name"), "out_of_plane_perp")
p <- ggplot(data=dens) + plot_parts +
	geom_line(aes(x=x, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_grid(don_chem_type_name ~ acc_chem_type_name) +
	ggtitle("H-Bond Out of Perpendicular Plane by Chemical Type, B-Factor < 30\n(normalized for equal volume per unit distance)")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)




ddply(f, .(sample_source, hybrid), function(df){
	ss <- as.character(df$sample_source[1])
	hybrid <- as.character(df$hybrid[1])

	plot_id = paste("hbond_out_of_plane_chem_type_by_hybrid", hybrid, ss, sep="_")
	dens <- estimate_density_1d(
		df, c("acc_chem_type_name", "don_chem_type_name"), "out_of_plane")

	p <- ggplot(data=dens) + plot_parts +
		geom_line(aes(x=x, y=y, colour=acc_chem_type_name)) +
		geom_indicator(aes(indicator=counts, colour=acc_chem_type_name, group=acc_chem_type_name)) +
		facet_wrap(~don_chem_type_name) +
		ggtitle(paste("H-Bond Acc Hybrid: ", hybrid, " Out of Plane by Donor Chemical type\nB-Factor < 30 SampleSource: ", ss, sep=""))
	if(nrow(sample_sources) <= 3){
		p <- p + theme(legend.position="bottom", legend.direction="horizontal")
	}

	save_plots(self, plot_id, sample_sources, output_dir, output_formats)



	plot_id = paste("hbond_out_of_plane_perp_chem_type_by_hybrid", hybrid, ss, sep="_")
	dens <- estimate_density_1d(
		df, c("acc_chem_type_name", "don_chem_type_name"), "out_of_plane_perp")

	p <- ggplot(data=dens) + plot_parts +
		geom_line(aes(x=x, y=y, colour=acc_chem_type_name)) +
		geom_indicator(aes(indicator=counts, colour=acc_chem_type_name, group=acc_chem_type_name)) +
		facet_wrap(~don_chem_type_name) +
		ggtitle(paste("H-Bond Acc Hybrid: ", hybrid, " Out of Perpendicular Plane by Donor Chemical type\nB-Factor < 30 SampleSource: ", ss, sep=""))
	if(nrow(sample_sources) <= 3){
		p <- p + theme(legend.position="bottom", legend.direction="horizontal")
	}

	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})

f$all_sp2 <- as.character(f$hybrid) == "sp2" | as.character(f$hybrid) == "bb_sp2"


ddply(f[f$all_sp2,], .(sample_source), function(df){
	ss <- as.character(df$sample_source[1])
	hybrid <- "All_Sp2"

	plot_id = paste("hbond_out_of_plane_chem_type_by_hybrid", hybrid, ss, sep="_")
	dens <- estimate_density_1d(
		df, c("acc_chem_type_name", "don_chem_type_name"), "out_of_plane")

	p <- ggplot(data=dens) + plot_parts +
		geom_line(aes(x=x, y=y, colour=acc_chem_type_name)) +
		geom_indicator(aes(indicator=counts, colour=acc_chem_type_name, group=acc_chem_type_name)) +
		facet_wrap(~don_chem_type_name) +
		ggtitle(paste("H-Bond Acc Hybrid: Sp2 Out of Plane by Donor Chemical type\nB-Factor < 30 SampleSource: ", ss, sep=""))
	if(nrow(sample_sources) <= 3){
		p <- p + theme(legend.position="bottom", legend.direction="horizontal")
	}

	save_plots(self, plot_id, sample_sources, output_dir, output_formats)


	plot_id = paste("hbond_out_of_plane_perp_chem_type_by_hybrid", hybrid, ss, sep="_")
	dens <- estimate_density_1d(
		df, c("acc_chem_type_name", "don_chem_type_name"), "out_of_plane_perp")

	p <- ggplot(data=dens) + plot_parts +
		geom_line(aes(x=x, y=y, colour=acc_chem_type_name)) +
		geom_indicator(aes(indicator=counts, colour=acc_chem_type_name, group=acc_chem_type_name)) +
		facet_wrap(~don_chem_type_name) +
		ggtitle(paste("H-Bond Acc Hybrid: Sp2 Out of Perpendicular Plane by Donor Chemical type\nB-Factor < 30 SampleSource: ", ss, sep=""))
	if(nrow(sample_sources) <= 3){
		p <- p + theme(legend.position="bottom", legend.direction="horizontal")
	}

	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})


})) # end FeaturesAnalysis
