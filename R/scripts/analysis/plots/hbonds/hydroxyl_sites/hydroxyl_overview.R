# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

library(ggplot2)


feature_analyses <- c(feature_analyses, methods::new("FeaturesAnalysis",
id = "hydroxyl_overview",
author = "Matthew O'Meara",
brief_description = "",
long_description = "

-----------------------------------------------------------------------
|A                   |E                      |I                      |
|										 |											 |										 	 |
|			HXL OH chi		 |		DON HXL AHDist  	 |			DON HXL AHD    	 |
|										 |											 |										 	 |
|										 |											 |										 	 |
|										 |											 |										 	 |
|										 |											 |										 	 |
---------------------------------------------------------------------|
|B									 |F											 |J										 	 |
|										 |											 |										 	 |
|			AHX OH chi		 |	 	 DON AHX AHDist		 |			DON AHX  AHD	 	 |
|										 |											 |										 	 |
|										 |											 |										 	 |
|										 |											 |										 	 |
|										 |											 |										 	 |
----------------------------------------------------------------------
|C                   |G                      |K                      |
|                    |                       |                       |
|     HXL BA-Chi     |    ACC HXL AHDist     |       ACC HXL AHD     |
|                    |                       |                       |
|                    |                       |                       |
|                    |                       |                       |
|                    |                       |                       |
|--------------------|-----------------------------------------------|
|D                   |H                      |L                      |
|                    |                       |                       |
|     AHX BA-Chi     |     ACC AHX AHDist    |       ACC AHX AHD     |
|                    |                       |                       |
|                    |                       |                       |
|                    |                       |                       |
|                    |                       |                       |
|---------------------------------------------------------------------

",
feature_reporter_dependencies = c("ResidueFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){



panelsAB <- function(narrow_output_formats){
  # both donotes and accepts, where both partners are neither patner
  # is ASN, IMD, IME, HXL, AHX make sure the hydrogen is well placed
  sele <- "
	create index if not exists hbond_sites_struct_id_resNum on hbond_sites (struct_id, resNum );
	create index if not exists hbonds_struct_id_don_id on hbonds ( struct_id, don_id );
	create index if not exists hbonds_struct_id_acc_id on hbonds ( struct_id, acc_id );

	SELECT
		res_conf.chi2 AS chi2,
		res_conf.chi3 AS chi3,
		CASE main_don.HBCHemType
			WHEN 'hbdon_HXL' THEN res_conf.chi2
			WHEN 'hbdon_AHX' THEN res_conf.chi3 END AS chi,
		main_don.HBChemType AS chem_type,
		'chi' AS dof
	FROM
		hbond_sites AS main_don CROSS JOIN
		hbond_sites_pdb AS main_don_pdb CROSS JOIN
		hbonds AS main_as_don_hb CROSS JOIN
		hbond_sites AS partner_acc CROSS JOIN
		protein_residue_conformation AS res_conf
	WHERE
		main_don_pdb.struct_id = main_don.struct_id AND
		main_don_pdb.site_id = main_don.site_id AND
		main_don_pdb.heavy_atom_temperature < 30 AND

		main_as_don_hb.struct_id = main_don.struct_id AND
		main_as_don_hb.don_id = main_don.site_id AND
		partner_acc.struct_id = main_as_don_hb.struct_id AND
		partner_acc.site_id = main_as_don_hb.acc_id AND
		partner_acc.HBChemType != 'hbacc_ASN' AND
		partner_acc.HBChemType != 'hbacc_IMD' AND
		partner_acc.HBChemType != 'hbacc_IME' AND
		partner_acc.HBChemType != 'hbacc_HXL' AND
		partner_acc.HBChemType != 'hbacc_AHX' AND
		ABS(partner_acc.resNum - main_don.resNum) > 5 AND

		res_conf.struct_id = main_don.struct_id AND
		res_conf.seqpos = main_don.resNum;"
	f <- query_sample_sources(sample_sources, sele)
	f <- na.omit(f, method="r")
	f$chi <- f$chi

	dens <- estimate_density_1d_wrap(
		f, c("sample_source", "chem_type"), "chi", xlim=c(-180, 180), adjust=.6)
	summary(dens)

	dens$chem_type <- factor(
		dens$chem_type,
		levels=c("hbdon_HXL", "hbdon_AHX", "hbacc_HXL", "hbacc_AHX"),
		labels=c("hbdon_HXL", "hbdon_AHX", "hbacc_HXL", "hbacc_AHX"))

	dens$dof <- 'chi'

	plot_id <- "hydroxyl_overview_AB"
	p <- ggplot(data=dens) + theme_bw() +
		facet_grid(chem_type ~ dof, drop=F) +
#		coord_equal(ratio=1) +
		geom_line(aes(x=x, y=y*1000, colour=sample_source), show_guide=F) +
		scale_x_continuous('Protein Chi (degrees)', limits=c(-180, 180), breaks=c(-120, -60, 0, 60, 120)) +
		scale_y_continuous('Feature Density', limits=c(0, 9.5), breaks=c(0,2,4,6,8))
	save_plots(self, plot_id, sample_sources, output_dir, narrow_output_formats)
}


panelsCD <- function(narrow_output_formats) {

	sele <-"
	SELECT
		geom.chi,
		geom.cosBAH,
		acc_atoms.base_x  AS abx,  acc_atoms.base_y  AS aby,  acc_atoms.base_z  AS abz,
	  acc_atoms.base2_x AS ab2x, acc_atoms.base2_y AS ab2y, acc_atoms.base2_z AS ab2z,
	  acc_atoms.atm_x   AS ax,  acc_atoms.atm_y   AS ay,  acc_atoms.atm_z   AS az,
	  don_atoms.atm_x   AS hx,  don_atoms.atm_y   AS hy,  don_atoms.atm_z   AS hz,
		acc.HBChemType AS chem_type,
		'chi' AS dof
	FROM
		hbond_geom_coords AS geom,
		hbonds AS hb,
		hbond_sites AS don,
		hbond_sites AS acc,
		hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb,
		hbond_site_atoms AS don_atoms,
		hbond_site_atoms AS acc_atoms
	WHERE
		hb.struct_id = geom.struct_id AND hb.hbond_id = geom.hbond_id AND
		hb.struct_id = don.struct_id AND hb.don_id = don.site_id AND
		hb.struct_id = acc.struct_id AND hb.acc_id = acc.site_id AND
		don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
		don_pdb.heavy_atom_temperature < 30 AND
		acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
		acc_pdb.heavy_atom_temperature < 30 AND
		don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
		acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id AND
	  ABS(don.resNum - acc.resNum ) > 5 AND
		(acc.HBChemType == 'hbacc_HXL' OR acc.HBChemType == 'hbacc_AHX');"


	f <- query_sample_sources(sample_sources, sele)

	f <- transform(f,
		cosvBAH = vector_dotprod(
	    vector_normalize(cbind(ax-(abx+ab2x)/2, ay-(aby+ab2y)/2, az-(abz+ab2z)/2)),
			vector_normalize(cbind(hx-ax, hy-ay, hz-az))),
		vBAchi = vector_dihedral(
			cbind(ab2x, ab2y, ab2z), cbind((abx+ab2x)/2, (aby+ab2y)/2, (abz+ab2z)/2),
			cbind(ax, ay, az), cbind(hx, hy, hz)))

	#equal area projection
	f <- transform(f,
		capx = 2*sin(acos(cosvBAH)/2)*cos(vBAchi),
		capy = 2*sin(acos(cosvBAH)/2)*sin(vBAchi))

	capx_limits <- c(-1.5,1.5)
	capy_limits <- c(-1.5,1.5)

	max_BAH_angle = 110 # in degrees
	abs_cap_limit <- 2*sin((max_BAH_angle*pi/180)/2)
	capx_limits <- c(-abs_cap_limit, abs_cap_limit)
	capy_limits <- c(-abs_cap_limit, abs_cap_limit)

	plot_parts <- list(
		theme_bw(),
		theme(panel.background=element_rect(fill="#00007F", colour="#00007F")),
		geom_raster(aes(x=x, y=y, fill=z), show_guide=F),
		facet_grid(chem_type ~ dof),
		polar_equal_area_grids_bw(),
		coord_equal(ratio=1),
		scale_fill_viridis("Density"),
		scale_x_continuous('<spacer>', limits=capx_limits, breaks=c()),
		scale_y_continuous('<spacer>', limits=capy_limits, breaks=c()),
		theme(
#			axis.text.x=element_blank(),
#			axis.text.y=element_blank(),
#			axis.title.x=element_blank(),
#			axis.title.y=element_blank(),
			axis.ticks.x = element_blank(),
			axis.ticks.y = element_blank()))

	dens <- estimate_density_2d(
		f[f$sample_source == sample_sources$sample_source[1],],
		c("chem_type"), "capx", "capy", n_pts=300)

	dens$chem_type <- factor(
		dens$chem_type,
		levels=c("hbdon_HXL", "hbdon_AHX", "hbacc_HXL", "hbacc_AHX"),
		labels=c("hbdon_HXL", "hbdon_AHX", "hbacc_HXL", "hbacc_AHX"))

	dens$dof <- 'chi'
	plot_id <- "hydroxyl_overview_CD"
	ggplot(data=dens) +
		theme_bw() +
		theme(panel.background=element_rect(fill="#00007F", colour="#00007F")) +
		geom_raster(aes(x=x, y=y, fill=z)) +
		facet_grid(chem_type ~ dof, drop=F) +
		polar_equal_area_grids_bw() +
		coord_equal(ratio=1) +
		scale_fill_viridis("Density") +
		scale_x_continuous('<spacer>', limits=capx_limits, breaks=c()) +
		scale_y_continuous('<spacer>', limits=capy_limits, breaks=c()) +
		theme(
#			axis.text.x=element_blank(),
#			axis.text.y=element_blank(),
#			axis.title.x=element_blank(),
#			axis.title.y=element_blank(),
			axis.ticks.x = element_blank(),
			axis.ticks.y = element_blank())
	save_plots(self, plot_id, sample_sources,	output_dir, narrow_output_formats)
}

panelsEtoL <- function(narrow_output_formats) {
	sele <-"
	SELECT
		geom.AHdist AS AHdist,
		geom.cosAHD AS cosAHD,
		don.HBChemType AS chem_type
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
	  ABS(don.resNum - acc.resNum ) > 5 AND
		(don.HBChemType == 'hbdon_HXL' OR don.HBChemType == 'hbdon_AHX')
	UNION
	SELECT
		geom.AHdist,
		geom.cosAHD,
		acc.HBChemType AS chem_type
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
	  ABS(don.resNum - acc.resNum ) > 5 AND
		(acc.HBChemType == 'hbacc_HXL' OR acc.HBChemType == 'hbacc_AHX');"
	f <- query_sample_sources(sample_sources, sele)
	f <- na.omit(f, method="r")
	dens_AHdist <- estimate_density_1d(
		f, c("sample_source", "chem_type"),
		"AHdist", weight_fun = radial_3d_normalization)
	dens_AHdist$dof <- "AHdist"
	dens_AHdist$chem_type <- factor(
		dens_AHdist$chem_type,
		levels=c("hbdon_HXL", "hbdon_AHX", "hbacc_HXL", "hbacc_AHX"),
		labels=c("hbdon_HXL", "hbdon_AHX", "hbacc_HXL", "hbacc_AHX"))

	plot_id <- "hydroxyl_overview_EtoH"
	p <- ggplot(data=dens_AHdist) + theme_bw() +
		facet_grid(chem_type ~ dof) +
		geom_line(aes(x=x, y=y, colour=sample_source), show_guide=F) +
		scale_y_continuous("FeatureDensity", limits=c(0,9.5), breaks=c(0,2,4,6,8)) +
		scale_x_continuous(
			expression(paste('A-H Dist. (', ring(A), ')')), limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.5))
	save_plots(self, plot_id, sample_sources, output_dir, narrow_output_formats)


	f$AHD <- acos(f$cosAHD) * 180/pi
	cdf_AHD <- compute_quantiles(
		f, c("sample_source", "chem_type"), "AHD")
	cdf_AHD$dof <- "AHD"
	cdf_AHD$y <- 100 * cdf_AHD$probs
	cdf_AHD$x <- cdf_AHD$quantiles
	cdf_AHD$chem_type <- factor(
		cdf_AHD$chem_type,
		levels=c("hbdon_HXL", "hbdon_AHX", "hbacc_HXL", "hbacc_AHX"),
		labels=c("Don Hydroxyl: S,T", "Don Aromatic Hydroxyl: Y", "Acc Hydroxyl: S,T", "Acc Aromatic Hydroxyl: Y"))


	plot_id <- "hydroxyl_overview_ItoL"
	p <- ggplot(data=cdf_AHD) + theme_bw() +
		facet_grid(chem_type ~ dof, drop=F) +
		geom_line(aes(x=x, y=y, colour=sample_source), show_guide=F) +
		scale_y_continuous("Quantiles", limit=c(0,1), breaks=c(0,.25, .50, .75)) +
		scale_x_continuous("AHD Angle (degrees)")
	save_plots(self, plot_id, sample_sources, output_dir, narrow_output_formats)

}

narrow_output_formats <- output_formats
narrow_output_formats$width <- narrow_output_formats$width/3


})) # end FeaturesAnalysis
