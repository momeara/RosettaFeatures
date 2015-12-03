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
id = "sp2_BAH_chi_polar_density_beta_sheet",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


#

sele <-"
CREATE TEMPORARY TABLE ee_bb_bb_hbonds AS
SELECT
  hb.struct_id,
	hb.hbond_id,
  acc.resNum AS acc_resNum,
  don.resNum AS don_resNum
FROM
  hbonds AS hb,
  hbond_sites AS acc,
  hbond_sites AS don,
	hbond_sites_pdb AS don_pdb,
	hbond_sites_pdb AS acc_pdb,
  residue_secondary_structure AS r1ss,
  residue_secondary_structure AS r2ss
WHERE
  acc.HBChemType = 'hbacc_PBA' AND
  acc.struct_id = hb.struct_id AND
  acc.site_id   = hb.acc_id AND
  don.site_id   = hb.don_id AND
  don.HBChemType = 'hbdon_PBA' AND
  don.struct_id = acc.struct_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
  r1ss.struct_id     = acc.struct_id AND
  r1ss.resNum        = acc.resNum AND
  r1ss.dssp          = 'E' AND
  r2ss.struct_id     = don.struct_id AND
  r2ss.resNum        = don.resNum AND
  r2ss.dssp          = 'E';

CREATE TEMPORARY TABLE antiparallel_close_contact_hbonds AS
SELECT
  hb1.struct_id,
	hb1.hbond_id,
  hb1.acc_resNum AS resNum_i,
  hb1.don_resNum AS resNum_j
FROM
	ee_bb_bb_hbonds AS hb1,
  ee_bb_bb_hbonds AS hb2
WHERE
  hb1.struct_id  = hb2.struct_id AND
  hb1.acc_resNum = hb2.don_resNum AND
  hb1.don_resNum = hb2.acc_resNum;

CREATE TEMPORARY TABLE parallel_close_contact_hbonds AS
SELECT
  hb1.struct_id,
	hb1.hbond_id,
  hb1.acc_resNum AS resNum_i,
  hb1.don_resNum AS resNum_j
FROM
  ee_bb_bb_hbonds AS hb1,
  ee_bb_bb_hbonds AS hb2
WHERE
  hb1.struct_id   = hb2.struct_id AND
  (( hb1.don_resNum  = hb2.acc_resNum AND hb1.acc_resNum +2 = hb2.don_resNum     ) OR
	 ( hb1.acc_resNum  = hb2.don_resNum AND hb1.don_resNum    = hb2.acc_resNum + 2 ));

CREATE TEMPORARY TABLE antiparallel_geoms AS
SELECT geom.cosBAH, geom.chi, 'antiparallel' AS strand_orientation
FROM	 hbond_geom_coords AS geom,
			 antiparallel_close_contact_hbonds AS aphbond
WHERE	 aphbond.struct_id = geom.struct_id AND aphbond.hbond_id = geom.hbond_id;


CREATE TEMPORARY TABLE parallel_geoms AS
SELECT geom.cosBAH, geom.chi, 'parallel' AS strand_orientation
FROM   hbond_geom_coords as geom,
			 parallel_close_contact_hbonds AS phbond
WHERE	 geom.struct_id =  phbond.struct_id AND geom.hbond_id = phbond.hbond_id;

SELECT * FROM antiparallel_geoms UNION
SELECT * FROM     parallel_geoms;"

f <- RosettaFeatures::query_sample_sources(sample_sources, sele)

#equal area projection
f <- transform(f,
	capx = 2*sin(acos(cosBAH)/2)*cos(chi),
	capy = 2*sin(acos(cosBAH)/2)*sin(chi))

capx_limits <- c(-1.5,1.5); capy_limits <- capx_limits;


narrow_output_formats <- transform(output_formats, width=height/1.5)


est_dens <- function(f, id.vars, n=100, h=1){
	ddply(f, id.vars, function(sub_f){
		print(summary(sub_f))
		dens <- MASS::kde2d(sub_f$capx, sub_f$capy, n=n, h=h)
		densdf <- data.frame(expand.grid(x=dens$x, y=dens$y), z=as.vector(dens$z)/max(dens$z))
		densdf$counts <- nrow(sub_f)
		densdf
	})
}

plot_parts <- function(){
	list(
		theme_bw(),
		theme(panel.background=element_rect(fill="#00007F", colour="#00007F")),
		geom_raster(aes(x=x, y=y, fill=z)),
		polar_equal_area_grids_bw(scale=.3, label_scale=.5),
		geom_indicator(data=counts, aes(indicator=counts), color="white", group=1),
		scale_x_continuous('', limits=capx_limits, breaks=c()),
		scale_y_continuous('', limits=capy_limits, breaks=c()),
#		scale_fill_gradient('Scaled Density', low="white", high="black"),
		scale_fill_viridis("Density"),
		coord_equal(ratio=1))
}

raster_n <- 500

alt_output_formats <- output_formats
alt_output_formats$height <- 6
alt_output_formats$width <- 8
alt_output_formats$extension <- ".pdf"

plot_id = "sp2_hbond_BAH_BAchi_equal_area_log_scale_beta_sheet_by_strand_orientation_sample_source"
counts <- group_counts(f, c("strand_orientation", "sample_source"))
dens <- estimate_density_2d(f, c("strand_orientation", "sample_source"), "capx", "capy", min_count=20, n=raster_n, h=.08, scaled=TRUE)
ggplot(data=dens) + plot_parts() +
	facet_grid(strand_orientation ~ sample_source) +
	ggtitle("Beta-Sheet H-Bonds B-Factor < 30")
save_plots(self, plot_id, sample_sources, output_dir, alt_output_formats)

})) # end FeaturesAnalysis
