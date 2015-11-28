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
id = "charge_charge_distances",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "ResidueTypeFeatures", "ChargeChargeFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")

sele <-"
SELECT
	CASE WHEN s1.HBChemType < s2.HBChemType
	THEN s1.HBChemType ELSE s2.HBChemType END AS s1_chem_type,

	CASE WHEN s1.HBChemType < s2.HBChemType
	THEN s2.HBChemType ELSE s1.HBChemType END AS s2_chem_type,

	cc.B1q1q2_angle,
	cc.B2q2q1_angle,
	cc.q1q2_distance,
	cc.B1q1_torsion,
	cc.B2q2_torsion
FROM
	charge_charge_pairs AS cc,
	hbond_sites AS s1,
	hbond_sites AS s2,
	hbond_sites_pdb AS s1_pdb,
	hbond_sites_pdb AS s2_pdb
WHERE
	s1.struct_id = cc.struct_id AND
	s1.site_id = cc.q1_site_id AND
	s2.struct_id = cc.struct_id AND
	s2.site_id = cc.q2_site_id AND
	s1_pdb.struct_id = cc.struct_id AND
	s1_pdb.site_id = cc.q1_site_id AND
	s2_pdb.struct_id = cc.struct_id AND
	s2_pdb.site_id = cc.q2_site_id AND
	s1_pdb.heavy_atom_temperature < 30 AND
	s2_pdb.heavy_atom_temperature < 30;"

f <-  query_sample_sources(sample_sources, sele)

f$s1_chem_type_name <- chem_type_name_linear(f$s1_chem_type)
f$s2_chem_type_name <- chem_type_name_linear(f$s2_chem_type)
f <- na.omit(f, method="r")

dens <- estimate_density_1d(f,
	c("sample_source", "s1_chem_type_name", "s2_chem_type_name"),
	"q1q2_distance", weight_fun=radial_3d_normalization)

plot_id <- "charge_charge_pair_distances"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source, group=sample_source), size=.8) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	facet_grid(s1_chem_type_name ~ s2_chem_type_name) +
	ggtitle("Charge-Charge Pair Distances; B-Fact < 30") +
	scale_x_continuous("Distance Between Polar Sites (A)", limits=c(1.5, 8)) +
	scale_y_continuous("Feature Density")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
