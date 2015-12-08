# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

library(plyr)
source("../../plots/hbonds/hbond_geo_dim_scales.R")

feature_analyses <- c(feature_analyses, methods::new("FeaturesAnalysis",
id = "main_geom_dendrogram",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
	don_atoms.base_x AS dx, don_atoms.base_y AS dy, don_atoms.base_z AS dz,
	don_atoms.atm_x  AS hx,  don_atoms.atm_y  AS hy,  don_atoms.atm_z  AS hz,
	acc_atoms.atm_x  AS ax, acc_atoms.atm_y  AS ay, acc_atoms.atm_z  AS az,
 	acc_atoms.base_x AS abx, acc_atoms.base_y AS aby, acc_atoms.base_z AS abz,
	acc_atoms.base2_x AS ab2x, acc_atoms.base2_y AS ab2y, acc_atoms.base2_z AS ab2z,
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
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
	acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id;";

f <- query_sample_sources(sample_sources, sele)

f <- transform(f,
	don_chem_type_name = factor(don_chem_type,
		levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
			"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
		labels = c("dIMD_h", "dIME_h", "dGDE_r", "dGDH_r",
			"dAHX_y", "dHXL_st", "dIND_w", "dAMO_k", "dCXA_nq", "dPBA_bb")),
	acc_chem_type_name = 	factor(acc_chem_type,
		levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
			"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
		labels = c("aIMD_h", "aIME_h", "aAHX_y", "aHXL_st",
			"aCXA_nq", "aCXL_de", "aPBA_bb")),
	ADdist = vector_distance(cbind(dx, dy, dz), cbind(ax, ay, az)),
	AHdist = vector_distance(cbind(ax, ay, az), cbind(hx, hy, hz)),
	AHD = acos(vector_dotprod(
		vector_normalize(cbind(hx-ax, hy-ay, hz-az)),
		vector_normalize(cbind(dx-hx, dy-hy, dz-hz)))) * 180/pi,
	BAH = acos(vector_dotprod(
		vector_normalize(cbind(ax-abx, ay-aby, az-abz)),
		vector_normalize(cbind(hx-ax, hy-ay, hz-az)))) * 180/pi,
	BAchi = vector_dihedral(
		cbind(abx, aby, abz), cbind(ab2x, ab2y, ab2z),
		cbind(ax, ay, az), cbind(hx, hy, hz)) * 180/pi)

f<- na.omit(f, method="r")

dms <- dlply(f, .(sample_source), function(df) {
	distance_matrix(df, "don_chem_type_name", "AHdist", kolmogorov_smirnov_test)
})
tree_id <- "AHdist_don_chem_type"
save_trees(self, dms, tree_id, sample_sources, output_dir, output_formats)

dms <- dlply(f, .(sample_source), function(df) {
	distance_matrix(df, "acc_chem_type_name", "AHdist", kolmogorov_smirnov_test)
})
tree_id <- "AHdist_acc_chem_type"
save_trees(self, dms, tree_id, sample_sources, output_dir, output_formats)


f$chem_type <- interaction(f$don_chem_type_name, f$acc_chem_type_name, sep="_")
f_prune <- ddply(f, .(chem_type), function(df){
	if(nrow(df) < 50){
		print(paste("pruning ", df$chem_type, " because it only has ", nrow(df), " rows.", sep=""))
		return(data.frame())
	} else {
		return(df)
	}
})
f_prune$chem_type <- factor(f_prune$chem_type)

dms <- dlply(f_prune[,c("sample_source", "chem_type", "AHdist")], .(sample_source), function(df) {
	print(as.character(df$sample_source[1]))
	distance_matrix(df, "chem_type", "AHdist", kolmogorov_smirnov_test)
})
tree_id <- "AHdist_chem_type"
save_trees(self, dms, tree_id, sample_sources, output_dir, output_formats)




dms <- dlply(f, .(sample_source), function(df) {
	distance_matrix(df, "don_chem_type_name", "AHD", kolmogorov_smirnov_test)
})
tree_id <- "AHD_don_chem_type"
save_trees(self, dms, tree_id, sample_sources, output_dir, output_formats)

dms <- dlply(f, .(sample_source), function(df) {
	distance_matrix(df, "acc_chem_type_name", "AHD", kolmogorov_smirnov_test)
})
tree_id <- "AHD_acc_chem_type"
save_trees(self, dms, tree_id, sample_sources, output_dir, output_formats)

dms <- dlply(f_prune[,c("sample_source", "chem_type", "AHD")], .(sample_source), function(df) {
	print(as.character(df$sample_source[1]))
	distance_matrix(df, "chem_type", "AHD", kolmogorov_smirnov_test)
})
tree_id <- "AHD_chem_type"
save_trees(self, dms, tree_id, sample_sources, output_dir, output_formats)



dms <- dlply(f, .(sample_source), function(df) {
	distance_matrix(df, "don_chem_type_name", "BAH", kolmogorov_smirnov_test)
})
tree_id <- "BAH_don_chem_type"
save_trees(self, dms, tree_id, sample_sources, output_dir, output_formats)

dms <- dlply(f, .(sample_source), function(df) {
	distance_matrix(df, "acc_chem_type_name", "BAH", kolmogorov_smirnov_test)
})
tree_id <- "BAH_acc_chem_type"
save_trees(self, dms, tree_id, sample_sources, output_dir, output_formats)

dms <- dlply(f_prune[,c("sample_source", "chem_type", "BAH")], .(sample_source), function(df) {
	print(as.character(df$sample_source[1]))
	distance_matrix(df, "chem_type", "BAH", kolmogorov_smirnov_test)
})
tree_id <- "BAH_chem_type"
save_trees(self, dms, tree_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis

