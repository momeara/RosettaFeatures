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

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "main_geom_conditional_comparison",
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
	CASE acc.HBChemType
		WHEN 'hbacc_IMD' THEN 'Ring' WHEN 'hbacc_IME' THEN 'Ring'
		WHEN 'hbacc_AHX' THEN 'Sp3'  WHEN 'hbacc_HXL' THEN 'Sp3'
		WHEN 'hbacc_CXA' THEN 'Sp2'  WHEN 'hbacc_CXL' THEN 'Sp2'
		WHEN 'hbacc_PBA' THEN 'Sp2'  ELSE NULL END AS acc_hybrid,
	CASE don.HBChemType
		WHEN 'hbdon_IMD' THEN 'Ring' WHEN 'hbdon_IME' THEN 'Ring'
		WHEN 'hbdon_IND' THEN 'Ring'
		WHEN 'hbdon_AHX' THEN 'Sp3'  WHEN 'hbdon_HXL' THEN 'sp3'
		WHEN 'hbdon_AMO' THEN 'Sp3+'
		WHEN 'hbdon_GDE' THEN 'Sp2+' WHEN 'hbdon_GDH' THEN 'Sp2+'
		WHEN 'hbdon_CXA' THEN 'Sp2'  WHEN 'hbdon_PBA' THEN 'Sp2'
		ELSE NULL END AS don_hybrid,
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
	don_chem_type_name = don_chem_type_name_linear(don_chem_type),
	acc_chem_type_name = acc_chem_type_name_linear(acc_chem_type),
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

do_comp <- function(
	sample_sources,
	f,
	id.vars_str,
	measure.var,
	dim_type,
	comp_funs
) {

	if(id.vars_str == ""){
		id.vars <- c()
	} else {
		id.vars <- strsplit(as.character(id.vars_str), " ")[[1]]
	}

        print(id.vars)
        print(measure.var)
        sub_f <- f[,c("sample_source", id.vars, measure.var)]
	comp_stats <- comparison_statistics(
		sample_sources,
		na.omit(sub_f, method="r"),
		id.vars, measure.var, comp_funs)

	table_id <- paste(measure.var, paste(id.vars, collapse="_"), paste(comp_funs, collapse="_"), "comparison", sep="_")
	table_title <- paste("H-Bond Distribution Comparison B-Factor < 30\nMeasure=", measure.var, " ", dim_type, " conditional on (", paste(id.vars, collapse=", "), ")", sep="")

	save_tables(self,
		comp_stats, table_id,
		sample_sources, output_dir, output_formats,
		caption=table_title, caption.placement="top")
}


measure.vars <- data.frame(
	var=c(
		"ADdist", "AHdist",
		"AHD", "BAH", "BAchi"),
	dim=c(
		"Distance", "Distance",
		"Angle", "Angle", "Torsion"))

id.vars <- data.frame(
	id.vars = c(
		"",
		"don_hybrid", "don_chem_type_name",
		"acc_hybrid", "acc_chem_type_name"))

comp_funs <- c("ref_count", "new_count", "kolmogorov_smirnov_test", "histogram_kl_divergence")

all_comps <- merge(id.vars, measure.vars, all=T)
a_ply(all_comps, 1, function(test){
	do_comp(
		sample_sources,
		f,
		as.character(test$id.vars[1]),
		as.character(test$var[1]),
		as.character(test$dim[1]),
		comp_funs)
}, .parallel=use_parallel)

})) # end FeaturesAnalysis

