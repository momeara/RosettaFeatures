# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "main_geom_conditional_statistics",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){



source("../../plots/hbonds/hbond_geo_dim_scales.R")

sele <-"
SELECT
	don_atoms.base_x AS dx, don_atoms.base_y AS dy, don_atoms.base_z AS dz,
	don_atoms.atm_x  AS hx,  don_atoms.atm_y  AS hy,  don_atoms.atm_z  AS hz,
	acc_atoms.atm_x  AS ax, acc_atoms.atm_y  AS ay, acc_atoms.atm_z  AS az,
 	acc_atoms.base_x AS abx, acc_atoms.base_y AS aby, acc_atoms.base_z AS abz,
	acc_atoms.base2_x AS ab2x, acc_atoms.base2_y AS ab2y, acc_atoms.base2_z AS ab2z,
	CASE WHEN acc.HBChemType == 'hbacc_PBA' THEN 1 ELSE 0 END AS acc_bb,
	CASE WHEN don.HBChemType == 'hbdon_PBA' THEN 1 ELSE 0 END AS don_bb,
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
	don.HBChemType AS don_chem_type,
	hb.DonRank,
	hb.AccRank,
	CASE
		WHEN (don_env.sasa_r140 == 0 AND acc_env.sasa_r140 == 0) THEN 'Buried'
		WHEN (don_env.sasa_r140 > 0 AND acc_env.sasa_r140 > 0) THEN 'Exposed'
		ELSE 'Partial' END AS sasa_burial,
	CASE don.resNum - acc.resNum
		WHEN -1 THEN '-1' WHEN -2 THEN '-2' WHEN -3 THEN '-3' WHEN -4 THEN '-4'
		WHEN 1 THEN '1' WHEN 2 THEN '2' WHEN 3 THEN '3' WHEN 4 THEN '4'
		ELSE 'Long' END AS seq_sep,
		don_ss_code.label AS don_dssp,
		acc_ss_code.label AS acc_dssp,
		CASE WHEN don_ss_code.code = acc_ss_code.code THEN don_ss_code.label ELSE NULL END AS dssp
FROM
	hbonds AS hb cross join
	hbond_sites AS don cross join
	hbond_sites AS acc cross join
	hbond_sites_pdb AS don_pdb cross join
	hbond_sites_pdb AS acc_pdb cross join
	hbond_site_atoms AS don_atoms cross join
	hbond_site_atoms AS acc_atoms cross join
	hbond_site_environment AS don_env cross join
	hbond_site_environment AS acc_env cross join
	residue_secondary_structure AS don_ss cross join
	dssp_codes AS don_ss_code cross join
	residue_secondary_structure AS acc_ss cross join
	dssp_codes AS acc_ss_code
WHERE
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
	acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id AND
	don_env.struct_id = hb.struct_id AND don_env.site_id = hb.don_id AND
	acc_env.struct_id = hb.struct_id AND acc_env.site_id = hb.acc_id AND
	acc_ss.struct_id = acc.struct_id AND acc_ss.resNum = acc.resNum AND
	acc_ss_code.code = acc_ss.dssp AND
	don_ss.struct_id = don.struct_id AND don_ss.resNum = don.resNum AND
	don_ss_code.code = don_ss.dssp;";

f <- query_sample_sources(sample_sources, sele)

f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type)

f <- transform(f,
	don_chem_type_name = don_chem_type_name_linear(don_chem_type),
	acc_chem_type_name = acc_chem_type_name_linear(acc_chem_type),
	AHdist = vector_distance(cbind(ax, ay, az), cbind(hx, hy, hz)),
  AHD = acos(vector_dotprod(
    vector_normalize(cbind(hx-ax, hy-ay, hz-az)),
		vector_normalize(cbind(dx-hx, dy-hy, dz-hz)))) * 180/pi,
  BAH = acos(vector_dotprod(
    vector_normalize(cbind(ax-abx, ay-aby, az-abz)))))

f <- na.omit(f, mehtod="r")


})) # end FeaturesAnalysis

