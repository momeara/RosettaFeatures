# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "rotamer_recovery_by_salt_bridge_geometry",
author = "Matthew O'Meara",
brief_description = "",
long_description = "
This features analysis requires the comparison to be the RRComparerRotBins
",

feature_reporter_dependencies = c("ResidueFeatures", "RotamerRecoveryFeatures", "PdbDataFeatures", "SaltBridgeFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
CREATE TABLE IF NOT EXISTS CXL_ARG_salt_bridges AS SELECT
	struct.tag AS tag,
	sb.struct_id AS struct_id,
	sb.don_resNum AS don_residue_number,
	sb.acc_id AS acc_id,
	acc.resNum AS acc_residue_number,
	sb.psi * 180/3.14159 AS psi
FROM
	structures AS struct,
	salt_bridges AS sb,
	residues AS don,
	hbond_sites AS acc
WHERE
	sb.struct_id = struct.struct_id AND
	don.struct_id = sb.struct_id AND
	don.resNum = sb.don_resNum AND
	don.name3 = 'ARG' AND
	acc.struct_id = sb.struct_id AND
	acc.site_id = sb.acc_id AND
	acc.HBChemType = 'hbacc_CXL' AND
	ABS(don.resNum - acc.resNum) > 4;


SELECT
	sb.struct_id,
	sb.don_residue_number,
	sb.acc_residue_number,
	CASE WHEN -120 < sb.psi AND sb.psi < -116 THEN 1 ELSE 0 END AS split,
	don_rr.divergence AS don_rr,
	acc_rr.divergence AS acc_rr
FROM
	CXL_ARG_salt_bridges AS sb,
	rotamer_recovery AS don_rr,
	rotamer_recovery AS acc_rr,
	residue_pdb_confidence AS don_conf,
	residue_pdb_confidence AS acc_conf
WHERE
	don_rr.struct_id = sb.struct_id AND
	don_rr.resNum = sb.don_residue_number AND
	acc_rr.struct_id = sb.struct_id AND
	acc_rr.resNum = sb.acc_residue_number AND
	don_conf.struct_id = don_rr.struct_id AND
	don_conf.residue_number = don_rr.resNum AND
	don_conf.max_sc_temperature < 30 AND
	acc_conf.struct_id = acc_rr.struct_id AND
	acc_conf.residue_number = acc_rr.resNum AND
	acc_conf.max_sc_temperature < 30;"


f <- query_sample_sources(sample_sources, sele)

print(summary(f))

cast(f, sample_source + split ~ . , value="don_rr")

})) # end FeaturesAnalysis
