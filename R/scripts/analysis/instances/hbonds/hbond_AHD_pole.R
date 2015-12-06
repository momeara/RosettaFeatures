# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

library(reshape2)


feature_analyses <- c(feature_analyses, methods::new("FeaturesAnalysis",
id = "hbond_AHD_pole",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("StructureFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
	geom.cosAHD,
	hb.hbond_id,
	lower(hex(struct.struct_id)) AS struct_id,
	hb.hbond_id,
	struct.tag,
	acc_res_pdb.chain_id || acc_res_pdb.pdb_residue_number ||
		(CASE WHEN acc_res_pdb.insertion_code = ' ' THEN '' ELSE acc_res_pdb.insertion_code END) ||
		'_' ||
		don_res_pdb.chain_id || don_res_pdb.pdb_residue_number ||
		(CASE WHEN don_res_pdb.insertion_code = ' ' THEN '' ELSE don_res_pdb.insertion_code END)
		AS id,
	don_res_pdb.chain_id AS don_chain,
	don_res_pdb.residue_number AS don_residue_number,
	(CASE WHEN don_res_pdb.insertion_code = ' ' THEN '' ELSE don_res_pdb.insertion_code END) AS don_insertion_code,
	'H' AS don_atom1, 'N' AS don_atom2, 'CA' AS don_atom3,
	acc_res_pdb.chain_id AS acc_chain,
	acc_res_pdb.residue_number AS acc_residue_number,
	(CASE WHEN acc_res_pdb.insertion_code = ' ' THEN '' ELSE acc_res_pdb.insertion_code END) AS acc_insertion_code,
 	'O' AS acc_atom1, 'C'  AS acc_atom2, 'CA'  AS acc_atom3
FROM
	structures as struct,
	hbonds AS hb,
	hbond_geom_coords AS geom,
	hbond_sites AS don,
	hbond_sites AS acc,
	residue_pdb_identification AS don_res_pdb,
	residue_pdb_identification AS acc_res_pdb,
	hbond_sites_pdb AS don_pdb,
	hbond_sites_pdb AS acc_pdb
WHERE
	hb.struct_id = struct.struct_id AND
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc_res_pdb.struct_id = acc.struct_id AND
	acc_res_pdb.residue_number = acc.resNum AND
	don_res_pdb.struct_id = don.struct_id AND
	don_res_pdb.residue_number = don.resNum AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 20 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 20 AND
	geom.cosAHD > .999997
LIMIT
	20;"

f <- query_sample_sources(sample_sources, sele)

print(f)

if(nrow(f) == 0){
	cat("WARNING: Query returned no rows. Skipping rest of features analysis.\n")
	return()
}

ss_ids <- as.character(unique(f$sample_source))


don_atoms <-
	reshape2::melt(f,
		id.vars=
			c("id", "sample_source",
				"tag", "don_chain", "don_residue_number", "don_insertion_code"),
		measure.vars=c("don_atom1", "don_atom2", "don_atom3"),
		value.xname = "atom_name")
names(don_atoms)[4] <- "chain"
names(don_atoms)[5] <- "residue_number"
names(don_atoms)[6] <- "insertion_code"
names(don_atoms)[8] <- "atom"

acc_atoms <-
	reshape2::melt(f,
		id.vars=
			c("id", "sample_source", "tag",
				"acc_chain", "acc_residue_number", "acc_insertion_code"),
		measure.vars=c("acc_atom1", "acc_atom2", "acc_atom3"),
		value.name = "atom_name")
names(acc_atoms)[4] <- "chain"
names(acc_atoms)[5] <- "residue_number"
names(acc_atoms)[6] <- "insertion_code"
names(acc_atoms)[8] <- "atom"

instance_atoms <- rbind(don_atoms, acc_atoms)

instances_id <- "hbond_AHD_pole"

print(summary(instance_atoms))

prepare_feature_instances(instances_id, sample_sources, instance_atoms, output_dir)

})) # end FeaturesAnalysis
