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
id = "hbond_aSER_dBB_seq_sep_-1_chi_170_190",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("StructureFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	acc_atoms.base_x  AS abx,  acc_atoms.base_y  AS aby,  acc_atoms.base_z  AS abz,
	acc_atoms.base2_x AS ab2x, acc_atoms.base2_y AS ab2y, acc_atoms.base2_z AS ab2z,
	acc_atoms.atm_x   AS ax,  acc_atoms.atm_y   AS ay,  acc_atoms.atm_z   AS az,
	don_atoms.atm_x   AS hx,  don_atoms.atm_y   AS hy,  don_atoms.atm_z   AS hz,
	hb.hbond_id,
	lower(hex(struct.struct_id)) AS struct_id,
	hb.hbond_id,
	struct.tag,
	lower(hex(struct.struct_id)) || '_' || acc.resNum || '_' || don.resNum AS id,
	'' AS chain,
	'CA' AS don_atom1, 'C' AS don_atom2, 'O' AS don_atom3,
	don.resNum AS don_resNum,
 	'CB' AS acc_atom1, 'OG'  AS acc_atom2, 'HG'  AS acc_atom3,
	acc.resNum AS acc_resNum
FROM
	structures as struct,
	hbonds AS hb,
	hbond_sites AS don,
	hbond_sites AS acc,
	residues AS acc_res,
	hbond_sites_pdb AS don_pdb,
	hbond_sites_pdb AS acc_pdb,
	hbond_site_atoms AS don_atoms,
	hbond_site_atoms AS acc_atoms
WHERE
	hb.struct_id = struct.struct_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.HBChemType = 'hbacc_HXL' AND
	acc_res.struct_id = acc.struct_id AND acc_res.resNum = acc.resNum AND
	acc_res.name3 = 'SER' AND
	don.HBChemType = 'hbdon_PBA' AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 20 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 20 AND
	don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
	acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id AND
	acc.resNum - don.resNum = -1;"

f <- query_sample_sources(sample_sources, sele)

print(summary(f))

f <- transform(f,
	vBAH = acos(vector_dotprod(
		vector_normalize(cbind(ax-(abx+ab2x)/2, ay-(aby+ab2y)/2, az-(abz+ab2z)/2)),
		vector_normalize(cbind(hx-ax, hy-ay, hz-az)))),
	vBAchi = vector_dihedral(
		cbind(ab2x, ab2y, ab2z), cbind((abx+ab2x)/2, (aby+ab2y)/2, (abz+ab2z)/2),
		cbind(ax, ay, az), cbind(hx, hy, hz)))

cat("Filter for 170 < vBAchi < 190:\n")

f <- f[f$vBAchi < pi/180 * 190 & f$vBAchi > pi/180 * 170,]

print(summary(f))

if(nrow(f) == 0){
	cat("WARNING: Query returned no rows. Skipping rest of features analysis.\n")
	return()
}

f_orig <- f
f_orig$sample_source <- "NativeRelaxDun10Score12"
f <- head(f_orig, n=10)

ss_ids <- as.character(unique(f$sample_source))


# f:
#
#           sample_source, struct_id, id, chain, CA, C, O, acc_resNum, don_resNum
# instance1
#   ...
#

# instance_atoms:
#
#                  sample_source, struct_id, id, chain, atom
# instance1, atom1
# instance1, atom2
#   ...



don_atoms <-
	reshape2::melt(f,
		id.vars=c("id", "sample_source", "struct_id", "chain", "don_resNum"),
		measure.vars=c("don_atom1", "don_atom2", "don_atom3"),
		value.name = "atom_name")
names(don_atoms)[5] <- "resNum"
names(don_atoms)[7] <- "atom"

acc_atoms <-
	reshape2::melt(f,
		id.vars=c("id", "sample_source", "struct_id", "chain", "acc_resNum"),
		measure.vars=c("acc_atom1", "acc_atom2", "acc_atom3"),
		value.name = "atom_name")
names(acc_atoms)[5] <- "resNum"
names(acc_atoms)[7] <- "atom"

instance_atoms <- rbind(don_atoms, acc_atoms)

instances_id <- "hbond_aSER_dBB_seq_sep_-1_chi_170_190"

print(summary(instance_atoms))

prepare_feature_instances(instances_id, sample_sources, instance_atoms, output_dir)

})) # end FeaturesAnalysis
