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
id = "sequence_recovery",
author = "Matthew O'Meara",
brief_description = "",
long_description = "",
feature_reporter_dependencies = c("StructureFeatures", "ResidueFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


do_analysis <- function(sele, id, title){
	f <- query_sample_sources_against_ref(sample_sources, sele)

	#reorder the residue types by functional group similarity
	f$ref_res_type <- factor(f$ref_res_type,
		levels=c("MET", "CYS", "PRO", "TRP", "ARG", "LYS", "HIS", "ASP", "GLU", "ASN", "GLN", "PHE", "TYR", "THR", "SER", "ALA", "VAL", "ILE", "LEU", "GLY"),
		labels=c("MET", "CYS", "PRO", "TRP", "ARG", "LYS", "HIS", "ASP", "GLU", "ASN", "GLN", "PHE", "TYR", "THR", "SER", "ALA", "VAL", "ILE", "LEU", "GLY"))

	#reorder the residue types by functional group similarity
	f$new_res_type <- factor(f$new_res_type,
		levels=c("MET", "CYS", "PRO", "TRP", "ARG", "LYS", "HIS", "ASP", "GLU", "ASN", "GLN", "PHE", "TYR", "THR", "SER", "ALA", "VAL", "ILE", "LEU", "GLY"),
		labels=c("MET", "CYS", "PRO", "TRP", "ARG", "LYS", "HIS", "ASP", "GLU", "ASN", "GLN", "PHE", "TYR", "THR", "SER", "ALA", "VAL", "ILE", "LEU", "GLY"))

	f <- f[complete.cases(f),]

	wide_f <- cast(f, sample_source + ref_res_type ~ new_res_type, value="count")
	save_tables(self, wide_f, id, sample_sources, output_dir, output_formats, caption=title, caption.placement="top")

	d_ply(f, .(sample_source), function(df){
		new_sample_source <- df$new_sample_source[1]
		p_recovered <- round(sum(df[df$ref_res_type == df$new_res_type, "count"]) / sum(df$count) * 100, 2)
		plot_id <- paste(
			id, new_sample_source, sep="_")
		p <- ggplot(data=df) + theme_bw() +
			ggtitle(
				paste(title,"\nSample Source: ", new_sample_source, " Percent Recovered: ", p_recovered, "%", sep="")) +
			geom_tile(aes(x=ref_res_type, y=new_res_type, fill=log(count))) +
			scale_fill_gradientn('Log(Counts)', colours=jet.colors(15))
		save_plots(self, plot_id, sample_sources, output_dir, output_formats)
	})
}


sele <- "
SELECT
	ref_res.name3 AS ref_res_type,
	new_res.name3 AS new_res_type,
	count(*) AS count
FROM
	ref.structures AS ref_struct, new.structures AS new_struct,
	ref.residues AS ref_res, new.residues AS new_res,
	ref.residue_pdb_confidence AS ref_conf
WHERE
	ref_res.struct_id = ref_struct.struct_id AND
	ref_conf.struct_id = ref_struct.struct_id AND
	ref_conf.residue_number = ref_res.resNum AND
	ref_conf.max_temperature < 30 AND
	new_struct.tag = ref_struct.tag AND
	new_res.struct_id = new_struct.struct_id AND
	new_res.resNum = ref_res.resNum
GROUP BY
	ref_res_type, new_res_type;"
do_analysis(
	sele, "residue_mutation_counts", "Residue Mutation Counts; B-Factor < 30")

sele <- "
SELECT
	ref_res.name3 AS ref_res_type,
	new_res.name3 AS new_res_type,
	count(*) AS count
FROM
	ref.structures AS ref_struct, new.structures AS new_struct,
	ref.residues AS ref_res, new.residues AS new_res,
	ref.residue_pdb_confidence AS ref_conf,
	ref.residue_burial AS ref_res_burial
WHERE
	ref_res.struct_id = ref_struct.struct_id AND
	ref_conf.struct_id = ref_struct.struct_id AND
	ref_conf.residue_number = ref_res.resNum AND
	ref_conf.max_temperature < 30 AND
	ref_res_burial.struct_id = ref_struct.struct_id AND
	ref_res_burial.resNum = ref_res.resNum AND
	ref_res_burial.ten_a_neighbors >= 24 AND
	new_struct.tag = ref_struct.tag AND
	new_res.struct_id = new_struct.struct_id AND
	new_res.resNum = ref_res.resNum
GROUP BY
	ref_res_type, new_res_type;"
do_analysis(
	sele,
	"residue_mutation_counts_buried",
	"Residue Mutation Counts; Buried (10A Neighbors >= 24); B-Factor < 30")

sele <- "
SELECT
	ref_res.name3 AS ref_res_type,
	new_res.name3 AS new_res_type,
	count(*) AS count
FROM
	ref.structures AS ref_struct, new.structures AS new_struct,
	ref.residues AS ref_res, new.residues AS new_res,
	ref.residue_pdb_confidence AS ref_conf,
	ref.residue_burial AS ref_res_burial
WHERE
	ref_res.struct_id = ref_struct.struct_id AND
	ref_conf.struct_id = ref_struct.struct_id AND
	ref_conf.residue_number = ref_res.resNum AND
	ref_conf.max_temperature < 30 AND
	ref_res_burial.struct_id = ref_struct.struct_id AND
	ref_res_burial.resNum = ref_res.resNum AND
	ref_res_burial.ten_a_neighbors < 24 AND
	ref_res_burial.ten_a_neighbors >= 17 AND
	new_struct.tag = ref_struct.tag AND
	new_res.struct_id = new_struct.struct_id AND
	new_res.resNum = ref_res.resNum
GROUP BY
	ref_res_type, new_res_type;"
do_analysis(
	sele,
	"residue_mutation_counts_partial",
	"Residue Mutation Counts; Partially Exposed (17 <= 10A Neighbors < 24); B-Factor < 30")

sele <- "
SELECT
	ref_res.name3 AS ref_res_type,
	new_res.name3 AS new_res_type,
	count(*) AS count
FROM
	ref.structures AS ref_struct, new.structures AS new_struct,
	ref.residues AS ref_res, new.residues AS new_res,
	ref.residue_pdb_confidence AS ref_conf,
	ref.residue_burial AS ref_res_burial
WHERE
	ref_res.struct_id = ref_struct.struct_id AND
	ref_conf.struct_id = ref_struct.struct_id AND
	ref_conf.residue_number = ref_res.resNum AND
	ref_conf.max_temperature < 30 AND
	ref_res_burial.struct_id = ref_struct.struct_id AND
	ref_res_burial.resNum = ref_res.resNum AND
	ref_res_burial.ten_a_neighbors < 17 AND
	new_struct.tag = ref_struct.tag AND
	new_res.struct_id = new_struct.struct_id AND
	new_res.resNum = ref_res.resNum
GROUP BY
	ref_res_type, new_res_type;"
do_analysis(
	sele,
	"residue_mutation_counts_exposed",
	"Residue Mutation Counts; Exposed (10A Neighbors < 17); B-Factor < 30")

})) # end FeaturesAnalysis
