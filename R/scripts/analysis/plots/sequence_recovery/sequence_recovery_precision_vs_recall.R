# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

library(reshape2)
library(ggplot2)
library(plyr)

feature_analyses <- c(feature_analyses, methods::new("FeaturesAnalysis",
id = "sequence_recovery_precision_vs_recall",
author = "Matthew O'Meara",
brief_description = "",
long_description = "",
feature_reporter_dependencies = c("StructureFeatures", "ResidueFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <- "
SELECT
	CASE
		WHEN ref_res_burial.ten_a_neighbors < 17 THEN 'Exposed'
		WHEN ref_res_burial.ten_a_neighbors BETWEEN 17 AND 23 THEN 'Partial'
		ELSE 'Buried' END AS burial,
	ref_res.name3 AS ref_res_type,
	new_res.name3 AS new_res_type
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
	new_struct.tag = ref_struct.tag AND
	new_res.struct_id = new_struct.struct_id AND
	new_res.resNum = ref_res.resNum AND
	ref_res.name3 IN ('MET', 'CYS', 'PRO', 'TRP', 'ARG', 'LYS', 'HIS', 'ASP', 'GLU', 'ASN', 'GLN', 'PHE', 'TYR', 'THR', 'SER', 'ALA', 'VAL', 'ILE', 'LEU', 'GLY');"

f <- query_sample_sources_against_ref(sample_sources, sele)

f$burial <- factor(f$burial, levels=c("Buried", "Partial", "Exposed"), labels=c("Buried", "Partial", "Exposed"))

p_vs_v <- ddply(f, .(sample_source, ref_res_type), function(ref_f) {
	sample_source <- ref_f$sample_source[1]
	ref_res_type <- ref_f$ref_res_type[1]
	new_f <- f[f$sample_source == sample_source & f$new_res_type == ref_res_type,]
	data.frame(
		res_type = ref_res_type,
		precision = sum(ref_f$new_res_type == ref_res_type) / nrow(new_f),
		recall = sum(ref_f$new_res_type == ref_res_type) / nrow(ref_f))
})

table_id <- "sequence_recovery_precision_vs_recall_by_res_type"
title <- "Sequence Recovery Precision vs. Recall by Residue Type"
long_p_vs_v <- reshape2::melt(p_vs_v, c("sample_source", "res_type"), c("precision", "recall"))
wide_p_vs_v <- dcast(long_p_vs_v)
save_tables(self, wide_p_vs_v, table_id, sample_sources, output_dir, output_formats, caption=title, caption.placement="top")

plot_id <- "sequence_recovery_precision_vs_recall_by_res_type"
title <- "Sequence Recovery Precision vs. Recall by Residue Type"
p <- ggplot(p_vs_v) + theme_bw() +
	geom_text(aes(x=recall, y=precision, label=res_type), size=1.4) +
	ggtitle(title) +
	facet_wrap(~sample_source)
save_plots(self, plot_id, sample_sources, output_dir, output_formats)




p_vs_v <- ddply(f, .(sample_source, burial, ref_res_type), function(ref_f) {
	sample_source <- ref_f$sample_source[1]
	ref_burial <- ref_f$burial[1]
	ref_res_type <- ref_f$ref_res_type[1]
	new_f <- f[f$sample_source == sample_source & f$burial == ref_burial & f$new_res_type == ref_res_type,]
	print(paste(sample_source, ref_burial, ref_res_type))
	print(nrow(ref_f))
	print(nrow(new_f))
	print(sum(ref_f$new_res_type == ref_res_type))
	data.frame(
		res_type = ref_res_type,
		precision = sum(ref_f$new_res_type == ref_res_type) / nrow(new_f),
		recall = sum(ref_f$new_res_type == ref_res_type) / nrow(ref_f))
})

table_id <- "sequence_recovery_precision_vs_recall_by_burial_res_type"
title <- "Sequence Recovery Precision vs. Recall by Burial and Residue Type"
long_p_vs_v <- reshape2::melt(p_vs_v, c("sample_source", "burial", "res_type"), c("precision", "recall"))
wide_p_vs_v <- dcast(long_p_vs_v)
save_tables(self, wide_p_vs_v, table_id, sample_sources, output_dir, output_formats, caption=title, caption.placement="top")

plot_id <- "sequence_recovery_precision_vs_recall_by_burial_res_type"
title <- "Sequence Recovery Precision vs. Recall by Burial and Residue Type"
p <- ggplot(p_vs_v) + theme_bw() +
	geom_text(aes(x=recall, y=precision, label=res_type), size=1.4) +
	ggtitle(title) +
	facet_wrap(burial ~ sample_source)
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
