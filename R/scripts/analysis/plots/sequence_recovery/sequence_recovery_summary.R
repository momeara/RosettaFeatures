# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

library(reshape2)
library(plyr)


feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "sequence_recovery_summary",
author = "Matthew O'Meara",
brief_description = "",
long_description = "",
feature_reporter_dependencies = c("StructureFeatures", "ResidueFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
sele <- "
SELECT
	ref_res.name3 == new_res.name3 AS recovered,
	CASE
		WHEN ref_res_burial.ten_a_neighbors < 17 THEN 'Exposed'
		WHEN ref_res_burial.ten_a_neighbors BETWEEN 17 AND 23 THEN 'Partial'
		ELSE 'Buried' END AS burial,
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
	new_struct.tag = ref_struct.tag AND
	new_res.struct_id = new_struct.struct_id AND
	new_res.resNum = ref_res.resNum
GROUP BY
	recovered, burial;"

f <- query_sample_sources_against_ref(sample_sources, sele)

f$burial <- factor(f$burial, levels=c("Buried", "Partial", "Exposed"), labels=c("Buried", "Partial", "Exposed"))
z <- ddply(f, .(sample_source, burial), function(df) {
	data.frame(p_recovered = df[df$recovered == 1, "count"][1] * 1.0/ sum(df$count) * 100)
})

z <- rbind(z, ddply(f, .(sample_source), function(df) {
  data.frame(
    burial="All",
    p_recovered=sum(df[df$recovered == 1, "count"] * 1.0/ sum(df$count)) * 100)
}))


table_id <- "recovery_by_burial"
z_wide <- dcast(z, sample_source ~ burial, value="p_recovered")
save_tables(self,
  z_wide, table_id, sample_sources, output_dir, output_formats,
  caption="Sequence Recovery Summary By Burial", caption.placement="top")

})) # end FeaturesAnalysis
