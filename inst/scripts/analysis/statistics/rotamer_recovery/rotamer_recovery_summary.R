# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

feature_analyses <- c(feature_analyses, methods::new("FeaturesAnalysis",
id = "rotamer_recovery_summary",
author = "Matthew O'Meara",
brief_description = "",
long_description = "",

feature_reporter_dependencies = c("ResidueFeatures", "RotamerRecoveryFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
  COUNT(rr.recovered) AS total,
  SUM(rr.recovered) AS num_recovered,
	AVG(rr.recovered)*100 AS avg_recovery
FROM
	residues AS res,
	residue_pdb_confidence AS res_conf,
	rotamer_recovery AS rr
WHERE
	res.name3 != 'ALA' AND res.name3 != 'GLY' AND
	res_conf.struct_id = res.struct_id AND
	res_conf.residue_number = res.resNum AND
	res_conf.max_sc_temperature < 30 AND
	rr.struct_id = res.struct_id AND
	rr.resNum = res.resNum;"


f <- query_sample_sources(sample_sources, sele)


sele <- "
DROP TABLE nchi;
DROP TABLE chi_recovered;"
query_sample_sources(sample_sources, sele, warn_zero_rows=F)

f$avg_recovery <- round(as.numeric(f$avg_recovery), 2)

table_id <- "rotamer_recovery_summary"
table_title <- "Average Rotamer Recovery, BFactor < 30\n"
save_tables(self,
	f, table_id,
	sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")

})) # end FeaturesAnalysis
