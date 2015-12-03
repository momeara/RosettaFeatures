# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

library(plyr)


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
	ADdist = vector_distance(cbind(dx, dy, dz), cbind(ax, ay, az)),
	AHdist = vector_distance(cbind(ax, ay, az), cbind(hx, hy, hz)),
  AHD = acos(vector_dotprod(
    vector_normalize(cbind(hx-ax, hy-ay, hz-az)),
		vector_normalize(cbind(dx-hx, dy-hy, dz-hz)))) * 180/pi,
  BAH = acos(vector_dotprod(
    vector_normalize(cbind(ax-abx, ay-aby, az-abz)),
		vector_normalize(cbind(hx-ax, hy-ay, hz-az)))) * 180/pi,
  B2AH = acos(vector_dotprod(
    vector_normalize(cbind(ax-ab2x, ay-ab2y, az-ab2z)),
		vector_normalize(cbind(hx-ax, hy-ay, hz-az)))) * 180/pi,
	vBAH = acos(vector_dotprod(
    vector_normalize(cbind(ax-(abx+ab2x)/2, ay-(aby+ab2y)/2, az-(abz+ab2z)/2)),
		vector_normalize(cbind(hx-ax, hy-ay, hz-az)))) * 180/pi,
	B2Achi = vector_dihedral(
		cbind(ab2x, ab2y, ab2z), cbind(abx, aby, abz),
		cbind(ax, ay, az), cbind(hx, hy, hz)) * 180/pi,
	BAchi = vector_dihedral(
		cbind(abx, aby, abz), cbind(ab2x, ab2y, ab2z),
		cbind(ax, ay, az), cbind(hx, hy, hz)) * 180/pi,
	vBAchi = vector_dihedral(
		cbind(ab2x, ab2y, ab2z), cbind((abx+ab2x)/2, (aby+ab2y)/2, (abz+ab2z)/2),
		cbind(ax, ay, az), cbind(hx, hy, hz)) * 180/pi,
	AHchi = vector_dihedral(
		cbind(abx, aby, abz), cbind(ax, ay, az),
		cbind(hx, hy, hz), cbind(dx, dy, dz)) * 180/pi,
	HOutOfPlane = vector_dotprod(
		vector_normalize(cbind(hx, hy, hz) - cbind(ax, ay, az)),
		vector_crossprod(
			vector_normalize(cbind(ab2x, ab2y, ab2z) - cbind(abx, aby, abz)),
			vector_normalize(cbind(ax, ay, az) - cbind(abx, aby, abz)))))

f <- na.omit(f, mehtod="r")

do_stats <- function(
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

	cat("id.vars: ", id.vars_str, "\n", sep="")
	cat("measure.var: ", measure.var, " (", dim_type, ")\n", sep="")
	tryCatch({
		stats <- ddply(
			na.omit(f[,c("sample_source", id.vars, measure.var)], method="r"),
			c("sample_source", id.vars),
			function(df) {
				print(summary(df))
				mean_stat <- mean(df[,measure.var])
				median_stat <- median(df[,measure.var])
				primary_mode_stat <- NA
#				tryCatch({
#					primary_mode_stat <- estimate_primary_modes_1d(df, c(), measure.var)$primary_mode[1]
#				}, error=function(e) {
#					cat(paste(
#						"ERROR: Generating Pimary Mode:\n",
#						e, sep=""))
#					primary_mode_stat <- NA
#				})
				sd_stat <- sd(df[,measure.var])

				cat(mean_stat, median_stat, primary_mode_stat, sd_stat, "\n")
				print(paste(mean_stat, median_stat, primary_mode_stat, sd_stat, "\n"))
				cat(class(mean_stat), class(median_stat), class(primary_mode_stat), class(sd_stat), "\n")

				data.frame(
					mean = mean_stat,
					median = median_stat,
#					primary_mode = primary_mode_stat,
					sd = sd_stat)
		})

		table_id <- paste(measure.var, paste(id.vars, collapse="_"), "statistics", sep="_")
		table_title <- paste("H-Bond Distribution statistics B-Factor < 30\nMeasure=", measure.var, " ", dim_type, " conditional on (", paste(id.vars, collapse=", "), ")", sep="")

		save_tables(self,
			stats, table_id,
			sample_sources, output_dir, output_formats,
			caption=table_title, caption.placement="top")
	}, error=function(e){
		cat(paste(
			"ERROR: Generating statistics:\n",
			e, sep=""))
	})
}


measure.vars <- data.frame(
	var=c(
#		"ADdist",
#		"AHdist",
		"AHD", "BAH", "B2AH", "vBAH",
#		"B2Achi", "BAchi", "vBAchi", "AHchi",
		"HOutOfPlane"),
	dim=c(
#		"Distance",
#		"Distance",
		"Angle", "Angle", "Angle", "Angle",
#		"Torsion", "Torsion", "Torsion", "Torsion",
		"Distance"))

id.vars <- data.frame(
	id.vars = c(
#		"",
#		"don_bb", "don_hybrid",
		"don_chem_type_name",
#		"acc_bb", "acc_hybrid",
		"acc_chem_type_name",
		"don_hybrid acc_chem_type_name", "acc_hybrid don_chem_type_name",
		"don_chem_type_name acc_chem_type_name",
		"sasa_burial", "don_chem_type_name sasa_burial", "acc_chem_type_name sasa_burial", "don_chem_type_name acc_chem_type_name sasa_burial",
		"seq_sep", "don_chem_type_name seq_sep", "acc_chem_type_name seq_sep", "don_chem_type_name acc_chem_type_name seq_sep",
		"seq_sep", "don_chem_type_name seq_sep", "acc_chem_type_name seq_sep", "don_chem_type_name acc_chem_type_name seq_sep"))

comp_funs <- c("ref_count", "new_count", "kolmogorov_smirnov_test", "histogram_kl_divergence")

all_comps <- merge(id.vars, measure.vars, all=T)
a_ply(all_comps, 1, function(test){
	do_stats(
		sample_sources,
		f,
		as.character(test$id.vars[1]),
		as.character(test$var[1]),
		as.character(test$dim[1]),
		comp_funs)
}, .parallel=use_parallel)

})) # end FeaturesAnalysis

