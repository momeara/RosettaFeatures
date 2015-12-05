# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

library(ggplot2)
library(plyr)
library(reshape2)
source("../hbond_geo_dim_scales.R")

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "num_hbonds_normalized",
author = "Matthew O'Meara",
brief_description = "Count the number of hydrogen bonds formed conditional on the donor and acceptor chemical types.",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){



do_analysis <- function(group_sele, bond_sele, id, title){
	group_counts <- query_sample_sources(sample_sources, group_sele)
	summary(group_counts)
	group_counts[as.character(group_counts$chem_type) == "hbdon_AMO","count"] <-
			group_counts[as.character(group_counts$chem_type) == "hbdon_AMO","count"] / 3

 	group_counts[as.character(group_counts$chem_type) == "hbdon_GDH","count"] <-
		group_counts[as.character(group_counts$chem_type) == "hbdon_GDH","count"] / 2

 	group_counts[as.character(group_counts$chem_type) == "hbdon_CXA","count"] <-
		group_counts[as.character(group_counts$chem_type) == "hbdon_CXA","count"] / 2

	table_id <- paste(id, "group_counts", sep="_")
	table_title <- paste(title, "Group Counts", sep="\n")
	wide_group_counts <- dcast(group_counts, sample_source ~ chem_type, value="count")
	save_tables( self,
		wide_group_counts, table_id,
		sample_sources, output_dir, output_formats,
		caption=table_title, caption.placement="top")

	bond_counts <- query_sample_sources(sample_sources, bond_sele)

	#hacky filter for non-protein hbonds
	bond_counts <- bond_counts[bond_counts$count > 30,]

	table_id <- paste(id, "bond_counts", sep="_")
	table_title <- paste(title, "Bond Counts", sep="\n")
	wide_bond_counts <- dcast(
		bond_counts, sample_source + don_chem_type ~ acc_chem_type, value="count")
	save_tables( self,
		wide_bond_counts, table_id,
		sample_sources, output_dir, output_formats,
		caption=table_title, caption.placement="top")


	total_don <- sum(group_counts[group_counts$is_donor,"count"])
	total_acc <- sum(group_counts[!group_counts$is_donor,"count"])
	total_bonds <- sum(bond_counts$count)

	p_score <- ddply(bond_counts, .(sample_source, acc_chem_type, don_chem_type), function(df){
		sample_source <- as.character(df$sample_source)
		don_chem_type <- as.character(df$don_chem_type[1])
		acc_chem_type <- as.character(df$acc_chem_type[1])
		p_bond <- df$count / total_bonds
		p_don <- group_counts[
			as.character(group_counts$sample_source) == sample_source &
			as.character(group_counts$chem_type) == don_chem_type, "count"] / total_don
		p_acc <- group_counts[
			as.character(group_counts$sample_source) == sample_source &
			as.character(group_counts$chem_type) == acc_chem_type, "count"] / total_acc
		data.frame(
			p_score = -log( p_bond / (p_don * p_acc) ) )
	})

	p_score$don_chem_type <- don_chem_type_name_linear(droplevels(p_score$don_chem_type))
	p_score$acc_chem_type <- acc_chem_type_name_linear(droplevels(p_score$acc_chem_type))
	p_score <- na.omit(p_score, method="r")




	table_id <- paste(id, "narrow", sep="_")
	save_tables( self,
		p_score, table_id,
		sample_sources, output_dir, output_formats,
		caption=title, caption.placement="top")

	wide_p_score <- dcast(p_score, sample_source + don_chem_type ~ acc_chem_type, value="p_score")
	table_id <- id
	save_tables( self,
		wide_p_score, table_id,
		sample_sources, output_dir, output_formats,
		caption=title, caption.placement="top")



	d_ply(p_score, .(sample_source), function(df){
		new_sample_source <- df$sample_source[1]
		plot_id <- paste(
			id, new_sample_source, sep="_")

		p <- ggplot(data=df) + theme_bw() +
			ggtitle(
				paste(title,"\nSample Source: ", new_sample_source, sep="")) +
			geom_tile(aes(x=don_chem_type, y=acc_chem_type, fill=p_score)) +
			scale_fill_viridis('P(Score)')
		save_plots(self, plot_id, sample_sources, output_dir, output_formats)
	})

}



#buried
# all hbond sites that have good B-Factor, are buried, and do not have
# any exposed partners with good B-Factor
#roup_sele <-"
#REATE TEMPORARY TABLE exposed_neighbor_count (
#	struct_id BLOB,
#	site_id INTEGER,
#	PRIMARY KEY ( struct_id, site_id ));
#
#REATE INDEX IF NOT EXISTS hbonds_struct_id_acc_id ON hbonds ( struct_id, acc_id );
#REATE INDEX IF NOT EXISTS hbonds_struct_id_don_id ON hbonds ( struct_id, don_id );
#
#NSERT INTO exposed_neighbor_count
#ELECT DISTINCT
#	site.struct_id,
#	site.site_id
#ROM
#	hbond_sites AS site,
#	hbond_sites_pdb AS pdb,
#	hbond_site_environment AS site_env,
#	hbond_sites AS alt_site,
#	hbond_sites_pdb AS alt_pdb,
#	hbond_site_environment AS alt_site_env,
#	hbonds AS hbond
#HERE
#	site.struct_id = hbond.struct_id AND
#	alt_site.struct_id = hbond.struct_id AND
#	site.is_donor == 1 AND alt_site.is_donor == 0 AND
#	hbond.don_id == site.site_id AND alt_site.site_id = hbond.acc_id AND
#	pdb.struct_id = site.struct_id AND
#	pdb.site_id = site.site_id AND
#	pdb.heavy_atom_temperature < 30 AND
#	site_env.struct_id = site.struct_id AND
#	site_env.site_id = site.site_id AND
#	site_env.sasa_r140 == 0 AND
#	alt_pdb.struct_id = alt_site.struct_id AND
#	alt_pdb.site_id = alt_site.site_id AND
#	alt_pdb.heavy_atom_temperature < 30 AND
#	alt_site_env.struct_id = alt_site.struct_id AND
#	alt_site_env.site_id = alt_site.site_id AND
#	alt_site_env.sasa_r140 > 0
#ROUP BY
#	site.struct_id, site.site_id;
#
#
#NSERT INTO exposed_neighbor_count
#ELECT DISTINCT
#	site.struct_id,
#	site.site_id
#ROM
#	hbond_sites AS site,
#	hbond_sites_pdb AS pdb,
#	hbond_site_environment AS site_env,
#	hbond_sites AS alt_site,
#	hbond_sites_pdb AS alt_pdb,
#	hbond_site_environment AS alt_site_env,
#	hbonds AS hbond
#HERE
#	site.struct_id = hbond.struct_id AND
#	alt_site.struct_id = hbond.struct_id AND
#	site.is_donor == 0 AND alt_site.is_donor == 1 AND
#	hbond.acc_id == site.site_id AND alt_site.site_id = hbond.don_id AND
#	pdb.struct_id = site.struct_id AND
#	pdb.site_id = site.site_id AND
#	pdb.heavy_atom_temperature < 30 AND
#	site_env.struct_id = site.struct_id AND
#	site_env.site_id = site.site_id AND
#	site_env.sasa_r140 == 0 AND
#	alt_pdb.struct_id = alt_site.struct_id AND
#	alt_pdb.site_id = alt_site.site_id AND
#	alt_pdb.heavy_atom_temperature < 30 AND
#	alt_site_env.struct_id = alt_site.struct_id AND
#	alt_site_env.site_id = alt_site.site_id AND
#	alt_site_env.sasa_r140 > 0
#ROUP BY
#	site.struct_id, site.site_id;
#
#ELECT
#	site.HBChemType AS chem_type,
#	site.is_donor,
#	count( * ) AS count
#ROM
#	hbond_sites AS site,
#	hbond_sites_pdb AS site_pdb,
#	hbond_site_environment AS site_env
#	LEFT JOIN
#		exposed_neighbor_count
#	ON
#		exposed_neighbor_count.struct_id == site.struct_id AND
#		exposed_neighbor_count.site_id == site.site_id
#HERE
#	site_env.struct_id = site.struct_id AND
#	site_env.site_id = site.site_id AND
#	site_env.sasa_r140 == 0 AND
#	site_pdb.struct_id = site.struct_id AND
#	site_pdb.site_id = site.site_id AND
#	site_pdb.heavy_atom_temperature < 30 AND
#	exposed_neighbor_count.struct_id IS NULL
#ROUP BY
#	chem_type;"
#
#ond_sele <-"
#ELECT
# acc_site.HBChemType AS acc_chem_type,
# don_site.HBChemType AS don_chem_type,
# count(*) AS count
#ROM
# hbonds AS hbond,
# hbond_sites AS acc_site,
# hbond_sites AS don_site,
#	hbond_sites_pdb AS don_pdb,
#	hbond_sites_pdb AS acc_pdb,
#	hbond_site_environment AS don_env,
#	hbond_site_environment AS acc_env
#HERE
#	don_site.is_donor == 1 AND acc_site.is_donor == 0 AND
# hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
# hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
#	don_pdb.struct_id = don_site.struct_id AND don_pdb.site_id = don_site.site_id AND
#	acc_pdb.struct_id = acc_site.struct_id AND acc_pdb.site_id = acc_site.site_id AND
#	don_pdb.heavy_atom_temperature < 30 AND
#	acc_pdb.heavy_atom_temperature < 30 AND
#	don_env.struct_id = don_site.struct_id AND don_env.site_id = don_site.site_id AND
#	acc_env.struct_id = acc_site.struct_id AND acc_env.site_id = acc_site.site_id AND
#	don_env.sasa_r140 == 0 AND acc_env.sasa_r140 == 0
#ROUP BY
# acc_site.HBChemType, don_site.HBChemType;"
#o_analysis(
#	group_sele, bond_sele,
#	"hbond_interaction_score_buried", "H-Bond (Buried) Interaction Score; B-Factor < 30")
#

#buried
group_sele <-"
SELECT
  site.HBChemType AS chem_type,
	site.is_donor AS is_donor,
  count(site.HBChemType) AS count
FROM
  hbond_sites AS site,
	hbond_sites_pdb AS site_pdb,
	hbond_site_environment AS site_env
WHERE
	site_pdb.struct_id = site.struct_id AND
	site_pdb.site_id = site.site_id AND
--	site_pdb.heavy_atom_temperature < 30 AND
	site_env.struct_id = site.struct_id AND
	site_env.site_id = site.site_id AND
	site_env.sasa_r140 == 0
GROUP BY
	site.HBChemType;"

bond_sele <-"
SELECT
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type,
  count(*) AS count
FROM
  hbonds AS hbond,
  hbond_sites AS acc_site,
  hbond_sites AS don_site,
	hbond_sites_pdb AS don_pdb,
	hbond_sites_pdb AS acc_pdb,
	hbond_site_environment AS don_env,
	hbond_site_environment AS acc_env
WHERE
  hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
  hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
	don_pdb.struct_id = don_site.struct_id AND don_pdb.site_id = don_site.site_id AND
	acc_pdb.struct_id = acc_site.struct_id AND acc_pdb.site_id = acc_site.site_id AND
--	don_pdb.heavy_atom_temperature < 30 AND
--	acc_pdb.heavy_atom_temperature < 30 AND
	don_env.struct_id = don_site.struct_id AND don_env.site_id = don_site.site_id AND
	acc_env.struct_id = acc_site.struct_id AND acc_env.site_id = acc_site.site_id AND
	don_env.sasa_r140 == 0 AND acc_env.sasa_r140 == 0
GROUP BY
  acc_site.HBChemType, don_site.HBChemType;"
do_analysis(
	group_sele, bond_sele,
	"hbond_interaction_score_buried", "H-Bond (Buried) Interaction Score; B-Factor < 30")


#
##exposed
#group_sele <-"
#SELECT
#  site.HBChemType AS chem_type,
#	site.is_donor AS is_donor,
#  count(site.HBChemType) AS count
#FROM
#  hbond_sites AS site,
#	hbond_sites_pdb AS site_pdb,
#	hbond_site_environment AS site_env
#WHERE
#	site_pdb.struct_id = site.struct_id AND
#	site_pdb.site_id = site.site_id AND
#	site_pdb.heavy_atom_temperature < 30 AND
#	site_env.struct_id = site.struct_id AND
#	site_env.site_id = site.site_id AND
#	site_env.sasa_r100 > 0
#GROUP BY
#	site.HBChemType;"
#
#bond_sele <-"
#SELECT
#  acc_site.HBChemType AS acc_chem_type,
#  don_site.HBChemType AS don_chem_type,
#  count(*) AS count
#FROM
#  hbonds AS hbond,
#  hbond_sites AS acc_site,
#  hbond_sites AS don_site,
#	hbond_sites_pdb AS don_pdb,
#	hbond_sites_pdb AS acc_pdb,
#	hbond_site_environment AS don_env,
#	hbond_site_environment AS acc_env
#WHERE
#  hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
#  hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
#	don_pdb.struct_id = don_site.struct_id AND don_pdb.site_id = don_site.site_id AND
#	acc_pdb.struct_id = acc_site.struct_id AND acc_pdb.site_id = acc_site.site_id AND
#	don_pdb.heavy_atom_temperature < 30 AND
#	acc_pdb.heavy_atom_temperature < 30 AND
#	don_env.struct_id = don_site.struct_id AND don_env.site_id = don_site.site_id AND
#	acc_env.struct_id = acc_site.struct_id AND acc_env.site_id = acc_site.site_id AND
#	don_env.sasa_r100 > 0 AND acc_env.sasa_r100 > 0
#GROUP BY
#  acc_site.HBChemType, don_site.HBChemType;"
#do_analysis(
#	group_sele, bond_sele,
#	"hbond_interaction_score_exposed", "H-Bond (Exposed) Interaction Score; B-Factor < 30")
#
})) # end FeaturesAnalysis
