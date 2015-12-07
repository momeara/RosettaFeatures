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
library(viridis)

feature_analyses <- c(feature_analyses, methods::new("FeaturesAnalysis",
id = "VBAH_VBchi_equal_area",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
	acc_atoms.base_x  AS abx,  acc_atoms.base_y  AS aby,  acc_atoms.base_z  AS abz,
  acc_atoms.base2_x AS ab2x, acc_atoms.base2_y AS ab2y, acc_atoms.base2_z AS ab2z,
  acc_atoms.atm_x   AS ax,  acc_atoms.atm_y   AS ay,  acc_atoms.atm_z   AS az,
  don_atoms.atm_x   AS hx,  don_atoms.atm_y   AS hy,  don_atoms.atm_z   AS hz,
  CASE don.resNum - don.resNum
    WHEN -1 THEN '-1' WHEN -2 THEN '-2' WHEN -3 THEN '-3' WHEN -4 THEN '-4'
    WHEN 1 THEN '1,2,3,4' WHEN 2 THEN '1,2,3,4' WHEN 3 THEN '1,2,3,4' WHEN 4 THEN '1,2,3,4'
    ELSE 'long' END AS seq_sep,
	acc.HBChemType AS acc_chem_type,
	don.HBChemType AS don_chem_type,
	don_ss_code.label AS don_dssp
FROM
  hbonds AS hb cross join
  hbond_sites AS don cross join
  hbond_sites AS acc cross join
	hbond_sites_pdb AS don_pdb cross join
	hbond_sites_pdb AS acc_pdb cross join
  hbond_site_atoms AS don_atoms cross join
  hbond_site_atoms AS acc_atoms cross join
	residue_secondary_structure AS don_ss cross join
	dssp_codes AS don_ss_code
WHERE
  acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
  don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
  (acc.HBChemType = 'hbacc_HXL' OR acc.HBChemType = 'hbacc_AHX') AND
  don.HBChemType = 'hbdon_PBA' AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
  don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
  acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id AND
	don_ss.struct_id = don.struct_id AND don_ss.resNum = don.resNum AND
	don_ss_code.code = don_ss.dssp;"

f <- query_sample_sources(sample_sources, sele)

f$acc_chem_type_name <- factor(f$acc_chem_type,
	levels = c("hbacc_HXL", "hbacc_AHX"),
	labels = c("aHXL: s,t", "aAHX: y"))

f <- transform(f,
	vBAH = acos(vector_dotprod(
    vector_normalize(cbind(ax-(abx+ab2x)/2, ay-(aby+ab2y)/2, az-(abz+ab2z)/2)),
		vector_normalize(cbind(hx-ax, hy-ay, hz-az)))),
	vBAchi = vector_dihedral(
		cbind(ab2x, ab2y, ab2z), cbind((abx+ab2x)/2, (aby+ab2y)/2, (abz+ab2z)/2),
		cbind(ax, ay, az), cbind(hx, hy, hz)))



#equal area projection
f <- transform(f,
	x = 2*sin(vBAH/2)*cos(vBAchi),
	y = 2*sin(vBAH/2)*sin(vBAchi))

capx_limits <- c(-1.5,1.5); capy_limits <- capx_limits;

narrow_output_formats <- transform(output_formats, width=height)

plot_parts <- function(){list(
	theme_bw(),
	theme(panel.background=element_rect(fill="#00007F", colour="#00007F")),
	geom_raster(aes(x=x, y=y, fill=z)),
	geom_indicator(aes(indicator=counts), color="white", group=1),
	polar_equal_area_grids_bw(),
	coord_equal(ratio=1),
	scale_fill_viridis("Density"),
	scale_x_continuous('', limits=capx_limits, breaks=c()),
	scale_y_continuous('', limits=capy_limits, breaks=c()),
	theme(
		axis.text.x=element_blank(),
		axis.text.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.ticks.x = element_blank(),
		axis.ticks.y = element_blank()))
}

#plot_id = "hbond_VBAH_VBAchi_equal_area_log_scale_seq_sep"
#ggplot(data=f) + plot_parts() +
#	facet_grid(seq_sep ~ sample_source) +
#	ggtitle("Sp3 Acc H-Bonds Virtual-Base Acceptor Hydrogen projection\nEqual Coordinate Projection") +
#save_plots(self, plot_id, sample_sources, output_dir, narrow_output_formats)
#
d_ply(
	.data=f[f$don_chem_type == "hbdon_PBA", ],
	.variables=.(seq_sep),
	.fun=function(df)
{

	seq_sep <- as.character(df$seq_sep[1])

	dens <- estimate_density_2d(
		df, c("sample_source"), "x", "y")


	print(summary(dens))

	plot_id = paste("hbond_VBAH_VBAchi_equal_area_log_scale_seq_sep_", seq_sep, sep="")
	ggplot(data=dens) + plot_parts() +
		facet_wrap( ~ sample_source, nrow=ceiling(sqrt(nrow(sample_sources)))) +
		ggtitle(paste("Sp3 Acc H-Bonds Virtual-Base Acceptor Hydrogen projection\nEqual Coordinate Projection   Sequence Separation: ", seq_sep, sep="")) +
	save_plots(self, plot_id, sample_sources, output_dir, narrow_output_formats)

})

d_ply(
	.data=f[f$don_chem_type == "hbdon_PBA", ],
	.variables=.(seq_sep, acc_chem_type_name, sample_source),
	.fun=function(df){

	seq_sep <- as.character(df$seq_sep[1])
	acc_chem_type_name <- as.character(df$acc_chem_type_name[1])
	ss_id <- as.character(df$sample_source[1])
	ss <- sample_sources[sample_sources$sample_source == ss_id, ]

	dens <- estimate_density_2d(
		df, c(), "x", "y")

	plot_id = paste("hbond_VBAH_VBAchi_equal_area_log_scale_seq_sep_", seq_sep, "_", acc_chem_type_name, "_", ss_id, sep="")
	ggplot(data=dens) + plot_parts() +
		ggtitle(paste(acc_chem_type_name,  " Acc BB Don H-Bonds Base Acceptor Hydrogen projection\nSeqSep: ", seq_sep, " SS: ", ss_id, sep="")) +
	save_plots(self, plot_id, sample_sources, output_dir, narrow_output_formats)
})


#
#
#d_ply(
#	.data=f[f$don_chem_type == "hbdon_PBA", ],
#	.variables=.(sample_source),
#	.fun=function(df){
#	ss_id <- as.character(df$sample_source[1])
#	ss <- sample_sources[sample_sources$sample_source == ss_id, ]
#	plot_id = paste("hbond_VBAH_VBAchi_equal_area_log_scale_seq_sep_ss_id_", ss_id, sep="")
#	counts <- group_counts(df, c("seq_sep"))
#	ggplot(data=df) + plot_parts() +
#		facet_wrap( ~ seq_sep, ncol=3) +
#		ggtitle(paste("Sp3 Acc H-Bonds Virtual-Base Acceptor Hydrogen projection\nEqual Coordinate Projection By Sequence Separation  Sample Source: ", ss_id, sep="")) +
#	save_plots(self, plot_id, ss, output_dir, output_formats)
#
#})
#
#plot_id = "hbond_VBAH_VBAchi_equal_area_log_scale_don_dssp"
#ggplot(data=f) + plot_parts() +
#	facet_grid( don_dssp ~ sample_source ) +
#	ggtitle("Sp3 Acc - BB Don H-Bonds Virtual-Base Acceptor Hydrogen projection\nEqual Coordinate Projection") +
#save_plots(self, plot_id, sample_sources, output_dir, narrow_output_formats)
#
#d_ply(
#	.data=f[as.character(f$don_dssp) == "hbdon_PBA",],
#	.variables=.(don_dssp),
#	.fun=function(df){
#	don_dssp <- as.character(df$don_dssp[1])
#	plot_id = paste("hbond_VBAH_VBAchi_equal_area_log_scale_don_dssp_", don_dssp, sep="")
#	ggplot(data=df) + plot_parts() +
#		facet_wrap( ~ sample_source, nrow=ceiling(sqrt(nrow(sample_sources)))) +
#		ggtitle(paste("Sp3 Acc - BB Don H-Bonds Virtual-Base Acceptor Hydrogen projection\nEqual Coordinate Projection   Donor DSSP: ", don_dssp, sep="")) +
#	save_plots(self, plot_id, sample_sources, output_dir, narrow_output_formats)
#})
#
#d_ply(
#	.data=f[as.character(f$don_dssp) == "hbdon_PBA",],
#	.variables=.(sample_source),
#	.fun=function(df){
#	ss_id <- as.character(df$sample_source[1])
#	ss <- sample_sources[sample_sources$sample_source == ss_id, ]
#	plot_id = paste("hbond_VBAH_VBAchi_equal_area_log_scale_don_dssp_ss_id_", ss_id, sep="")
#	ggplot(data=df) + plot_parts() +
#		facet_wrap( ~ don_dssp) +
#		ggtitle(paste("Sp3 Acc - BB Don H-Bonds Virtual-Base Acceptor Hydrogen projection\nEqual Coordinate Projection  by Donor DSSP  Sample Source: ", ss_id, sep="")) +
#	save_plots(self, plot_id, ss, output_dir, output_formats)
#})



})) # end FeaturesAnalysis
