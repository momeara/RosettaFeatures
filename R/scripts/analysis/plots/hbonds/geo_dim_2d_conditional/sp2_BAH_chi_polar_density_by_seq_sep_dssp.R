# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "sp2_BAH_chi_polar_density_by_seq_sep_dssp",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("../hbond_geo_dim_scales.R")

sele <-"
SELECT
	acc_atoms.base_x  AS abx,  acc_atoms.base_y  AS aby,  acc_atoms.base_z  AS abz,
  acc_atoms.base2_x AS ab2x, acc_atoms.base2_y AS ab2y, acc_atoms.base2_z AS ab2z,
  acc_atoms.atm_x   AS ax,  acc_atoms.atm_y   AS ay,  acc_atoms.atm_z   AS az,
  don_atoms.atm_x   AS hx,  don_atoms.atm_y   AS hy,  don_atoms.atm_z   AS hz,
  CASE don.resNum - acc.resNum
    WHEN -1 THEN '-1' WHEN -2 THEN '-2' WHEN -3 THEN '-3' WHEN -4 THEN '-4'
    WHEN 1 THEN '1,2,3,4' WHEN 2 THEN '1,2,3,4' WHEN 3 THEN '1,2,3,4' WHEN 4 THEN '1,2,3,4'
    ELSE 'long' END AS seq_sep,
	acc.HBChemType AS acc_chem_type,
	don.HBChemType AS don_chem_type,
	acc_ss_code.label AS acc_ss
FROM
  hbonds AS hb cross join
  hbond_sites AS don cross join
  hbond_sites AS acc cross join
	hbond_sites_pdb AS don_pdb cross join
	hbond_sites_pdb AS acc_pdb cross join
  hbond_site_atoms AS don_atoms cross join
  hbond_site_atoms AS acc_atoms cross join
	residue_secondary_structure AS acc_ss cross join
	dssp_codes AS acc_ss_code
WHERE
  acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
  don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
  (acc.HBChemType = 'hbacc_PBA' OR
	acc.HBChemType = 'hbacc_CXA' OR
 	acc.HBChemType = 'hbacc_CXL') AND
	don.HBChemType != 'hbdon_NONE' AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
  don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
  acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id AND
	acc_ss.struct_id = acc.struct_id AND acc_ss.resNum = acc.resNum AND
	acc_ss_code.code = acc_ss.dssp;"

f <- query_sample_sources(sample_sources, sele)

f <- transform(f,
	acc_chem_type_name = factor(acc_chem_type,
	levels = c("hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
	labels = c("aCXA: n,q", "aCXL: d,e", "aPBA: bb")),
	acc_chem_type_name = factor(
  	acc_chem_type,
		levels = c("hbacc_PBA", "hbacc_CXA", "hbacc_CXL"),
		labels = c("aPBA: bb", "aCXA: n,q", "aCXL: d,e")),
	acc_ss_wrap = factor(
		acc_ss,
		levels = c(
			"H: a-Helix", "E: b-Sheet", "T: HB Turn",
			"G: 3/10 Helix", "B: b-Bridge", "S: Bend",
			"I: pi-Helix", "Irregular"),
		labels = c(
			"H: a-Helix", "E: b-Sheet", "T: HB Turn",
			"G: 3/10 Helix", "B: b-Bridge", "S: Bend",
			"I: pi-Helix", "Irregular")),
	acc_ss_linear = factor(
		acc_ss,
		levels = c(
			"H: a-Helix", "G: 3/10 Helix", "I: pi-Helix",
			"E: b-Sheet", "B: b-Bridge",
			"T: HB Turn", "S: Bend", "Irregular"),
		labels = c(
			"H: a-Helix", "G: 3/10 Helix", "I: pi-Helix",
			"E: b-Sheet", "B: b-Bridge",
			"T: HB Turn", "S: Bend", "Irregular")))


f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type)
f <- na.omit(f, method="r")

f <- transform(f,
	BAH = acos(vector_dotprod(
    vector_normalize(cbind(ax-abx, ay-aby, az-abz)),
		vector_normalize(cbind(hx-ax, hy-ay, hz-az)))),
	BAchi = vector_dihedral(
		cbind(ab2x, ab2y, ab2z), cbind(abx, aby, abz),
		cbind(ax, ay, az), cbind(hx, hy, hz)))

#equal area projection
f <- transform(f,
	capx = 2*sin(BAH/2)*cos(BAchi),
	capy = 2*sin(BAH/2)*sin(BAchi))


f <- subset(
	f,
	select=c(
		"sample_source",
		"seq_sep",
		"acc_chem_type", "acc_chem_type_name",
		"don_chem_type", "don_chem_type_name",
		"acc_ss_linear", "acc_ss_wrap",
		"capx", "capy"))

n_pts <- 100

capx_limits <- c(-1.5,1.5); capy_limits <- capx_limits;


narrow_output_formats <- transform(output_formats, width=height/1.5)

plot_parts_grid_log <- function(){
	list(
		theme_bw(),
		theme(panel.background=element_rect(fill="#00007F", colour="#00007F")),
		geom_raster(aes(x=x, y=y, fill=log(z+1))),
		polar_equal_area_grids_bw(scale=.3, label_scale=.5),
		geom_indicator(aes(indicator=counts), color="white", group=1),
		scale_x_continuous('', limits=capx_limits, breaks=c()),
		scale_y_continuous('', limits=capy_limits, breaks=c()),
		scale_fill_gradientn('log(Density + 1)', colours=jet.colors(15)),
		coord_equal(ratio=1),
		theme(
			axis.text.x=element_blank(),
			axis.text.y=element_blank(),
			axis.title.x=element_blank(),
			axis.title.y=element_blank(),
			axis.ticks.x = element_blank(),
			axis.ticks.y = element_blank()))

}

plot_parts <- function(){
	list(
		theme_bw(),
		theme(panel.background=element_rect(fill="#00007F", colour="#00007F")),
		geom_raster(aes(x=x, y=y, fill=log(z+1))),
		polar_equal_area_grids_bw(),
		geom_indicator(aes(indicator=counts), color="white", group=1),
		scale_x_continuous('', limits=capx_limits, breaks=c()),
		scale_y_continuous('', limits=capy_limits, breaks=c()),
		scale_fill_gradientn('Density', colours=jet.colors(15)),
		coord_equal(ratio=1),
		theme(
			axis.text.x=element_blank(),
			axis.text.y=element_blank(),
			axis.title.x=element_blank(),
			axis.title.y=element_blank(),
			axis.ticks.x = element_blank(),
			axis.ticks.y = element_blank()))

}

est_dens <- function(z, ids, n_pts=400) {
	estimate_density_2d(z, ids, "capx", "capy", n_pts=n_pts, verbose=T, scaled=T)
}

plot_id = "sp2_hbond_BAH_BAchi_equal_area_log_scale_bb_bb_by_acc_ss"
sub_f <- f[
	f$don_chem_type == "hbdon_PBA" & f$acc_chem_type == "hbacc_PBA",
	c("sample_source", "acc_ss_linear", "capx", "capy")]
dens <- est_dens(sub_f, c("acc_ss_linear", "sample_source"))
ggplot(data=dens) + plot_parts_grid_log() +
	facet_grid(acc_ss_linear ~ sample_source) +
	ggtitle("BB-BB H-Bonds Base-Acceptor-Hydrogen Scaled Eq-Area Projection by Don-DSSP")
narrow_output_formats <- transform(output_formats, width=height/1.5)
save_plots(self, plot_id, sample_sources, output_dir, narrow_output_formats)

plot_id = "sp2_hbond_BAH_BAchi_equal_area_log_scale_cxl_bb_by_seq_sep"
sub_f <- f[
	f$don_chem_type == "hbdon_PBA" & f$acc_chem_type == "hbacc_CXL",
	c("sample_source", "seq_sep", "capx", "capy")]
dens <- est_dens(sub_f, c("seq_sep", "sample_source"))
ggplot(data=dens) + plot_parts_grid_log() +
	facet_grid(seq_sep ~ sample_source) +
	ggtitle("ASP/GLU-BB H-Bonds Base-Acceptor-Hydrogen Scaled Eq-Area Projection")
narrow_output_formats <- transform(output_formats, width=height/1.5)
save_plots(self, plot_id, sample_sources, output_dir, narrow_output_formats)

plot_id = "sp2_hbond_BAH_BAchi_equal_area_log_scale_cxa_bb_by_seq_sep"
sub_f <- f[
	f$don_chem_type == "hbdon_PBA" & f$acc_chem_type == "hbacc_CXA",
	c("sample_source", "seq_sep", "capx", "capy")]
dens <- est_dens(sub_f, c("seq_sep", "sample_source"))
ggplot(data=dens) + plot_parts_grid_log() +
	facet_grid(seq_sep ~ sample_source) +
	ggtitle("ASN/GLN-BB H-Bonds Base-Acceptor-Hydrogen Scaled Eq-Area Projection")
narrow_output_formats <- transform(output_formats, width=height/1.5)
save_plots(self, plot_id, sample_sources, output_dir, narrow_output_formats)

plot_id = "sp2_hbond_BAH_BAchi_equal_area_log_scale_bb_sc_by_don_chem_type"
sub_f <- f[
	f$don_chem_type != "hbdon_PBA" & f$acc_chem_type == "hbacc_PBA",
	c("sample_source", "don_chem_type_name", "capx", "capy")]
dens <- est_dens(sub_f, c("don_chem_type_name", "sample_source"))
ggplot(data=dens) + plot_parts_grid_log() +
	facet_grid(don_chem_type_name ~ sample_source) +
	ggtitle("aBB-dSC H-Bonds Base-Acceptor-Hydrogen Scaled Eq-Area Projection")
narrow_output_formats <- transform(output_formats, width=height/2)
save_plots(self, plot_id, sample_sources, output_dir, narrow_output_formats)

plot_id = "sp2_hbond_BAH_BAchi_equal_area_log_scale_cxl_sc_by_don_chem_type"
sub_f <- f[
	f$don_chem_type != "hbdon_PBA" & f$acc_chem_type == "hbacc_CXL",
	c("sample_source", "don_chem_type_name", "capx", "capy")]
dens <- est_dens(sub_f, c("don_chem_type_name", "sample_source"))
ggplot(data=dens) + plot_parts_grid_log() +
	facet_grid(don_chem_type_name ~ sample_source) +
	ggtitle("aCXL-dSC H-Bonds Base-Acceptor-Hydrogen Scaled Eq-Area Projection")
narrow_output_formats <- transform(output_formats, width=height/2)
save_plots(self, plot_id, sample_sources, output_dir, narrow_output_formats)

plot_id = "sp2_hbond_BAH_BAchi_equal_area_log_scale_cxa_sc_by_don_chem_type"
sub_f <- f[
	f$don_chem_type != "hbdon_PBA" & f$acc_chem_type == "hbacc_CXA",
	c("sample_source", "don_chem_type_name", "capx", "capy")]
dens <- est_dens(sub_f, c("don_chem_type_name", "sample_source"))
ggplot(data=dens) + plot_parts_grid_log() +
	facet_grid(don_chem_type_name ~ sample_source) +
	ggtitle("aCXA-dSC H-Bonds Base-Acceptor-Hydrogen Scaled Eq-Area Projection")
narrow_output_formats <- transform(output_formats, width=height/2)
save_plots(self, plot_id, sample_sources, output_dir, narrow_output_formats)








plot_id = "sp2_hbond_BAH_BAchi_equal_area_log_scale_seq_sep"
dens <- est_dens(f[f$acc_chem_type != "hbacc_PBA",], c("seq_sep", "sample_source"))
ggplot(data=dens) + plot_parts_grid_log() +
	facet_grid(seq_sep ~ sample_source) +
	ggtitle("Sp2 Sc Acc H-Bonds Base Acceptor Hydrogen Scaled EQ-Area Projection") +
save_plots(self, plot_id, sample_sources, output_dir, narrow_output_formats)



d_ply(
	.data=f[f$acc_chem_type != "hbacc_PBA" & f$don_chem_type == "hbdon_PBA", ],
	.variables=.(seq_sep, acc_chem_type_name, sample_source),
	.parallel=use_parallel,
	.fun=function(df){
	seq_sep <- as.character(df$seq_sep[1])
	acc_chem_type_name <- as.character(df$acc_chem_type_name[1])
	ss_id <- as.character(df$sample_source[1])
	ss <- sample_sources[sample_sources$sample_source == ss_id, ]
	plot_id = paste("hbond_BAH_BAchi_equal_area_log_scale_seq_sep_", seq_sep, "_", acc_chem_type_name, "_", ss_id, sep="")
	dens <- est_dens(df,	c())
	ggplot(data=dens) + plot_parts() +
		ggtitle(paste(acc_chem_type_name,  " Acc BB Don H-Bonds Base Acceptor Hydrogen projection\nSeqSep: ", seq_sep, " SS: ", ss_id, sep="")) +
	save_plots(self, plot_id, sample_sources, output_dir, narrow_output_formats)
})

d_ply(
	.data=f,
	.variables=.(sample_source),
	.parallel=use_parallel,
	.fun=function(df){
	ss_id <- as.character(df$sample_source[1])
	ss <- sample_sources[sample_sources$sample_source == ss_id, ]
	plot_id = paste("hbond_VBAH_VBAchi_equal_area_log_scale_seq_sep_ss_id_", ss_id, sep="")
	dens <- est_dens(df,	c("seq_sep"))
	ggplot(data=dens) + plot_parts() +
		facet_wrap( ~ seq_sep, ncol=3) +
		ggtitle(paste("Sp2 Acc H-Bonds Base Acceptor Hydrogen projection\nEqual Coordinate Projection By Sequence Separation  Sample Source: ", ss_id, sep="")) +
	save_plots(self, plot_id, ss, output_dir, output_formats)
})

plot_id = "hbond_VBAH_VBAchi_equal_area_log_scale_acc_ss"
dens <- est_dens(df,	c("acc_ss_linear", "sample_source"))
ggplot(data=dens) + plot_parts_grid_log() +
	facet_grid( acc_ss_linear ~ sample_source ) +
	ggtitle("Sp2 Acc - BB Don H-Bonds Base Acceptor Hydrogen projection\nEqual Coordinate Projection") +
save_plots(self, plot_id, sample_sources, output_dir, narrow_output_formats)

d_ply(
	.data=f[as.character(f$don_chem_type) == "hbdon_PBA",],
	.variables=.(acc_ss_linear),
	.parallel=use_parallel,
	.fun=function(df){
	acc_ss_linear <- as.character(df$acc_ss_linear[1])
	plot_id = paste("hbond_VBAH_VBAchi_equal_area_log_scale_acc_ss_", acc_ss_linear, sep="")
		dens <- est_dens(df,	c("sample_source"))
	ggplot(data=dens) + plot_parts() +
		facet_wrap( ~ sample_source, nrow=ceiling(sqrt(nrow(sample_sources)))) +
		ggtitle(paste("Sp2 Acc - BB Don H-Bonds Base Acceptor Hydrogen projection\nEqual Coordinate Projection   Donor DSSP: ", acc_ss_linear, sep="")) +
	save_plots(self, plot_id, sample_sources, output_dir, narrow_output_formats)
})

d_ply(
	.data=f[as.character(f$don_chem_type) == "hbdon_PBA",],
	.variables=.(sample_source),
	.parallel=use_parallel,
	.fun=function(df){
	ss_id <- as.character(df$sample_source[1])
	ss <- sample_sources[sample_sources$sample_source == ss_id, ]
	plot_id = paste("hbond_VBAH_VBAchi_equal_area_log_scale_acc_ss_ss_id_", ss_id, sep="")
	dens <- est_dens(df,	c("acc_ss_wrap"))
	ggplot(data=dens) + plot_parts_grid_log() +
		facet_wrap( ~ acc_ss_wrap) +
		ggtitle(paste("Sp2 Acc - BB Don H-Bonds Base Acceptor Hydrogen projection\nEqual Coordinate Projection  by Donor DSSP  Sample Source: ", ss_id, sep="")) +
	save_plots(self, plot_id, ss, output_dir, output_formats)
})




d_ply(
	.data=f[f$don_chem_type == "hbdon_PBA",],
	.variables=.(sample_source),
	.parallel=use_parallel,
	.fun=function(df){
	ss_id <- as.character(df$sample_source[1])
	plot_id = paste("hbond_VBAH_VBAchi_equal_area_log_scale_bb_don_ss_id_", ss_id, sep="")
	dens <- est_dens(df,	c("acc_ss_linear", "acc_chem_type_name"))
	ggplot(data=dens) + plot_parts_grid_log() +
		facet_grid(acc_ss_linear ~ acc_chem_type_name) +
		ggtitle(paste("HBond ", ss_id, " to BB Don by Acc DSSP", sep=""))
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})

d_ply(
	.data=f[f$don_chem_type != "hbdon_PBA",],
	.variables=.(sample_source),
	.parallel=use_parallel,
	.fun=function(df){
	ss_id <- as.character(df$sample_source[1])
	plot_id = paste("hbond_VBAH_VBAchi_equal_area_log_scale_sc_don_don_chem_type_ss_id_", ss_id, sep="")
	dens <- est_dens(df,	c("don_chem_type_name", "acc_chem_type_name"))
	ggplot(data=dens) + plot_parts_grid_log() +
		facet_grid(don_chem_type_name ~ acc_chem_type_name) +
		ggtitle(paste("HBond ", ss_id, " By Don Chem Type ", sep=""))
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})



})) # end FeaturesAnalysis
