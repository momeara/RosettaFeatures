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


feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "helical_i_p3",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){



source("../hbonds/hbond_geo_dim_scales.R")

sele <-"
SELECT
	geom.AHdist, geom.cosBAH, geom.cosAHD, geom.chi
FROM
	hbonds AS hb cross join
	hbond_geom_coords AS geom cross join
	hbond_sites AS don cross join
	hbond_sites AS acc cross join
	hbond_sites_pdb AS don_pdb cross join
	hbond_sites_pdb AS acc_pdb cross join
	residue_secondary_structure AS don_ss cross join
	residue_secondary_structure AS acc_ss
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	geom.AHdist < 3.4 AND
	acc_ss.struct_id = acc.struct_id AND acc_ss.resNum = acc.resNum AND
	don_ss.struct_id = don.struct_id AND don_ss.resNum = don.resNum AND
	don_ss.dssp == 'H' AND acc_ss.dssp == 'H' AND
	don.resNum - acc.resNum == 3;"

f <- query_sample_sources(sample_sources, sele)
f <- group_counts(f, id.vars=c("sample_source"))

#### SUMMARY TABLE #####
counts <- ddply(f, .(sample_source), function(df){
  data.frame(counts = nrow(df))
})

table_id <- "Alpha_helix_i_to_ip3_hbond_count"
table_title <- "Counts of a-Helix i->i+3 H-bonds; B-Factor < 30"
save_tables(self,
	counts, table_id,
	sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")


#### AHdist ######
dens <- estimate_density_1d(
	f, c("sample_source"),
	"AHdist", weight_fun = radial_3d_normalization)

plot_id <- "hbond_helical_i_p3_AHdist"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source), size=1.5) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	ggtitle("a-Helix i->i+3 H-Bonds: A-H Distance;\nB-Factor < 30 normalized for equal weight per unit distance") +
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(1.4,3.4), breaks=c(1.4, 1.8, 2.2, 2.8, 3.2))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


#### cosBAH vs AHdist3 ####
plot_id <- "helical_i_p3_cosBAH_AHdist3"
ggplot(data=f, aes(x=cosBAH, y=AHdist^3)) +
	theme_bw() +
	geom_point(size=.4) +
	stat_density2d(size=.2) +
	facet_grid(~sample_source) +
	geom_indicator(aes(indicator=counts), group=1) +
	scale_x_continuous(
		"cos(Base -- Acceptor -- Hydrogen)",
		limit=c(-.6,1), breaks=c(-.6, -.3, 0, .3, .6, .9)) +
	scale_y_continuous(
		expression(paste('(Acceptor -- Hydrogen Distance)^3')),
		limits=c(1.4^3, 3.4^3), breaks=c(1.4^3, 2^3, 3^3), labels=c(1.4, 2, 3)) +
	ggtitle("a-Helix i->i+3 H-Bonds: BAH vs AHdist\nB-Factor < 30 normalized for equal weight per unit distance")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

##### cosAHD vs AHdist3 ####
plot_id <- "helical_i_p3_cosAHD_AHdist3"
f <- group_counts(f, id.vars=c("sample_source"))
ggplot(data=f, aes(x=cosAHD, y=AHdist^3)) +
  	theme_bw() + 
	geom_point(size=.4) +
	stat_density2d(size=.2) +
	facet_grid(~sample_source) +
	geom_indicator(aes(indicator=counts), group=1) +
	scale_x_continuous(
		"cos(Acceptor -- Hydrogen -- Donor)",
		limit=c(0,1), breaks=c(.2, .4, .6, .8, 1)) +
	scale_y_continuous(
		expression(paste('(Acceptor -- Hydrogen Distance)^3')),
		limits=c(1.4^3, 3.4^3), breaks=c(1.4^3, 2^3, 3^3), labels=c(1.4, 2, 3)) +	
	ggtitle("a-Helix i->i+3 H-Bonds: AHD vs AHdist\nB-Factor < 30 normalized for equal weight per unit distance")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

#### BAH vs chi #####
#equal area projection
f <- transform(f,
	capx = 2*sin(acos(cosBAH)/2)*cos(chi),
	capy = 2*sin(acos(cosBAH)/2)*sin(chi))
capx_limits <- c(-1.5,1.5)
capy_limits <- capx_limits

plot_id <- "helical_i_p3_BAH_chi"
ggplot(data=f, aes(x=capx, y=capy)) +
	theme_bw() +
	geom_point(size=.4) +
	stat_density2d(size=.2) +
	geom_indicator(aes(indicator=counts), group=1) +
	polar_equal_area_grids_bw(line_color="gray50", box_bgcolor="gray50") +
	ggtitle("a-Helix i->i+3 H-Bonds: BAH vs chi\nB-Factor < 30 normalized for equal weight per unit distance") +
	scale_x_continuous('', limits=capx_limits, breaks=c()) +
	scale_y_continuous('', limits=capy_limits, breaks=c()) +
	facet_wrap( ~ sample_source, ncol=ceiling(sqrt(nrow(sample_sources)))) +
	coord_equal(ratio=1)
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
