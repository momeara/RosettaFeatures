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
id = "chi_BAH_eq_polar_density_PBAtoAMO_by_resolution",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	geom.cosBAH,
	geom.chi,
--	resolution.resolution,
	CASE
		WHEN resolution.resolution < 1.2 THEN 'HIGH'
		WHEN resolution.resolution > 1.8 THEN 'LOW' END AS resolution
FROM
	hbond_geom_coords AS geom,
	hbonds AS hbond,
	hbond_sites AS don_site,
	hbond_sites AS acc_site,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb,
	resolutions AS resolution
WHERE
  don_site.HBChemType == 'hbdon_AMO' AND acc_site.HBChemType == 'hbacc_PBA' AND
	hbond.struct_id = geom.struct_id AND hbond.hbond_id = geom.hbond_id AND
	hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
	hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
	ABS(don_site.resNum - acc_site.resNum) > 5 AND
  don_pdb.struct_id = hbond.struct_id AND don_pdb.site_id = don_site.site_id AND
  acc_pdb.struct_id = hbond.struct_id AND acc_pdb.site_id = acc_site.site_id AND
	don_pdb.heavy_atom_temperature < 30 AND acc_pdb.heavy_atom_temperature < 30 AND
	resolution.struct_id = hbond.struct_id AND
	(resolution.resolution < 1.2 OR resolution.resolution > 1.8);";

f <- query_sample_sources(sample_sources[sample_sources$reference==TRUE,], sele)

f <- ddply(f, c("sample_source", "resolution"),
	transform, counts = length(sample_source))

#equal area projection
f <- transform(f,
	capx = 2*sin(acos(cosBAH)/2)*cos(chi),
	capy = 2*sin(acos(cosBAH)/2)*sin(chi))

capx_limits <- c(-1.5,1.5)
capy_limits <- capx_limits

plot_parts <- list(
	theme_bw(),
	theme(panel.background=element_rect(fill="#00007F", colour="#00007F")),
	stat_density2d(
		aes(x=capx,y=capy, fill=..density..), geom="tile", contour=F),
	geom_point(aes(x=capx, y=capy)),
#	geom_indicator(aes(indicator=counts), color="white"),
	polar_equal_area_grids_bw(),
	coord_equal(ratio=1),
	scale_fill_gradientn('Density', colours=jet.colors(15)),
	scale_x_continuous(limits=capx_limits),
	scale_y_continuous(limits=capy_limits),
	theme(
		axis.text.x=element_blank(),
		axis.text.y=element_blank(),
		axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		axis.ticks.x = element_blank(),
		axis.ticks.y = element_blank()))


d_ply(f, .(sample_source, resolution), function(sub_f){
	ss_id <- sub_f$sample_source[1]
	ss = sample_sources[sample_sources$sample_source == ss_id,]
	sub_f$counts <- nrow(sub_f)

	print(summary(sub_f))

	plot_id = paste("chi_BAH_eq_polar_density_PBAtoAMO_by_resolution_bfact30", ss_id, sub_f$resolution[1], sep="_")
	ggplot(data=sub_f) + plot_parts +
		geom_indicator(aes(indicator=counts), color="white", group=1) +
		ggtitle(paste("Hydrogen Bonds chi vs BAH Angles\nBackbone/Lysine Hydrogen Bonds\nEqual Coordinate Projection   Sample Source: ", ss_id, sep="")) +
#		facet_wrap(~resolution)
	save_plots(self, plot_id, ss, output_dir, output_formats)

})
 

})) # end FeaturesAnalysis
