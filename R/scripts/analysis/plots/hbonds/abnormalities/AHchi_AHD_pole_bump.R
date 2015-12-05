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
source("../hbond_geo_dim_scales.R")

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "AHchi_AHD_pole_bump",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){



sele <-"
SELECT
	geom.cosAHD,
 	acc_atoms.base_x AS abx, acc_atoms.base_y AS aby, acc_atoms.base_z AS abz, -- acceptor base atom
	acc_atoms.atm_x  AS ax,  acc_atoms.atm_y  AS ay,  acc_atoms.atm_z  AS az,  -- acceptor atom
	don_atoms.atm_x  AS hx,  don_atoms.atm_y  AS hy,  don_atoms.atm_z  AS hz,  -- hydrogen atom
	don_atoms.base_x AS dx,  don_atoms.base_y AS dy,  don_atoms.base_z AS dz   -- donor atom
FROM
	hbonds AS hbond cross join
	hbond_sites_pdb AS don_pdb cross join
	hbond_sites_pdb AS acc_pdb cross join
	hbond_sites AS don_site cross join
	hbond_sites AS acc_site cross join
	hbond_site_atoms AS don_atoms cross join
	hbond_site_atoms AS acc_atoms cross join
	hbond_geom_coords AS geom
WHERE
	don_pdb.struct_id = hbond.struct_id AND don_pdb.site_id = hbond.don_id AND
	acc_pdb.struct_id = hbond.struct_id AND acc_pdb.site_id = hbond.acc_id AND
	don_pdb.heavy_atom_temperature < 30 AND acc_pdb.heavy_atom_temperature < 30 AND
	hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
	hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
	ABS(don_site.resNum - acc_site.resNum) > 5 AND
	don_atoms.struct_id = hbond.struct_id AND don_atoms.site_id = hbond.don_id AND
	acc_atoms.struct_id = hbond.struct_id AND acc_atoms.site_id = hbond.acc_id AND
	hbond.struct_id = geom.struct_id AND hbond.hbond_id = geom.hbond_id;"
f <- query_sample_sources(sample_sources, sele)

f <- transform(f,
	AHD = acos(cosAHD),
	AHchi = vector_dihedral(
		cbind(abx, aby, abz), cbind(ax, ay, az),
		cbind(hx, hy, hz), cbind(dx, dy, dz)))

#equal area projection
f <- transform(f,
	capx = 2*sin(acos(cosAHD)/2)*cos(AHchi),
	capy = 2*sin(acos(cosAHD)/2)*sin(AHchi))

zoom_range <- .5
f.zoom <- f[
	f$capx > -zoom_range & f$capx < zoom_range &
	f$capy > -zoom_range & f$capy < zoom_range,]
capx_limits <- c(-zoom_range,zoom_range)
capy_limits <- capx_limits

dens <- ddply(f.zoom, .(sample_source), function(sub_f) {
	dens <- MASS::kde2d(sub_f$capx, sub_f$capy, n=500, h=.01)
	densdf <- data.frame(expand.grid(x=dens$x, y=dens$y), z=as.vector(dens$z)/max(dens$z))
	densdf$counts <- nrow(sub_f)
	densdf
})

facet_labels <- data.frame(
	sample_source = sample_sources$sample_source,
	label=toupper(letters[1:nrow(sample_sources)]))


breaks_map <- function(x){2*sin((180-x)*pi/360)}
break_labels <- c(160,180,160)
breaks <- c(-0.3472964, 0.0000000, 0.3472964)

narrow_output_formats <- transform(output_formats, width=height/2.8)

plot_id <- "AHchi_AHD_pole_bump"
ggplot(data=dens) +
	theme_bw() +
	geom_raster(aes(x=x, y=y, fill=z)) +
#	polar_equal_area_grids_bw(bgcolor="white") +
	geom_indicator(aes(indicator=counts), group=1, color="black", ypos=.99) +
	geom_indicator(
		data=facet_labels,
		aes(indicator=label),
		group=1,
		color="black",
		xpos=.03,
		ypos=.03,
		size=20) +
	coord_equal(ratio=1) +
	scale_fill_gradient('Scaled Density', low="white", high="black") +
	theme(legend.position="bottom", legend.direction="horizontal") +
	scale_x_continuous('Acceptor -- Hydrogen -- Donor Angle (degrees)', limits=capx_limits, breaks=breaks, labels=break_labels) +
	scale_y_continuous('Acceptor -- Hydrogen -- Donor Angle (degrees)', limits=capy_limits, breaks=breaks, labels=break_labels) +
	facet_wrap(~sample_source)
save_plots(self, plot_id, sample_sources, output_dir, narrow_output_formats)


qs <- compute_quantiles(f, c("sample_source"), "AHD")

plot_id = "hbond_AHD_qq"
p <- ggplot(data=qs) + theme_bw() +
	geom_line(aes(y=probs, x=quantiles*180/pi, colour=sample_source)) +
	ggtitle(paste("HBond AHD Angle Quantiles, B-Fact < 30", sep="")) +
	scale_y_continuous("Fraction") +
	scale_x_continuous("Quantiles (degrees)")

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)



qs <- compute_quantiles(f, c("sample_source"), "AHD")

plot_id = "hbond_AHD_cdf"
p <- ggplot(data=qs) + theme_bw() +
	geom_line(aes(y=probs, x=quantiles*180/pi, colour=sample_source)) +
	ggtitle(paste("H-Bond AHD Angle Cumulative Distribution Function, B-Fact < 30", sep="")) +
	scale_y_continuous("Fraction") +
	scale_x_continuous("A-H-D Angle (degrees)")

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)



})) # End FeaturesAnalysis
