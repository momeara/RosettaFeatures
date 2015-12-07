# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

library(rgl)
library(misc3d)
library(mvtnorm)
library(ks)
library(plyr)
source("../hbond_geo_dim_scales.R")

feature_analyses <- c(feature_analyses, methods::new("FeaturesAnalysis",
id = "AHdist_sp2_dPBA_aPBA",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
	geom.AHdist,
	geom.chi,
	acc_atoms.base_x  AS bx,  acc_atoms.base_y  AS by,  acc_atoms.base_z  AS bz,
	acc_atoms.base2_x AS b2x, acc_atoms.base2_y AS b2y, acc_atoms.base2_z AS b2z,
	acc_atoms.atm_x   AS ax,  acc_atoms.atm_y   AS ay,  acc_atoms.atm_z   AS az,
	don_atoms.atm_x   AS hx,  don_atoms.atm_y   AS hy,  don_atoms.atm_z   AS hz,
	acc_site.HBChemType AS acc_chem_type,
	don_site.HBChemType AS don_chem_type,
	CASE acc_site.HBChemType
		WHEN 'hbacc_IMD' THEN 'ring' WHEN 'hbacc_IME' THEN 'ring'
		WHEN 'hbacc_AHX' THEN 'sp3'  WHEN 'hbacc_HXL' THEN 'sp3'
		WHEN 'hbacc_CXA' THEN 'sp2'  WHEN 'hbacc_CXL' THEN 'sp2'
		WHEN 'hbacc_PBA' THEN 'bb_sp2'  END AS hybrid
FROM
	hbond_geom_coords AS geom,
	hbonds AS hb,
	hbond_sites AS don_site,
	hbond_sites AS acc_site,
	hbond_site_atoms AS don_atoms,
	hbond_site_atoms AS acc_atoms
WHERE
	hb.struct_id = geom.struct_id AND hb.hbond_id = geom.hbond_id AND
	hb.struct_id = don_site.struct_id AND hb.don_id = don_site.site_id AND
	hb.struct_id = acc_site.struct_id AND hb.acc_id = acc_site.site_id AND
	don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
	acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id AND
	ABS(don_site.resNum - acc_site.resNum) > 5;"

f <- query_sample_sources(sample_sources, sele)

coords <- f %>%
	dplyr::mutate(
		don_chem_type_name = don_chem_type_name_linear(don_chem_type),
		acc_chem_type_name = acc_chem_type_name_linear(acc_chem_type),
		cosBAH = ifelse(hybrid %in% c("sp3", "ring"),
			vector_dotprod(
				vector_normalize(cbind(ax-(bx+b2x)/2, ay-(by+b2y)/2, az-(bz+b2z)/2)),
				vector_normalize(cbind(hx-ax, hy-ay, hz-az))),
			vector_dotprod(
				vector_normalize(cbind(ax-bx, ay-by, az-bz)),
				vector_normalize(cbind(hx-ax, hy-ay, hz-az)))),
		bah = acos(bah),
		chi = ifelse(hybrid %in% c("sp3", "ring"),
			vector_dihedral(
				cbind(b2x, b2y, b2z),
				cbind((bx+b2x)/2, (by+b2y)/2, (bz+b2z)/2),
				cbind(ax, ay, az),
				cbind(hx, hy, hz)),
			vector_dihedral(
				cbind(b2x, b2y, b2z),
				cbind(bx, by, bz),
				cbind(ax, ay, az),
				cbind(hx, hy, hz))),
		x = AHdist * sin(bah) * cos(chi),
		y = AHdist * sin(bah) * sin(chi),
		z = AHdist * cos(bah)) %>%
	dplyr::select(
		sample_source,
		don_chem_type_name,
		acc_chem_type_name,
		x, y, z)

#save memory
f <- NULL
gc()

make_acc_plot_3d <- function(dens, full_path, n=200, alo=0.1, ahi=0.5) {
	cat("Saving plot: ", full_path)
	timing <- system.time({
		plot(dens, xlab="", ylab="", zlab="", box=FALSE, axes=FALSE, add=FALSE)

		lim <- sqrt(r_max^2/3)
		th <- c(-1.2, -1.0, -.78 )
		col <- c("#00007F", "#7F0000", "#7FFF7F" )
		#  col <- rev(cmap(length(th)))
		al <- seq(alo, ahi, len=length(th))
		x <- seq(-lim, lim, len=n)

		atm <- function(a,b,c, color, m=0, s=0.2) {
			atm_fn <- function(x,y,z) {
				dnorm( x-a, m, s) * dnorm(y-b, m, s) * dnorm(z-c, m, s)
			}
			contour3d(atm_fn, c(1), x,x,x, color=color, smooth=5, engine="rgl", add=TRUE)
		}

		bond <- function(vx, vy, vz, wx, wy, wz, color, width=.1) {
			bond_fn <- function(x,y,z) {
				p <- cbind(x,y,z)
				v <- rep(c(vx, vy, vz), each=length(x))
				w <- rep(c(wx, wy, wz), each=length(x))
				dim(v) <- c(length(x), 3)
				dim(w) <- c(length(x), 3)
				norm <- function(a) {sqrt(a[,1]*a[,1] + a[,2]*a[,2] + a[,3]*a[,3])}
				dot <- function(a,b) { a[,1] * b[,1] + a[,2] * b[,2] + a[,3] * b[,3]}
				dist_to_vw <- sin(acos(dot((v-w)/norm(v-w), (p-w)/norm(p-w)))) * norm(p-w)
				ifelse(
					(dist_to_vw > width)  | (dot(p-v,w-v) < 0) | (dot(p-w,v-w) < 0),
					0, 1)
			}
			contour3d(bond_fn, c(.01), x,x,x, color=color, smooth=5, engine="rgl", add=TRUE)
		}

		a_atm <- atm(0, 0, 0, "red")
		b_atm <- atm(0, 0, 2, "gray")
		bb_atm <- atm(sqrt(2), 0, sqrt(2)+2, "gray")

		a_b_bond <- bond(0,0,0, 0,0,2, "blue")
		b_bb_bond <- bond(0,0,2, sqrt(2), 0, sqrt(2)+2, "blue")

		rgl.bg(col="white")
		rgl.viewpoint(theta=0,phi=60,fov=30,zoom=1)
		rgl.snapshot(full_path, fmt="png")

		rgl.viewpoint(theta=170,phi=-40,fov=30,zoom=1)
		rgl.snapshot(full_path, fmt="png")

	})
	cat(" ... ", as.character(round(timing[3],2)), "s\n", sep="")
}

ddply(coords, .(sample_source, don_chem_type, acc_chem_type), function(sub_coords){
	dm <- ks::kde(as.matrix(sub_coords[,c("x", "y", "z")]))
	plot_id <- paste("acc_contour_contour", sub_coords$sample_source[1], sep="_")
	full_output_dir <- file.path(output_dir, self@id, "output_print_raster")
	if(!file.exists(full_output_dir)){
		dir.create(full_output_dir, recursive=TRUE)
	}
	full_path <- file.path(full_output_dir, paste(plot_id, ".png", sep=""))

	make_acc_plot_3d(dm, full_path)
})


})) # end FeaturesAnalysis
