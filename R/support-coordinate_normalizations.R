# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

radial_3d_normalization <- function(x){ 1/(x^2*sum(1/x^2))}

conical_3d_normalization <- function(x){(1/sin(x))*1/sum(1/(sin(x)))}

uniform_normalization <- function(x){ rep(1/length(x),length(x)) }

# normalize according to spherical volume unit
# r, radial distance and p, polar angle
spherical_normalization <- function(r,p){ 1/( r^2 * sin(p) ) * 1/sum( r^2 * sin(p) )  }

no_normalization <- function(x) rep(1,length(x))




primary_modes_radial_3d_norm_diff <- function(a, b) {
	primary_modes_diff(a, b, weight_fun=radial_3d_normalization)
}

primary_modes_conical_3d_norm_diff <- function(a, b) {
	primary_modes_diff(a, b, weight_fun=radial_3d_normalization)
}




# These the coordinates for the Lambert-Azmuthal plots:
#   longitude (around) 30, 60, 90, 120 degrees
#   latitude  (in-out) 30, 60, 90, 120 degrees (starting from the positive x-axis)

# Equal Area Coordinate Grids
major_long_coords <- function(
	lat=c(pi/6, pi/2)
){
	transform(
		expand.grid(long=seq(0, 2*pi, length.out=200), lat=lat),
		capx = 2*sin(lat/2)*cos(long),
		capy = 2*sin(lat/2)*sin(long))
}

minor_long_coords <- function(
	lat=c(pi/3, pi*2/3)
) {
	transform(
		expand.grid(long=seq(0, 2*pi, length.out=200), lat=lat),
		capx = 2*sin(lat/2)*cos(long),
		capy = 2*sin(lat/2)*sin(long))
}

major_lat_coords <- function(
	n_long=5,
	lat=c(pi/3, pi*2/3)
) {
	transform(
		expand.grid(
			long=seq(0,2*pi, length.out=n_long), lat=c(pi/3, pi*2/3, length.out=200)),
		capx = 2*sin(lat/2)*cos(long),
		capy = 2*sin(lat/2)*sin(long))
}

minor_lat_coords <- function() {
	transform(
		expand.grid(long=seq(0,2*pi, length.out=5) + pi/4, lat=c(pi/3, pi*2/3, length.out=200)),
		capx = 2*sin(lat/2)*cos(long),
		capy = 2*sin(lat/2)*sin(long))
}

long_labels <- function() {
	transform(
		expand.grid(long=pi*3/2, lat=c(pi/6, pi/3, pi/2, pi*2/3)),
		label = as.character(round(180-180/pi*lat, 0)),
		capx = 2*sin(lat/2)*cos(long),
		capy = 2*sin(lat/2)*sin(long))
}

polar_equal_area_grids_bw <- function(
	scale=1,
	label_scale=1,
	box_bgcolor="#00007F",
	label_color="white",
	line_color="gray98",
	...
) {
	list(
		geom_path(
			data=minor_long_coords(),
			aes(x=capx, y=capy, group=lat),
			size=scale * .5, colour=line_color, ...),
		geom_path(
			data=minor_lat_coords(),
			aes(x=capx, y=capy, group=long),
			size=scale * .5, colour=line_color, ...),
		geom_path(
			data=major_long_coords(),
			aes(x=capx, y=capy, group=lat),
			size=scale * .2, colour=line_color, ...),
		geom_path(
			data=major_lat_coords(),
			aes(x=capx, y=capy, group=long),
			size=scale * .2, colour=line_color, ...),
		geom_label(
			data=long_labels(),
			aes(x=capx, y=capy, label=paste(label, sep="")),
			parse=T,
			size=label_scale*3, colour=label_color, fill=box_bgcolor),
    theme(
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank()))
}
