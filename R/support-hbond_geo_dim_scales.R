# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


scale_x_AHdist <- scale_x_continuous(
	expression(paste('Acceptor -- Hydrogen Distance (', ring(A), ')')),
	limits=c(1.4, 3), breaks=c(1.4, 1.8, 2.2, 2.6, 3))

scale_y_AHdist <- scale_y_continuous(
	expression(paste('Acceptor -- Hydrogen Distance (', ring(A), ')')),
	limits=c(1.4, 3), breaks=c(1.4, 1.8, 2.2, 2.6, 3))

scale_x_AHdist3 <- scale_x_continuous(
	expression(paste('Acceptor -- Hydrogen Distance (', ring(A), ')')),
	limits=c(1.4^3, 3.2^3), breaks=c(1.5^3, 2^3, 3^3), labels=c(1.5, 2, 3))

scale_y_AHdist3 <- scale_y_continuous(
	expression(paste('Acceptor -- Hydrogen Distance (', ring(A), ')')),
	limits=c(1.4^3, 3.2^3), breaks=c(1.5^3, 2^3, 2.5^3, 3^3), labels=c(1.5, 2, 2.5, 3))

scale_x_ADdist <- scale_x_continuous(
	expression(paste('Acceptor -- Donor Distance (', ring(A), ')')),
	limits=c(2.4, 3.3), breaks=c(2.6, 2.9, 3.2))

scale_y_ADdist <- scale_y_continuous(
	expression(paste('Acceptor -- Donor Distance (', ring(A), ')')),
	limits=c(2.4, 3.3), breaks=c(2.3, 2.8, 3.3))

scale_x_ADdist3 <- scale_x_continuous(
	expression(paste('(Acceptor -- Donor Distance)^3')),
	limits=c(2.4^3, 3.3^3), breaks=c(2.6^3, 2.9^3, 3.2^3))

scale_y_ADdist3 <- scale_y_continuous(
	expression(paste('(Acceptor -- Donor Distance)^3')),
	limits=c(2.4^3, 3.3^3), breaks=c(2.3^3, 2.8^3, 3.3^3))

scale_x_cosBAH <- scale_x_continuous(
	"cos(Base -- Acceptor -- Hydrogen)",
	limit=c(-.3,1), breaks=c(-.2, .2, .6, 1))

scale_y_cosBAH <- scale_y_continuous(
	"Base -- Acceptor -- Hydrogen",
	limit=c(-.3,1), breaks=c(-.2, .2, .6, 1))

scale_x_cosAHD <- scale_x_continuous(
	"cos(Acceptor -- Hydrogen -- Donor)",
	limit=c(0,1), breaks=c(.2, .4, .6, .8, 1))

scale_y_cosAHD <- scale_y_continuous(
	"cos(Acceptor -- Hydrogen -- Donor)",
	limit=c(0,1), breaks=c(.2, .4, .6, .8, 1))

scale_x_chi <- scale_x_continuous(
	"Base -- Acceptor Torsion (Radians)",
	limit=c(-pi,pi), breaks=c(-pi, -pi*2/3, -pi/3, 0, pi/3, pi*2/3, pi))

scale_y_chi <- scale_y_continuous(
	"Base -- Acceptor Torsion (Radians)",
	limit=c(-pi,pi), breaks=c(-pi, -pi*2/3, -pi/3, 0, pi/3, pi*2/3, pi))

scale_x_chi_degrees <- scale_x_continuous(
	"Base -- Acceptor Torsion (Degrees)",
	limit=c(-180,180), breaks=c(-180, -90, 0, 90, 180))

scale_y_chi_degrees <- scale_y_continuous(
	"Base -- Acceptor Torsion (Degrees)",
	limit=c(-180,180), breaks=c(-180, -90, 0, 90, 180))


scale_x_energy <- scale_x_continuous(
	"Unweighted HBond Energy",
	limits=c(-1.5, 0), breaks = c(-1.5, -1, -.5, 0))

scale_y_energy <- scale_y_continuous(
	"Unweighted HBond Energy",
	limits=c(-1.5, 0), breaks = c(-1.5, -1, -.5, 0))

hbond_scale_x_by_name <- function(dim){
	if(as.character(dim)=="AHdist") scale_x_AHdist
	else if(as.character(dim)=="AHdist") scale_x_AHdist
	else if(as.character(dim)=="cosBAH") scale_x_cosBAH
	else if(as.character(dim)=="cosAHD") scale_x_cosAHD
	else if(as.character(dim)=="chi") scale_x_chi_degrees
}

hbond_scale_y_by_name <- function(dim){
	if(as.character(dim)=="AHdist") scale_y_AHdist
	else if(as.character(dim)=="AHdist") scale_y_AHdist
	else if(as.character(dim)=="cosBAH") scale_y_cosBAH
	else if(as.character(dim)=="cosAHD") scale_y_cosAHD
	else if(as.character(dim)=="chi") scale_y_chi_degrees
}


don_chem_type_name_linear <- function(don_chem_type){
	factor(don_chem_type,
		levels = c(
			"hbdon_GDE", "hbdon_GDH", "hbdon_AMO",
			"hbdon_IMD", "hbdon_IME", "hbdon_IND",
			"hbdon_HXL", "hbdon_AHX", "hbdon_PBA", "hbdon_CXA"),
		labels = c(
			"dGDE: R", "dGDH: R", "dAMO: K",
			"dIMD: H", "dIME: H", "dIND: W",
			"dHXL: S,T", "dAHX: Y", "dPBA: bb", "dCXA: N,Q"))
}

don_chem_type_name_wrap <- function(don_chem_type){
	factor(don_chem_type,
		levels = c(
			"hbdon_GDE", "hbdon_IMD", "hbdon_HXL", "hbdon_PBA",
			"hbdon_GDH", "hbdon_IME", "hbdon_AHX",  "hbdon_CXA",
			"hbdon_AMO", "hbdon_IND"),
		labels = c(
			"dGDE: R", "dIMD: H", "dHXL: S,T", "dPBA: bb",
			"dGDH: R", "dIME: H", "dAHX: Y", "dCXA: N,Q",
			"dAMO: K", "dIND: W"))
}

acc_chem_type_name_linear <- function(acc_chem_type){
	factor(acc_chem_type,
		levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
			"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
		labels = c("aIMD: H", "aIME: H", "aAHX: Y", "aHXL: S,T",
			"aCXA: N,Q", "aCXL: D,E", "aPBA: bb"))
}

acc_chem_type_name_wrap <- function(acc_chem_type){
	factor(acc_chem_type,
		levels = c("hbacc_HXL", "hbacc_CXL", "hbacc_IMD", "hbacc_AHX",
			"hbacc_CXA", "hbacc_IME",  "hbacc_PBA"),
		labels = c("aHXL: S,T", "aCXL: D,E", "aIMD: H",   "aAHX: Y",
			"aCXA: N,Q", "aIME: H",    "aPBA: bb"))
}

chem_type_name_linear <- function(chem_type){
	factor(chem_type,
		levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
			"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA",
			"hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
			"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
		labels = c("dIMD: H", "dIME: H", "dGDE: R", "dGDH: R",
			"dAHX: Y", "dHXL: S,T", "dIND: W", "dAMO: K", "dCXA: N,Q", "dPBA: bb",
			"aIMD: H", "aIME: H", "aAHX: Y", "aHXL: S,T",
			"aCXA: N,Q", "aCXL: D,E", "aPBA: bb"))
}

