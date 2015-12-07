# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

library(ggplot2)


source("../hbonds/hbond_geo_dim_scales.R")

feature_analyses <- c(feature_analyses, methods::new("FeaturesAnalysis",
id = "charge_charge_energies",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "ResidueTypeFeatures", "ChargeChargeFeatures"),
run=function(self, sample_sources, output_dir, output_formats){




base_charge <- Vectorize(function(hb_chem_type){
	switch(as.character(hb_chem_type),
		hbdon_PBA=-0.35, # ALA   ATOM  N   Nbb  NH1  -0.47  -0.350
		hbdon_CXA=-0.78, # ASN   ATOM  ND2 NH2O NH2  -0.62  -0.780
		hbdon_IMD=-0.40, # HIS_D ATOM  ND1 Ntrp NR1  -0.36  -0.400
		hbdon_IME=-0.40, # HIS   ATOM  NE2 Ntrp NR1  -0.36  -0.400
		hbdon_IND=-0.40, # TRP   ATOM  NE1 Ntrp NY   -0.61  -0.400
		hbdon_AMO=-0.32, # LYS   ATOM  NZ  Nlys NH3  -0.3   -0.320
		hbdon_GDE=-0.35, # ARG   ATOM  NE  Narg NC2  -0.7   -0.350
		hbdon_GDH=-0.70, # ARG   ATOM  NH2 Narg NC2  -0.8   -0.700
		hbdon_AHX=-0.49, # TYR   ATOM  OH  OH   OH1  -0.54  -0.490
		hbdon_HXL=-0.49, # SER   ATOM  OG  OH   OH1  -0.66  -0.490
		hbacc_PBA= 0.55, # ALA   ATOM  C   CObb C     0.51   0.550
		hbacc_CXA= 0.55, # ASN   ATOM  CG  CNH2 CC    0.55   0.550
		hbacc_CXL= 0.10, # ASP   ATOM  CG  COO  CC    0.62   0.100
		hbacc_IMD= 0.155,# HIS   ATOM  CG  aroC CPH1  0.22   0.155
		hbacc_IME= 0.155,# HIS_D ATOM  CE1 aroC CPH2  0.25   0.155
		hbacc_AHX= 0.055,# TYR   ATOM  CZ  aroC CA    0.11   0.055
		hbacc_HXL= 0)    # SER   ATOM  CB  CH2  CT2   0.05   0.000 # THR has 0.055
})

atom_charge <- Vectorize(function(hb_chem_type){
	switch(as.character(hb_chem_type),
		hbdon_PBA= 0.25, # ALA   ATOM  H   HNbb H     0.31   0.250
		hbdon_CXA= 0.39, # ASN   ATOM 1HD2 Hpol H     0.32   0.390
		hbdon_IMD= 0.40, # HIS_D ATOM  HD1 Hpol H     0.32   0.400
		hbdon_IME= 0.40, # HIS   ATOM  HE2 Hpol H     0.32   0.400
		hbdon_IND= 0.40, # TRP   ATOM  HE1 Hpol H     0.38   0.400
		hbdon_AMO= 0.33, # LYS   ATOM 3HZ  Hpol HC    0.33   0.330
		hbdon_GDE= 0.45, # ARG   ATOM  HE  Hpol HC    0.44   0.450
		hbdon_GDH= 0.40, # ARG   ATOM 1HH1 Hpol HC    0.46   0.400
		hbdon_AHX= 0.43, # TYR   ATOM  HH  Hpol H     0.43   0.435
		hbdon_HXL= 0.49, # SER   ATOM  HG  Hpol H     0.43   0.490
		hbacc_PBA=-0.55, # ALA   ATOM  O   OCbb O    -0.51  -0.550
		hbacc_CXA=-0.55, # ASN   ATOM  OD1 ONH2 O    -0.55  -0.550
		hbacc_CXL=-0.55, # ASP   ATOM  OD1 OOC  OC   -0.76  -0.550
		hbacc_IMD=-0.56, # HIS   ATOM  ND1 Nhis NR2  -0.7   -0.560
		hbacc_IME=-0.56, # HIS_D ATOM  NE2 Nhis NR2  -0.7   -0.560
		hbacc_AHX=-0.49, # TYR   ATOM  OH  OH   OH1  -0.54  -0.490
		hbacc_HXL=-0.49) # SER   ATOM  OG  OH   OH1  -0.66  -0.490
})

coulomb_energy <- function(d, q1, q2, max_dist=5.5){
	die <- 10
	C0 <- 322.0637
	C1 <- C0 / die
	C2 <- C1 / (max_dist * max_dist)
	q1 * q2 * (C1 / (d*d) - C2)
}

sele <-"
CREATE TEMPORARY TABLE hb_chem_type_charges (
	HBChemType TEXT,
	base_charge REAL,
	atom_charge REAL,
	PRIMARY KEY (HBChemType));

INSERT INTO hb_chem_type_charges VALUES ('hbdon_PBA', -0.35 , 0.25);
INSERT INTO hb_chem_type_charges VALUES ('hbdon_CXA', -0.78 , 0.39);
INSERT INTO hb_chem_type_charges VALUES ('hbdon_IMD', -0.40 , 0.40);
INSERT INTO hb_chem_type_charges VALUES ('hbdon_IME', -0.40 , 0.40);
INSERT INTO hb_chem_type_charges VALUES ('hbdon_IND', -0.40 , 0.40);
INSERT INTO hb_chem_type_charges VALUES ('hbdon_AMO', -0.32 , 0.33);
INSERT INTO hb_chem_type_charges VALUES ('hbdon_GDE', -0.35 , 0.45);
INSERT INTO hb_chem_type_charges VALUES ('hbdon_GDH', -0.70 , 0.40);
INSERT INTO hb_chem_type_charges VALUES ('hbdon_AHX', -0.49 , 0.43);
INSERT INTO hb_chem_type_charges VALUES ('hbdon_HXL', -0.49 , 0.49);
INSERT INTO hb_chem_type_charges VALUES ('hbacc_PBA',  0.55 ,-0.55);
INSERT INTO hb_chem_type_charges VALUES ('hbacc_CXA',  0.55 ,-0.55);
INSERT INTO hb_chem_type_charges VALUES ('hbacc_CXL',  0.10 ,-0.55);
INSERT INTO hb_chem_type_charges VALUES ('hbacc_IMD',  0.155,-0.56);
INSERT INTO hb_chem_type_charges VALUES ('hbacc_IME',  0.155,-0.56);
INSERT INTO hb_chem_type_charges VALUES ('hbacc_AHX',  0.055,-0.49);
INSERT INTO hb_chem_type_charges VALUES ('hbacc_HXL',  0    ,-0.49);

CREATE INDEX IF NOT EXISTS hbonds_struct_id_don_id_acc_id ON hbonds (
	struct_id, don_id, acc_id);

SELECT
	cc_charge_coords.s1_chem_type, cc_charge_coords.s2_chem_type,
	cc_charge_coords.s1bq AS s1bq, cc_charge_coords.s1aq AS s1aq,
	cc_charge_coords.s2bq AS s2bq, cc_charge_coords.s2aq AS s2aq,
 	cc_charge_coords.s1bx AS s1bx, cc_charge_coords.s1by AS s1by, cc_charge_coords.s1bz AS s1bz,
	cc_charge_coords.s1ax AS s1ax, cc_charge_coords.s1ay AS s1ay, cc_charge_coords.s1az AS s1az,
 	cc_charge_coords.s2bx AS s2bx, cc_charge_coords.s2by AS s2by, cc_charge_coords.s2bz AS s2bz,
	cc_charge_coords.s2ax AS s2ax, cc_charge_coords.s2ay AS s2ay, cc_charge_coords.s2az AS s2az,
	hb.energy AS hb_energy
FROM
	(SELECT
		cc.struct_id, s1.site_id AS s1_id, s2.site_id AS s2_id,
		s1.HBChemType AS s1_chem_type,
		s2.HBChemType AS s2_chem_type,
		s1q.base_charge AS s1bq, s1q.atom_charge AS s1aq,
		s2q.base_charge AS s2bq, s2q.atom_charge AS s2aq,
	 	s1_atoms.base_x AS s1bx, s1_atoms.base_y AS s1by, s1_atoms.base_z AS s1bz,
		s1_atoms.atm_x  AS s1ax, s1_atoms.atm_y  AS s1ay, s1_atoms.atm_z  AS s1az,
	 	s2_atoms.base_x AS s2bx, s2_atoms.base_y AS s2by, s2_atoms.base_z AS s2bz,
		s2_atoms.atm_x  AS s2ax, s2_atoms.atm_y  AS s2ay, s2_atoms.atm_z  AS s2az
	FROM
		charge_charge_pairs AS cc,
		hbond_sites AS s1,
		hbond_sites AS s2,
		hb_chem_type_charges AS s1q,
		hb_chem_type_charges AS s2q,
--		hbond_sites_pdb AS s1_pdb,
--		hbond_sites_pdb AS s2_pdb,
		hbond_site_atoms AS s1_atoms,
		hbond_site_atoms AS s2_atoms
	WHERE
		s1.struct_id = cc.struct_id AND
		s1.site_id = cc.q1_site_id AND
		s2.struct_id = cc.struct_id AND
		s2.site_id = cc.q2_site_id AND
		s1q.HBChemType = s1.HBChemType AND
		s2q.HBChemType = s2.HBChemType AND
--		s1_pdb.struct_id = cc.struct_id AND
--		s1_pdb.site_id = cc.q1_site_id AND
--		s2_pdb.struct_id = cc.struct_id AND
--		s2_pdb.site_id = cc.q2_site_id AND
--		s1_pdb.heavy_atom_temperature < 30 AND
--		s2_pdb.heavy_atom_temperature < 30 AND
		s1_atoms.struct_id = cc.struct_id AND
		s1_atoms.site_id = cc.q1_site_id AND
		s2_atoms.struct_id = cc.struct_id AND
		s2_atoms.site_id = cc.q2_site_id) AS cc_charge_coords LEFT JOIN
		hbonds AS hb ON
			hb.struct_id = cc_charge_coords.struct_id AND
			((hb.don_id = cc_charge_coords.s1_id AND hb.acc_id = cc_charge_coords.s2_id) OR
			(hb.don_id = cc_charge_coords.s2_id AND hb.acc_id = cc_charge_coords.s1_id));"

f <-  query_sample_sources(sample_sources, sele)

f <- transform(f,
	s1b_s2b_dist=vector_distance(cbind(s1bx, s1by, s1bz), cbind(s2bx, s2by, s2bz)),
	s1b_s2a_dist=vector_distance(cbind(s1bx, s1by, s1bz), cbind(s2ax, s2ay, s2az)),
	s1a_s2b_dist=vector_distance(cbind(s1ax, s1ay, s1az), cbind(s2bx, s2by, s2bz)),
	s1a_s2a_dist=vector_distance(cbind(s1ax, s1ay, s1az), cbind(s2ax, s2ay, s2az)))

z <- with(f,
	data.frame(
		sample_source = sample_source,
		hb_energy = hb_energy,
		s1_chem_type = s1_chem_type,
		s2_chem_type = s2_chem_type,
		coulomb =
			coulomb_energy(s1b_s2b_dist, s1bq, s2bq) +
			coulomb_energy(s1b_s2a_dist, s1bq, s2aq) +
			coulomb_energy(s1a_s2b_dist, s1aq, s2bq) +
			coulomb_energy(s1a_s2a_dist, s1aq, s2aq),
		neg_coulomb =
			pmin(0, coulomb_energy(s1b_s2b_dist, s1bq, s2bq)) +
			pmin(0, coulomb_energy(s1b_s2a_dist, s1bq, s2aq)) +
			pmin(0, coulomb_energy(s1a_s2b_dist, s1aq, s2bq)) +
			pmin(0, coulomb_energy(s1a_s2a_dist, s1aq, s2aq)),
		pos_coulomb =
			pmax(0, coulomb_energy(s1b_s2b_dist, s1bq, s2bq)) +
			pmax(0, coulomb_energy(s1b_s2a_dist, s1bq, s2aq)) +
			pmax(0, coulomb_energy(s1a_s2b_dist, s1aq, s2bq)) +
			pmax(0, coulomb_energy(s1a_s2a_dist, s1aq, s2aq))))

dens <- estimate_density_1d(z[z$coulomb < .2 & z$coulomb > -.2 & (z$coulomb > .000000001 | z$coulomb < -.0000000001),],
	c("sample_source"), "coulomb", weight_fun=uniform_normalization)

plot_id <- "charge_charge_pair_energies"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Charge-Charge Energies; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)") +
	scale_y_continuous("Feature Density")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d(z[z$coulomb < .3 & z$coulomb > -3 & (z$coulomb > .000000001 | z$coulomb < -.0000000001),],
	c("sample_source"), "coulomb", weight_fun=uniform_normalization)

plot_id <- "charge_charge_pair_energies_all_hb_range"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=log(y+1), colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Charge-Charge Energies HB-Range; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)", limits=c(-3, .3)) +
	scale_y_continuous("Feature Density")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d(z[!is.na(z$hb_energy) & z$coulomb < .3 & z$coulomb > -3 & (z$coulomb > .000000001 | z$coulomb < -.0000000001),],
	c("sample_source"), "coulomb", weight_fun=uniform_normalization)
plot_id <- "charge_charge_pair_energies_hb"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Charge-Charge Energies; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)",limits=c(-3, .3)) +
	scale_y_continuous("Feature Density")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

dens <- estimate_density_1d(z[is.na(z$hb_energy) & z$coulomb < .3 & z$coulomb > -3 & (z$coulomb > .000000001 | z$coulomb < -.0000000001),],
	c("sample_source"), "coulomb", weight_fun=uniform_normalization)
plot_id <- "charge_charge_pair_energies_no_hb"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=log(y+1), colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Charge-Charge Energies; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)", limits=c(-3, .3)) +
	scale_y_continuous("Feature Density")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d(z[z$neg_coulomb > -.18 & z$neg_coulomb < -.0000001,],
	c("sample_source"), "neg_coulomb", weight_fun=uniform_normalization)

plot_id <- "charge_charge_pair_neg_energies"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Attractive Charge-Charge Energies; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)") +
	scale_y_continuous("Feature Density")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d(z[z$pos_coulomb < .07 & z$pos_coulomb > .00000001,],
	c("sample_source"), "pos_coulomb", weight_fun=uniform_normalization)

plot_id <- "charge_charge_pair_pos_energies"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Repulsive Charge-Charge Energies; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)") +
	scale_y_continuous("Feature Density")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



dens <- estimate_density_1d(z[!is.na(z$hb_energy) & z$neg_coulomb > -3 & z$neg_coulomb < -.0000001,],
	c("sample_source"), "neg_coulomb", weight_fun=uniform_normalization)

plot_id <- "charge_charge_pair_neg_energies_hbonded"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Attractive Charge-Charge Energies for H-Bonded Groups; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)") +
	scale_y_continuous("Feature Density")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d(z[!is.na(z$hb_energy) & z$pos_coulomb < 2 & z$pos_coulomb > .0000001,],
	c("sample_source"), "pos_coulomb", weight_fun=uniform_normalization)

plot_id <- "charge_charge_pair_pos_energies_hbonded"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Repulsive Charge-Charge Energies for H-Bonded Groups; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)") +
	scale_y_continuous("Feature Density")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

#################
dens <- estimate_density_1d(z[z$coulomb < .2 & z$coulomb > -.2 & (z$coulomb > .000000001 | z$coulomb < -.0000000001),],
	c("sample_source", "s1_chem_type"), "coulomb", weight_fun=uniform_normalization)

plot_id <- "charge_charge_pair_energies_s1_chem_type"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Charge-Charge Energies; B-Fact < 30") +
	facet_wrap(~s1_chem_type) +
	scale_x_continuous("Coulomb Between Polar Sites (REU)") +
	scale_y_continuous("Feature Density")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d(z[z$coulomb < .3 & z$coulomb > -3 & (z$coulomb > .000000001 | z$coulomb < -.0000000001),],
	c("sample_source", "s1_chem_type"), "coulomb", weight_fun=uniform_normalization)

plot_id <- "charge_charge_pair_energies_all_hb_range_s1_chem_type"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=log(y+1), colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Charge-Charge Energies HB-Range; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)", limits=c(-3, .3)) +
	scale_y_continuous("Feature Density") +
	facet_wrap(~s1_chem_type)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d(z[!is.na(z$hb_energy) & z$coulomb < .3 & z$coulomb > -3 & (z$coulomb > .000000001 | z$coulomb < -.0000000001),],
	c("sample_source", "s1_chem_type"), "coulomb", weight_fun=uniform_normalization)
plot_id <- "charge_charge_pair_energies_hb_s1_chem_type"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Charge-Charge Energies; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)",limits=c(-3, .3)) +
	scale_y_continuous("Feature Density") +
	facet(~s1_chem_type)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

dens <- estimate_density_1d(z[is.na(z$hb_energy) & z$coulomb < .3 & z$coulomb > -3 & (z$coulomb > .000000001 | z$coulomb < -.0000000001),],
	c("sample_source", "s1_chem_type"), "coulomb", weight_fun=uniform_normalization)
plot_id <- "charge_charge_pair_energies_no_hb_s1_chem_type"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=log(y+1), colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Charge-Charge Energies; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)", limits=c(-3, .3)) +
	scale_y_continuous("Feature Density") +
	facet_wrap(~s1_chem_type)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d(z[z$neg_coulomb > -.18 & z$neg_coulomb < -.0000001,],
	c("sample_source", "s1_chem_type"), "neg_coulomb", weight_fun=uniform_normalization)

plot_id <- "charge_charge_pair_neg_energies_s1_chem_type"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Attractive Charge-Charge Energies; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)") +
	scale_y_continuous("Feature Density") +
	facet_wrap(~s1_chem_type)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d(z[z$pos_coulomb < .07 & z$pos_coulomb > .00000001,],
	c("sample_source", "s1_chem_type"), "pos_coulomb", weight_fun=uniform_normalization)

plot_id <- "charge_charge_pair_pos_energies_s1_chem_type"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Repulsive Charge-Charge Energies; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)") +
	scale_y_continuous("Feature Density") +
	facet_wrap(~s1_chem_type)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



dens <- estimate_density_1d(z[!is.na(z$hb_energy) & z$neg_coulomb > -3 & z$neg_coulomb < -.0000001,],
	c("sample_source", "s1_chem_type"), "neg_coulomb", weight_fun=uniform_normalization)

plot_id <- "charge_charge_pair_neg_energies_hbonded_s1_chem_type"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Attractive Charge-Charge Energies for H-Bonded Groups; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)") +
	scale_y_continuous("Feature Density") +
	facet_wrap(~s1_chem_type)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d(z[!is.na(z$hb_energy) & z$pos_coulomb < 2 & z$pos_coulomb > .0000001,],
	c("sample_source", "s1_chem_type"), "pos_coulomb", weight_fun=uniform_normalization)

plot_id <- "charge_charge_pair_pos_energies_hbonded_s1_chem"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Repulsive Charge-Charge Energies for H-Bonded Groups; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)") +
	scale_y_continuous("Feature Density") +
	facet_wrap(~s1_chem_type)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


#################

dens <- estimate_density_1d(z[z$coulomb < .2 & z$coulomb > -.2 & (z$coulomb > .000000001 | z$coulomb < -.0000000001),],
	c("sample_source", "s1_chem_type", "s2_chem_type"), "coulomb", weight_fun=uniform_normalization)

plot_id <- "charge_charge_pair_energies_s1_chem_type_s2_chem_type"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Charge-Charge Energies; B-Fact < 30") +
	facet_grid(s1_chem_type ~ s2_chem_type) +
	scale_x_continuous("Coulomb Between Polar Sites (REU)") +
	scale_y_continuous("Feature Density")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d(z[z$coulomb < .3 & z$coulomb > -3 & (z$coulomb > .000000001 | z$coulomb < -.0000000001),],
	c("sample_source", "s1_chem_type", "s2_chem_type"), "coulomb", weight_fun=uniform_normalization)

plot_id <- "charge_charge_pair_energies_all_hb_range_s1_chem_type_s2_chem_type"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=log(y+1), colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Charge-Charge Energies HB-Range; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)", limits=c(-3, .3)) +
	scale_y_continuous("Feature Density") +
	facet_grid(s1_chem_type ~ s2_chem_type)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d(z[!is.na(z$hb_energy) & z$coulomb < .3 & z$coulomb > -3 & (z$coulomb > .000000001 | z$coulomb < -.0000000001),],
	c("sample_source", "s1_chem_type", "s2_chem_type"), "coulomb", weight_fun=uniform_normalization)
plot_id <- "charge_charge_pair_energies_hb_s1_chem_type_s2_chem_type"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Charge-Charge Energies; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)",limits=c(-3, .3)) +
	scale_y_continuous("Feature Density") +
	facet_grid(s1_chem_type ~ s2_chem_type)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

dens <- estimate_density_1d(z[is.na(z$hb_energy) & z$coulomb < .3 & z$coulomb > -3 & (z$coulomb > .000000001 | z$coulomb < -.0000000001),],
	c("sample_source", "s1_chem_type", "s2_chem_type"), "coulomb", weight_fun=uniform_normalization)
plot_id <- "charge_charge_pair_energies_no_hb_s1_chem_type_s2_chem_type"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=log(y+1), colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Charge-Charge Energies; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)", limits=c(-3, .3)) +
	scale_y_continuous("Feature Density") +
	facet_grid(s1_chem_type ~ s2_chem_type)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d(z[z$neg_coulomb > -.18 & z$neg_coulomb < -.0000001,],
	c("sample_source", "s1_chem_type", "s2_chem_type"), "neg_coulomb", weight_fun=uniform_normalization)

plot_id <- "charge_charge_pair_neg_energies_s1_chem_type_s2_chem_type"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Attractive Charge-Charge Energies; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)") +
	scale_y_continuous("Feature Density") +
	facet_grid(s1_chem_type ~ s2_chem_type)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d(z[z$pos_coulomb < .07 & z$pos_coulomb > .00000001,],
	c("sample_source", "s1_chem_type", "s2_chem_type"), "pos_coulomb", weight_fun=uniform_normalization)

plot_id <- "charge_charge_pair_pos_energies_s1_chem_type_s2_chem_type"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Repulsive Charge-Charge Energies; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)") +
	scale_y_continuous("Feature Density") +
	facet_grid(s1_chem_type ~ s2_chem_type)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



dens <- estimate_density_1d(z[!is.na(z$hb_energy) & z$neg_coulomb > -3 & z$neg_coulomb < -.0000001,],
	c("sample_source", "s1_chem_type", "s2_chem_type"), "neg_coulomb", weight_fun=uniform_normalization)

plot_id <- "charge_charge_pair_neg_energies_hbonded_s1_chem_type_s2_chem_type"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Attractive Charge-Charge Energies for H-Bonded Groups; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)") +
	scale_y_continuous("Feature Density") +
	facet_grid(s1_chem_type ~ s2_chem_type)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d(z[!is.na(z$hb_energy) & z$pos_coulomb < 2 & z$pos_coulomb > .0000001,],
	c("sample_source", "s1_chem_type", "s2_chem_type"), "pos_coulomb", weight_fun=uniform_normalization)

plot_id <- "charge_charge_pair_pos_energies_hbonded_s1_chem_type_s2_chem_type"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=counts, group=sample_source, colour=sample_source)) +
	ggtitle("Repulsive Charge-Charge Energies for H-Bonded Groups; B-Fact < 30") +
	scale_x_continuous("Coulomb Between Polar Sites (REU)") +
	scale_y_continuous("Feature Density") +
	facet_grid(s1_chem_type ~ s2_chem_type)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



})) # end FeaturesAnalysis
