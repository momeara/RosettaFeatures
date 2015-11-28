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
id = "hydroxyl_proton_chi",
author = "Matthew O'Meara",
brief_description = "Measure the angle of the dihedral angle of the hydroxyl hydrogen and the best hydrogen bonded hydrogen relative to the plane of the aromatic ring for tyrosine amino acid types or BB atom for serine and theonine",

long_description = "
The dihedral angle in tyrosine residues defined by the atoms (CE2, C, OH, HH) may be rotameric because the C-OH bond is has partial sp2 character because of the resonance involving the aromatic ring. It is challening to definitively place the HH atom. Therefore use the dihedral angle (CE2, C, OH, H) where the final H is the hydrogen atom of the best hydrogen bond donor making a hydrogen bond with the OH group.

Currently the atoms extracted for hydrogen bond sites are the atm, base, bbase, base2. For the OH acceptors on tyrosines have acceptor chemical type hbacc_AHX (for aromatic hydroxyl) and have the following atoms:

  atm:  OH
  base: CZ
  bbase: CE2
  base2: HH

  H         HH
   `       /
     `    /
       `OH
        |
        |
        CZ
       /  \
      /    \
     |      CE2
     |  ()  |
     |      |
      \    /
       \ _/

  proton_chi = DIHEDRAL(CE2, CZ, OH, HH)
  hbond_chi = DIHEDRAL(CE2, CZ, OH, H)
",


feature_reporter_dependencies = c("ResidueFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <- "
SELECT
	don_site.HBChemType AS don_chem_type,
	don_site.HBChemType AS acc_chem_type,
	don_conf.chi3 AS chi3,
	don_conf.chi2 AS chi2,
	dssp_code.label AS don_dssp,
	CASE don_site.resNum - acc_site.resNum
		WHEN -1 THEN '-1' WHEN -2 THEN '-2' WHEN -3 THEN '-3' WHEN -4 THEN '-4'
		WHEN 1 THEN '1' WHEN 2 THEN '2' WHEN 3 THEN '3' WHEN 4 THEN '4'
		ELSE 'long' END AS seq_sep
FROM
	hbonds AS hb cross join
	hbond_sites AS don_site, hbond_sites AS acc_site cross join
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb cross join
	protein_residue_conformation AS don_conf cross join
	residue_secondary_structure AS don_dssp cross join
	dssp_codes AS dssp_code
WHERE
	hb.struct_id = acc_site.struct_id AND hb.acc_id = acc_site.site_id AND
	hb.struct_id = don_site.struct_id AND hb.don_id = don_site.site_id AND
	(don_site.HBChemType = 'hbdon_AHX' OR don_site.HBChemType = 'hbdon_HXL') AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	don_conf.struct_id = don_site.struct_id AND
	don_conf.seqpos = don_site.resNum AND
	don_dssp.struct_id = hb.struct_id AND don_dssp.resNum = don_site.resNum AND
	dssp_code.code = don_dssp.dssp;"
f <- query_sample_sources(sample_sources, sele)

f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type)
f <- na.omit(f, method="r")

dens <- estimate_density_1d_wrap(
  f[f$don_chem_type == "hbdon_AHX",], c("sample_source"), "chi3", xlim=c(-180, 180), adjust=.3)

plot_id = "tyr_proton_chi"
ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=x, y=y, colour=sample_source)) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  ggtitle("Proton Chi Torsional angle on Tyrosine Residues") +
  scale_x_continuous('Protein Chi (degrees)', breaks=c(-180, -90, 0, 90, 180)) +
  scale_y_continuous('Feature Density')
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d_wrap(
  f[f$don_chem_type == "hbdon_HXL",], c("sample_source"), "chi2", xlim=c(-180, 180), adjust=.3)

plot_id = "ser_thr_proton_chi"
ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=x, y=y, colour=sample_source)) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  ggtitle("Proton Chi Torsional angle on Tyrosine Residues") +
  scale_x_continuous('Protein Chi (degrees)', breaks=c(-180, -90, 0, 90, 180)) +
  scale_y_continuous('Feature Density')
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

###################

dens <- estimate_density_1d_wrap(
  f[f$don_chem_type == "hbdon_AHX",], c("sample_source", "don_dssp"), "chi3", xlim=c(-180, 180), adjust=.6)

plot_id = "tyr_proton_chi_by_dssp"
ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=x, y=y, colour=sample_source)) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source), ypos=.17) +
  ggtitle("Proton Chi Torsional angle on Tyrosine Residues") +
  scale_x_continuous('Protein Chi (degrees)', breaks=c(-180, -90, 0, 90, 180)) +
  scale_y_continuous('Feature Density') +
	facet_wrap(~don_dssp) +
	theme(legend.position=c(.8,.2))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d_wrap(
  f[f$don_chem_type == "hbdon_HXL",], c("sample_source", "don_dssp"), "chi2", xlim=c(-180, 180), adjust=.6)

plot_id = "ser_thr_proton_chi_by_dssp"
ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=x, y=y, colour=sample_source)) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  ggtitle("Proton Chi Torsional Angle on Serine/Theonine Residues") +
  scale_x_continuous('Protein Chi (degrees)', breaks=c(-180, -90, 0, 90, 180)) +
  scale_y_continuous('Feature Density') +
	facet_wrap(~don_dssp) +
	theme(legend.position=c(.8,.2))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


###################

dens <- estimate_density_1d_wrap(
  f[f$don_chem_type == "hbdon_AHX",], c("sample_source", "seq_sep"), "chi3", xlim=c(-180, 180), adjust=.6)

plot_id = "tyr_proton_chi_by_seq_sep"
ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=x, y=y, colour=sample_source)) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source), ypos=.17) +
  ggtitle("Proton Chi Torsional angle on Tyrosine Residues by Sequence Separation") +
  scale_x_continuous('Protein Chi (degrees)', breaks=c(-180, -90, 0, 90, 180)) +
  scale_y_continuous('Feature Density') +
	facet_wrap(~don_dssp) +
	theme(legend.position=c(.8,.2))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d_wrap(
  f[f$don_chem_type == "hbdon_HXL",], c("sample_source", "seq_sep"), "chi2", xlim=c(-180, 180), adjust=.6)

plot_id = "ser_thr_proton_chi_by_seq_sep"
ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=x, y=y, colour=sample_source)) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  ggtitle("Proton Chi Torsional Angle on Serine/Theonine Residues by Sequence Separation") +
  scale_x_continuous('Protein Chi (degrees)', breaks=c(-180, -90, 0, 90, 180)) +
  scale_y_continuous('Feature Density') +
	facet_wrap(~don_dssp) +
	theme(legend.position=c(.8,.2))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

###################

dens <- estimate_density_1d_wrap(
  f[f$don_chem_type == "hbdon_AHX",], c("sample_source", "acc_chem_type"), "chi3", xlim=c(-180, 180), adjust=.6)

plot_id = "tyr_proton_chi_by_acc_chem_type"
ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=x, y=y, colour=sample_source)) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source), ypos=.17) +
  ggtitle("Proton Chi Torsional angle on Tyrosine Residues by Acceptor Chemical Type") +
  scale_x_continuous('Protein Chi (degrees)', breaks=c(-180, -90, 0, 90, 180)) +
  scale_y_continuous('Feature Density') +
	facet_wrap(~don_dssp) +
	theme(legend.position=c(.8,.2))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d_wrap(
  f[f$don_chem_type == "hbdon_HXL",], c("sample_source", "acc_chem_typ"), "chi2", xlim=c(-180, 180), adjust=.6)

plot_id = "ser_thr_proton_chi_by_acc_chem_type"
ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=x, y=y, colour=sample_source)) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  ggtitle("Proton Chi Torsional Angle on Serine/Theonine Residues by Acceptor Chemical Type") +
  scale_x_continuous('Protein Chi (degrees)', breaks=c(-180, -90, 0, 90, 180)) +
  scale_y_continuous('Feature Density') +
	facet_wrap(~don_dssp) +
	theme(legend.position=c(.8,.2))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



})) # end FeaturesAnalysis
