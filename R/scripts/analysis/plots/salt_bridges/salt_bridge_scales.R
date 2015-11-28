# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.




scale_x_rho <- list(
	scale_x_continuous(expression(
		paste('Cental Carbon -- Oxygen Distance (', ring(A), ')')), limit=c(2,6)))

scale_y_rho <- list(
	scale_y_continuous(expression(
		paste('Cental Carbon -- Oxygen Distance (', ring(A), ')')), limit=c(2,6)))

scale_x_psi <- list(scale_x_continuous("Angle Around Donor (Degrees)"))
scale_y_psi <- list(scale_y_continuous("Angle Around Donor (Degrees)"))

scale_x_theta <- list(scale_x_continuous("Angle Out of Donor Plane (Degrees)"))
scale_y_theta <- list(scale_y_continuous("Angle Out of Donor Plane (Degrees)"))
