#-*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

locate_rosetta_application <- function(
	app_name,
	rosetta_base_path = NULL,
	platform=NULL,
	extras="default",
	compiler="clang",
	mode="release") {

	# find application is located at
	# ${rosetta_base_dir}/rosetta_source/bin/${app_name}.${extras}.${platform}${compiler}${mode}
	# detect what is not specified

	if(is.null(rosetta_base_path)){
		rosetta_base_path <- file.path(base_dir, "..", "..")
	}

	if(is.null(platform)){
		sysname <- Sys.info()[1]
		if (sysname == "Linux") {
			platform = "linux"
		} else if (sysname == "Darwin") {
			platform = "mac"
		} else {
			stop(paste("Unable to determine platform in trying to locate the application '", app_name, ". Try specifying the platform explicitly. (e.g. platform='linux' or platform='mac')", sep=""))
		}
	}

	full_app_path = file.path(rosetta_base_path, "rosetta_source", "bin", paste(app_name, ".", extras, ".", platform, compiler, mode, sep=""))
	if(!file.exists(full_app_path)){
		stop(paste("Looking for application '", app_name, "' at the path '", full_app_path, "', but it does not exist.", sep=""))
	}
	full_app_path
}
