# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#' @export
compute_summary_statistics_1d <- function(
	data,
	ids,
	variable,
	statistic_functions,
	density_estimation_function=estimate_density_1d,
	density_estimation_args=list(),
	...
) {
	extra.args <- list(...)
	data <- as.data.frame(data)
	if(!(class(data) == "data.frame")){
		stop(paste("The data argument must be a data.frame, instead it is of class '", class(data), "'"))
	}
	if(nrow(data) == 0){
		stop(paste("Unable to compute density estimation because the data argument has no rows."))
	}
	for(id in ids){
		if(!(id %in% names(data))){
			stop(paste("The id variable '", id, "' is not a column name of the data. The ids are used to group the data instances for computing the density estimation.", sep=""))
		}
	}
	if(!(variable %in% names(data))){
		stop(paste("The value variable '", variable, "' is not a column name of the data. The value variable is used to compute the density estimation.", sep=""))
	}

	plyr::ldply(statistic_functions, function(stat_fun){
		plyr::ddply(data, ids, function(df){
			z <- stat_fun(df, variable, extra.args)
			data.frame(statistic=names(z), value=z[1])
		})
	})
}
