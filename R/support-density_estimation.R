# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


#' @export
estimate_density_1d <-function(
  data,
  ids,
  variable,
  weight_fun=uniform_normalization,
  min_count=20,
  n_pts=512,
  histogram=FALSE,
	sample_domain=NULL,
	bw = "nrd0",
	adjust=1,
  ...){
	density.args <- list(...)
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

	if(is.null(sample_domain)){
		sample_domain <- range(data[,variable])
	}

  compute_density <- function(factor_df){
    if (nrow(factor_df) < min_count){
			#mjo I was having issues with whole columns/rows
			#disappearing. This worked for ggplot2 version 0.8.9 but now
			#gives warnings: asking "Do you need to adjust the group
			#aesthetic?" with version 0.9.0. I'm trying just returning an
			#empty data.frame to see if that works alright.
			#return( data.frame(x=seq(sample_domain[1], sample_domain[2], n_pts), y=0))

			return(data.frame())

    } else {
      weights <- weight_fun(factor_df[,variable])
      if(histogram){
        breaks = seq(from=sample_domain[1], to=sample_domain[2], length=n_pts)
        d <- plotrix::weighted.hist(x=factor_df[,variable], w=weights, breaks=breaks, plot=FALSE)
        return(data.frame(x=d$mids, y=d$density, counts=nrow(factor_df)))
      } else {
				#TODO general_kernel_adjust should be able to be set from the configuration file
				#adjust <- adjust * general_kernel_adjust
				tryCatch({
	        d <- do.call(stats::density,
						c(list(x=factor_df[,variable], from=sample_domain[1], to=sample_domain[2], n=n_pts,
	          weights=weights, bw=bw, adjust=adjust), density.args))
	          return(data.frame(x=d$x, y=d$y, counts=nrow(factor_df)))
				}, error=function(e){
					cat(paste("ERROR computing density for ids=(", paste(plyr::laply(factor_df[1,ids], as.character), collapse=", "), "): ", e, "\n", sep=""))
					return(data.frame(x=c(), y=c(), counts=c()))
				})
      }

    }
  }
  plyr::ddply(data, ids, compute_density)
}

#' @export
estimate_density_1d_wrap <-function(
  data,
  ids,
  variable,
  weight_fun=uniform_normalization,
  min_count=20,
  n_pts=512,
	xlim=c(0,360),
	bw="nrd0",
	adjust=1,
  ...){
	density.args <- list(...)
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

	extended_data <- plyr::mdply(
		c(xlim[1]-xlim[2], 0, xlim[2]-xlim[1]),function(s){
		y <- data;
		y[,variable] <- y[,variable] + s;
		y
	})

	extended_n_pts <- n_pts*3
	extended_min_count <- min_count*3

	# there is three times the data and assume the the bin width scales linearly with the sample size
	adjust <- adjust/3
  compute_density <- function(factor_df){
    if (nrow(factor_df) < extended_min_count){
      return( data.frame(x=seq(xlim[1], xlim[2], length.out=n_pts), y=0))
    } else {
      weights <- weight_fun(factor_df[,variable])
			d <- do.call(stats::density,
				c(list(x=factor_df[,variable], from=xlim[1], to=xlim[2], n=extended_n_pts,
					weights=weights, bw=bw, adjust=adjust), density.args))
      return(data.frame(
				x=d$x[xlim[1] <= d$x & d$x <= xlim[2]],
				y=d$y[xlim[1] <= d$x & d$x <= xlim[2]]*3,
				counts=nrow(factor_df)/3))
    }
  }
  plyr::ddply(extended_data, ids, compute_density)
}

#' @export
estimate_density_1d_reflect_boundary <-function(
  data,
  ids,
  variable,
  weight_fun=uniform_normalization,
  min_count=20,
  n_pts=512,
	reflect_left=FALSE,
	reflect_right=FALSE,
  left_boundary=NULL,
	right_boundary=NULL,
	bw="nrd0",
	adjust=1,
	verbose=FALSE,
  ...){
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

	if(is.null(left_boundary)){
		left_boundary=min(data[,variable])
	}
	if(is.null(right_boundary)){
		right_boundary=max(data[,variable])
	}
	extended_factor <- 1

	# Do the density estimation on the covering variable but apply the
	# weight funtion to the original variable
	data$covering_variable <- data[,variable]

	if(reflect_left==TRUE){
		data_lower <- data
		data_lower$covering_variable <- 2*left_boundary - data[,variable]
		extended_factor = extended_factor + 1
	}
	if(reflect_right==TRUE){
		data_upper <- data
		data_upper$covering_variable <- 2*right_boundary - data[,variable]
		extended_factor = extended_factor + 1
	}

	if(reflect_left==TRUE){
		data <- rbind(data, data_lower)
	}

	if(reflect_right==TRUE){
		data <- rbind(data, data_upper)
	}

	#TODO general_kernel_adjust should be able to be set from the configuration file
	#adjust <- adjust * general_kernel_adjust
	adjust <- adjust/(1 + reflect_left + reflect_right)

	if(verbose){
		if(reflect_left) {
			cat("Reflect LEFT: left_boundary = ", left_boundary, "\n")
		}
		if(reflect_right) {
			cat("Reflect RIGHT: right_boundary = ", right_boundary, "\n")
		}

		cat("Extended_factor = ", extended_factor, "\n")

		cat("Adjust = ", adjust, "\n")
		cat("data:")
		print(plyr::ddply(data, ids, function(df) data.frame(x=nrow(df))))
	}

  compute_density <- function(factor_df){
			if (nrow(factor_df) < min_count*extended_factor){
				return(data.frame(x=seq(left_boundary, right_boundary, length.out=n_pts), y=0, counts=nrow(factor_df)))
    } else {
      weights <- weight_fun(factor_df[,variable])
			d <- stats::density(
				x=factor_df$covering_variable,
				adjust=adjust,
				weights=weights,
				n=n_pts*extended_factor,
				from=left_boundary,
				to=right_boundary,
        ...)

			if(verbose){
				cat("nrow(factor_df) = ", nrow(factor_df), "\n")
				cat("nrow(factor_df)/extended_factor = ", nrow(factor_df)/extended_factor, "\n")
			}

			return(data.frame(
				x=d$x[left_boundary <= d$x & d$x <= right_boundary],
				y=d$y[left_boundary <= d$x & d$x <= right_boundary]*extended_factor,
				counts=round(nrow(factor_df)/extended_factor),0))
    }
  }
	z <- plyr::ddply(data, ids, compute_density)
	z[,!(names(z) %in% "X0")]
}

#' @export
estimate_density_1d_logspline <-function(
  data,
  ids,
  variable,
  min_count=20,
  n_pts=512,
  weight_fun=NULL,
  ...){

	if(!(class(data) == "data.frame")){
		stop(paste("The data argument must be a data.frame, instead it is of class '", class(data), "'", sep=""))
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

  xlim <- range(data[,variable])
  if(!is.null(weight_fun)){
    xlim_transformed <- weight_fun(xlim)
  }

  compute_density <- function(factor_df){
    if (nrow(factor_df) < min_count){
      return( data.frame(x=seq(xlim[1], xlim[2], n_pts), y=0))
    } else {
      if(!is.null(weight_fun)){
        values_transformed <- weight_fun(factor_df[,variable])
        lgs <- logspline::logspline(
          values_transformed,
          lbound=xlim_transformed[1],
          ubound=xlim_transformed[2])
        x_transformed <- seq(
          from=xlim_transformed[1],
          to=xlim_transformed[2],
          length.out=n_pts)
        y <- logspline::dlogspline(x_transformed, lgs)
        x <- seq(from=xlim[1], to=xlim[2], length.out=n_pts)

      } else {
        lgs <- logspline::logspline(factor_df[,variable], lbound=xlim[1], ubound=xlim[2])
        x <- seq(from=xlim[1], to=xlim[2], length.out=n_pts)
        y <- logspline::dlogspline(x, lgs)
      }
      return(data.frame(x=x, y=y, counts=nrow(factor_df)))
    }
  }
  plyr::ddply(data, ids, compute_density)
}

#' @export
estimate_density_2d <-function(
	data,
	ids,
	xvariable,
	yvariable,
	min_count=20,
	n_pts=512,
	bandwidth=NA,
	method='MASS::kde2',
	verbose=FALSE,
	scaled=FALSE,
	...){
	if(!(class(data) == "data.frame")){
		stop(paste("The data argument must be a data.frame, instead it is of class '", class(data), "'", sep=""))
	}
	if(nrow(data) == 0){
		stop(paste("Unable to compute density estimation because the data argument has no rows."))
	}
	for(id in ids){
		if(!(id %in% names(data))){
			stop(paste("The id variable '", id, "' is not a column name of the data. The ids are used to group the data instances for computing the density estimation.", sep=""))
		}
	}
	if(!method %in% c("ks::kde", "MASS::kde2", "histogram")){
		stop(paste("The method variable '", method, "' must be either 'ks::kde', 'MASS::kde2', or 'histogram'.", sep=""))
	}
	if(!(xvariable %in% names(data))){
		stop(paste("The value variable '", xvariable, "' is not a column name of the data. The xvariable and yvariable are used to compute the density estimation.", sep=""))
	}
	if(!(yvariable %in% names(data))){
		stop(paste("The value variable '", yvariable, "' is not a column name of the data. The xvariable and yvariable are used to compute the density estimation.", sep=""))
	}
	dens <- plyr::ddply(
		.data=data,
		.variables=ids,
		.fun=function(df){
		if(verbose){
			print(paste(
					"Estimating density 2d: ",
					paste(plyr::llply(df[1,ids], as.character), collapse=", "), " for '", nrow(df), "' instances.", sep=""))
		}
	  if (nrow(df) < min_count){
      return(data.frame(x=NULL, y=NULL, z=NULL, counts=NULL))
    } else {
			if(method=="histogram"){
			  h <- gplots::hist2d(
          x=as.matrix(df[,c(xvariable, yvariable)]), nbins=n_pts, show=FALSE)
        d <- with(h, data.frame(expand.grid(x=x, y=y), z=as.vector(counts), density=as.vector(counts)/nrow(df), count=nrow(df)))
        d$z <- melt(h$counts)$value
				return(d)
			} else if( method=="MASS::kde2"){
				if(is.na(bandwidth)){
					h <- c(MASS::bandwidth.nrd(df[,xvariable]), MASS::bandwidth.nrd(df[,yvariable]))
				} else {
					h <- bandwidth
				}
				dm <- MASS::kde2d(x=df[,xvariable], y=df[,yvariable], n=n_pts, h=h)
				if(verbose){
					print(summary(dm))
				}
				dz <- with(dm, data.frame(expand.grid(x=x, y=y), z=as.vector(z), counts=nrow(df)))
				if(verbose) {
					print(summary(dz))
				}
				return(dz)
      } else if( method=="ks::kde"){
				dm <- ks::kde(df[,c(xvariable, yvariable)])
				return(with(dm,
					data.frame(expand.grid(x=eval.points[1][[1]], y=eval.points[2][[1]]), z=estimate[[1]], counts=nrow(df))))
			}
    }
    d
  })
	if(scaled){
		dens <- plyr::ddply(
			.data=dens,
			.variables=ids,
			.drop=FALSE,
			transform, z=as.vector(z) / max(z))
	}
	return(dens)
}

#'  Compute quantiles for specified 'variable' grouping by 'ids' columns
#'  in 'variable'. The quantile for a real numer x is the fraction of
#'  values for 'variable' that are less than x.
#'
#'  If lengths(probs)==1, then "Ordinates for Probability Plotting" with
#'  'probs' probability points are used.
#'
#'  Returns data.frame with the following columns
#'    ids, probs, quantiles, counts
#'
#'  # for example:
#'  q <- compute_quantiles(f, c("id_col1", "id_col2"), "var")
#'
#'  ggplot(q) + geom_line(aes(x=quantiles*180/pi, y=probs))
#' @export
compute_quantiles <- function(
	data,
	ids,
	variable,
	probs=ppoints(1000),
	min_count=15
) {
	if(length(probs) == 1){
		probs <- ppoints(probs)
	}
	plyr::ddply(
		.data=data,
		.variables=ids,
		.fun=function(df){
		if(nrow(df) < min_count){
			return(data.frame())
		} else {
			return(data.frame(
				probs=probs,
				quantiles=quantile(df[,variable], probs=probs),
				counts=nrow(df)))
		}
	})
}

#' @export
compute_qq <- function(
	ref_data,
	new_data,
	group_ids,
	merge_ids,
	variable,
	probs=ppoints(1000)
) {
	if(length(probs) == 1){
		probs <- ppoints(probs)
	}

	ref_q <- plyr::ddply(
		.data=ref_data,
		.variables=group_ids,
		.run=function(df){
		data.frame(
			probs=probs, quantiles=quantile(df[,variable], probs=probs), counts=nrow(df))
	})

	new_q <- plyr::ddply(
		.data=new_data,
		.variables=group_ids,
		.fun=function(df){
		data.frame(
			probs=probs, quantiles=quantile(df[,variable], probs=probs), counts=nrow(df))
	})

	merge(ref_q, new_q, by=c(merge_ids, "probs"), suffixes=c(".ref", ".new"))
}

#'  Extract multiple, overlapping subsets of the data based on quantiles
#'
#'  Output: The rows of the input data duplicated once for each window
#'  it is contained in. The window id is in a column 'window.id' numbered
#'  1, 2, 3, ...
#'
#' [--------------------------------]
#' [----------]
#'       [----------]
#'            [----------]
#'                  [----------]
#'                       [----------]
#'
#' @export
sliding_windows <- function(
	data,
	id.vars,
	measure.var,
	n_windows=5,
	overlap_fraction=.3,
	verbose=FALSE
) {
	if(!(class(data) == "data.frame")){
		stop(paste("The data argument must be a data.frame, instead it is of class '", class(data), "'"))
	}

	if(!(measure.var %in% names(data))){
		stop(paste("The value measure.var '", measure.var, "' is not a column name of the data. The value measure.var is used to compute the sliding windows.", sep=""))
	}

	if(!is.numeric(n_windows)){
		stop(paste("n_windows must be numeric. Instead it is of class,", class(n_windows), ".", sep=""))
	}

	if(n_windows <= 0 || is.na(n_windows)){
		stop(paste("n_windows must be a number greater than zero. Instead it has value '", n_windows, "'", sep=""))
	}

if(!is.numeric(overlap_fraction)){
		stop(paste("overlap_fraction must be numeric. Instead it is of class,", class(overlap_fraction), ".", sep=""))
	}

	if(overlap_fraction <= 0 || overlap_fraction >= 1 || is.na(n_windows)){
		stop(paste("overlap_fraction must be a number in (0, 1). Instead it has value '", overlap_fraction, "'", sep=""))
	}

	# total_quantiles = total_width_of_windows - double_counted_regions
	# 1 = width * n_windows - (n_windows - 1) * (width * overlap_fraction)
	# width = 1 / ( n_windows - (n_windows - 1) * overlap_fraction)
	# xmin = { i*(width - overlap_fraction) | i in [0, n_windows - 1] }
	# xmax = { width + i*(width - overlap_fraction) | i in [0, n_windows - 1] }

	width <- 1L / ( n_windows - (n_windows - 1L) * overlap_fraction )
	windows <- data.frame(
		window.id = factor(seq(1L, n_windows)),
		pmin = seq(0,n_windows-1L) * (width - width*overlap_fraction),
		pmax = width + seq(0,n_windows-1L) * (width - width*overlap_fraction))

	windows$rmin <- as.numeric(quantile(data[,measure.var], probs=windows$pmin))
	windows$rmax <- as.numeric(quantile(data[,measure.var], probs=windows$pmax))


	if(verbose){
		cat(paste("id.vars: [", paste(id.vars, collapse=", ", sep=""), "]\n", sep=""))
		cat(paste("n_windows: ", n_windows, "\n", sep=""))
		cat(paste("overlap_fraction: ", overlap_fraction, "\n", sep=""))
		cat(paste("width: ", width, "\n", sep=""))
		print(windows)
	}

	plyr::ddply(data, id.vars, function(df) {
		plyr::adply(windows, 1L, function(window) {
			sub_df <- df[
				df[,measure.var] >= window$rmin[1] &
				df[,measure.var] <= window$rmax[1],]
			sub_df$windows <-
				factor(paste(round(window$rmin[1],2), "-", round(window$rmax[1],2), sep=""))
			if(verbose){
				cat(paste(
					"window.id: '", window$window.id[1], "'\n",
					"(pmin, pmax): '", window[1,"pmin"], "', '", window[1,"pmax"], "'\n",
					"(xmin, xmax): '", window[1,"rmin"], "', '", window[1,"rmax"], "'\n",
					"count: '", nrow(sub_df), "'\n", sep=""))
			}
			sub_df
		})
	})
}

#' @export
distance_matrix <- function(
	data,
	id.var,
	measure.var,
	cmp_fun,
	verbose=FALSE
) {
	plyr::daply(data, id.var, function(da) {
		mid_z <- plyr::daply(data, id.var, function(db) {
			z <- cmp_fun(da[,measure.var], db[,measure.var])$statistic
			if(verbose){
				print(paste("    inner: ", z, sep=""))
			}
			z
		})
		if(verbose){
			print(paste("  mid: ", mid_z, sep=""))
		}
		mid_z
	})
}
