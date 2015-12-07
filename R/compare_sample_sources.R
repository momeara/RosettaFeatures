#!/usr/bin/env Rscript
# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



load_config_file <- function(config_filename, verbose=F){
	if(!file.exists(config_filename)){
		cat("ERROR: Config file '", config_filename, "' does not exist.\n", sep="")
		stop(1)
	}

	tryCatch({
		configuration <- jsonlite::fromJSON(config_filename)
		return(configuration)
	}, error=function(e){
		cat(
			"ERROR: Unable to parse configuration file '", as.character(config_filename), "'\n",
			"failed with the following error:\n",
			e, sep="")
	})
  return(NULL)
}

initialize_output_dir <- function(configuration, verbose=F){
	if(!("output_dir" %in% names(configuration)) || !is.character(configuration$output_dir)){
		configuration$output_dir <- "build"
		cat("WARNING: output_dir not specified. Using default of 'build'\n")
	} else {
		configuration$output_dir <- as.character(configuration$output_dir)
		cat("\toutput_dir: ", configuration$output_dir, "\n", sep="")
	}
	if(!file.exists(file.path(configuration$output_dir))){
		if(verbose){
			cat("output_dir does not exist creating it...\n")
		}
		dir.create(file.path(configuration$output_dir), recursive=TRUE)
	}

	configuration
}

initialize_db_engine <- function(db_engine=NULL, configuration, verbose=F){
	if(is.null(db_engine)){
		if(verbose){
			cat("Using SQLite as the database engine.\n")
		}
		db_engine <- RSQLite::SQLite()
	}
	db_engine
}


initialize_sample_sources <- function(
	configuration,
	db_engine,
	default_cache_size = 10000,
	verbose=F,
	initialize_databases=T
){

	sample_sources_error_msg <- "To specify sample sources, define a section in the config file like this:

  \"sample_sources\" : [{
    \"database_path\" : \"native_features/features.db3\",
    \"id\" : \"Native\",
    \"reference\" : true
  }, {
    \"database_path\" : \"relax_native_features_talaris2014/features.db3\",
    \"id\" : \"talaris2014\",
    \"reference\" : false
  }]

Each sample source needs at least a database_path, id, and reference fields. This is parsed into a data.frame that is passed to each analysis script.\n"

	if(!("sample_sources" %in% names(configuration)) ||
		!is.data.frame(configuration$sample_sources)){
		stop(paste0("ERROR: no sample sources were provided.\n", sample_sources_error_msg))
	}
	sample_sources <- configuration$sample_sources

	if(!("database_path" %in% names(sample_sources)) ||
		any(is.na(sample_sources$database_path))){
		stop(paste0(
			"ERROR: Each sample sources must have a field 'database_path'.\nspecified sample_sources:\n",
			sample_sources, "\n",
			sample_sources_error_msg))
	}
	if(!("id" %in% names(sample_sources)) ||
		any(is.na(sample_sources$reference))){
		stop(paste0(
			"ERROR: Each sample sources must have a field 'id'.\nspecified sample_sources:\n",
			sample_sources, "\n",
			sample_sources_error_msg))
	}
	if(anyDuplicated(sample_sources$id)){
		stop(paste0(
			"ERROR: Each sample source id must be distinct"))
	}
	# sample_source is a synonym for id
	sample_sources$sample_source <- sample_sources$id

	if(!("reference" %in% names(sample_sources)) ||
		any(is.na(sample_sources$reference))){
		stop(paste0(
			"ERROR: Each sample sources must have a true/false field 'reference' indicating if the sample source is used as a reference.\nspecified sample_sources:\n",
			sample_sources, "\n",
			sample_sources_error_msg))
	}

	if("db_cache_size" %in% configuration){
		cache_size <- as.numeric(configuration$db_cache_size)
	} else {
		cache_size <- default_cache_size
	}

	if(!initialize_databases){
		sample_sources$con <- NA
	} else {
		sample_sources$con <- plyr::llply(
			sample_sources$database_path,
			function(database_path){
				con <- NULL
				tryCatch({
					con <- DBI::dbConnect(db_engine, as.character(database_path))
					set_db_cache_size(con, cache_size)
				}, error=function(e){
					stop("Unable to connecto to database '", database_path, "' with error:\n",e, sep="")
				})
				con
			})
	}

	configuration$sample_sources <- sample_sources
	configuration
}

finalize_sample_sources <- function(configuration, verbose=T){
	plyr::a_ply(configuration$sample_sources, 1, function(ss){
		if(verbose){
			cat("Closing connection to database for sample source '", ss$id, "' ... \n", sep="")
		}
		tryCatch({
			DBI::dbDisconnect(ss$con[[1]])
		}, error=function(e){
			cat("ERROR: unable to close database connection to sample source '", ss$id, "' with error:\n", e, sep="")
		})
	})
}

initialize_analysis_scripts <- function(configuration, verbose=F){
	error_message <-  "Each should be a path to a R script that appends a FeaturesAnalysis S4 class to the features_analyses:

   feature_analyses <- c(feature_analyses, new(\"FeaturesAnalysis\",
   id = \"AHdist_chem_type\",
   filename=\"analysis/plots/hbonds/1d_geom/AHdist_chem_type.R\",
   author = \"Matthew O'Meara\",
   brief_description = \"Feature Distribution for AHdist H-bond angles\",
   feature_reporter_dependencies = c(\"StructureFeatures\", \"HBondFeatures\"),
   run=function(self, sample_sources, output_dir, output_formats){
   ...
   })"

	if(!("analysis_scripts" %in% names(configuration)) ||
		!is.character(configuration$analysis_scripts)){
			stop(paste0("ERROR: No analysis scripts were speficied.\n", error_message, "\n"))
	}

	package_scripts_base <- paste(path.package("RosettaFeatures"), sep="/")
	configuration$analysis_scripts <- plyr::llply(
		configuration$analysis_scripts,
		function(analysis_script){
		# parse all the analysis scripts
	if(file.exists(normalizePath(analysis_script, mustWork=F))){
				analysis_script <- normalizePath(analysis_script, mustWork=T)
		} else if(file.exists(paste(package_scripts_base, analysis_script, sep="/"))){
				analysis_script <- paste(package_scripts_base, analysis_script, sep="/")
		} else {
			stop(paste(
				"ERROR: The features analysis script '",
				analysis_script,"' does not exist at either of these locations:\n",
				"\t", normalizePath(analysis_script, mustWork=F), "\n",
				"\t", paste(package_scripts_base, analysis_script, sep="/"), "\n", sep=""))
		}
	})

	configuration
}

initialize_output_formats <- function(configuration, verbose=F){
	if(!("output_formats" %in% names(configuration))){
		cat("WARNING: no output formats specified, using defaults: output_print_pdf, output_csv\n")
		configuration$output_formats <- c('output_print_pdf', 'output_csv')
 	}

	if("add_footer" %in% names(configuration)){
		add_footer <- configuration$add_footer
	} else {
		add_footer <- TRUE
	}

	configuration$output_formats <- get_output_formats(
		configuration$output_formats,
		add_footer)

	if(nrow(configuration$output_formats) == 0){
		stop("ERROR: The output formats specified were not recognized.\n")
	}
	configuration
}

summarize_configuration <- function(configuration, verbose=F){
	if(verbose){
	 	cat(
	 		"Sample Source Comparison:\n",
	 		"  Output Directory <- '", configuration$output_dir, "'\n", sep="")
	 	cat(
	 		"  Output Formats <- ", paste(configuration$output_formats$id, collapse=", "), "\n\n",
	 		sep="")

	 	cat("  Sample Sources:\n")
	 	plyr::a_ply(configuration$sample_sources, 1, function(ss) {
	 		cat("  ", ss$id, " <- ", ss$database_path, "\n", sep="")
	 	})
	 	cat("\n  Analysis_scripts:\n")
	 	cat(" ", paste(configuration$analysis_scripts, sep="", colapse="\n "))
	 	cat("\n")
	}
	return(NULL)
}

parse_analysis_scripts <- function(configuration, verbose=F){

	# This is a vector of FeaturesAnalysis objects that are defined each features analysis script
	feature_analyses <- c()
	num_feature_analyses_before <- 0
	for(analysis_script in configuration$analysis_scripts){

		tryCatch({
			source(analysis_script, local=T, chdir=T)
		}, error=function(e){
			cat(paste(
				"ERROR: Failed to parse the Features Analysis '",
				analysis_script,"' with the following error message:\n", e, sep=""))
		})

		# assign the filename to each feature analysis
		num_new_feature_analyses = length(feature_analyses) - num_feature_analyses_before
		for(i in seq(num_feature_analyses_before+1, length.out=num_new_feature_analyses)){
			feature_analyses[i][[1]]@filename <- analysis_script
		}
		num_feature_analyses_before <- length(feature_analyses)
	}
	feature_analyses
}

#' @export
compare_sample_sources <- function(
	config_filename,
	db_engine=NULL,
	verbose=T,
	dry_run=F
){
	if(verbose){
		cat("Initializing configuration:\n")
	}
	configuration <- load_config_file(config_filename, verbose=verobse)
	configuration <- initialize_output_dir(configuration, verbose)

	db_engine <- initialize_db_engine(db_engine, configuration, verbose=verbose)
	configuration <- initialize_sample_sources(configuration, db_engine, verbose=verbose)
	configuration <- initialize_analysis_scripts(configuration, verbose)
	configuration <- initialize_output_formats(configuration, verbose)

	summarize_configuration(configuration, verbose=verbose)

	feature_analyses <- parse_analysis_scripts(configuration, verbose=verbose)

	# set current working directory the base_dir so that way scripts
	# can reference other scripts in a canonical way
	#setwd(base_dir)

	if(!dry_run){
		for(features_analysis in feature_analyses){
			cat(paste0("Features Analysis: ", features_analysis@id, "\n"))
			cat(paste0("Script path: ", features_analysis@filename, "\n"))

			tryCatch({
				features_analysis@run(
					features_analysis,
					configuration$sample_sources,
					configuration$output_dir,
					configuration$output_formats)
			}, error=function(e){
				cat(paste0(
					"ERROR: Failed to run the Features Analysis '", features_analysis@id, "' ",
					"with the following error message:\n", e))
			})
			cat("\n")
		}
	}
	finalize_sample_sources(configuration, verbose=verbose)
}

#' @export
compare_sample_sources_prepare_SGE_run <- function(
	config_filename,
	run_base,
	db_engine=NULL,
	verbose=T
){

	if(verbose){
		cat("Preparing SGE run to do features analysis:\n")
	}
	initial_configuration <- load_config_file(config_filename, verbose=verobse)
	configuration <- initialize_output_dir(initial_configuration, verbose)
	configuration <- initialize_analysis_scripts(configuration, verbose)

	submit.sh <- paste0("#/bin/bash

# this script is auto-generated by the compare_sample_sources_prepare_SGE_run
# function in the RosettaFeatures R package
#
# This script submits each analysis script as a separate job to a Sun Grid Engine cluster

qsub -t 1-3 -N features_analysis run_job.sh
#qsub -t 1-$(ls -l jobs/ | wc -l) -N features_analysis run_job.sh
")

	run_job.sh <-paste0("#/bin/bash
#$ -S /bin/bash
#$ -cwd

# this script is auto-generated by the compare_sample_sources_prepare_SGE_run
# function in the RosettaFeatures R package
#
# Use the submit.sh to submit this script as a job to a Sun Grid Engine cluster

TASK_INPUT=$(cd jobs && ls | sed -n ${SGE_TASK_ID}p )
cd jobs/$TASK_INPUT
echo \"In run_job.sh\" &> log
echo \"  SGE_TASK_ID: ${SGE_TASK_ID}\" &> log
echo \"  pwd: $(pwd)\" &> log
time /mnt/nfs/home/momeara/opt/bin/Rscript --default-packages=methods,utils analyze_features.R &> log
")

	features_analysis.R <- "#!/usr/bin/env Rscript

# this script is auto-generated by the compare_sample_sources_prepare_SGE_run
# function in the RosettaFeatures R package
#
# It loads the RosettaFeatures package and executes the features analysis

# load the packages that normally are loaded when R starts
library(grid)
library(stats)
library(graphics)
library(grDevices)
library(utils)
library(datasets)
library(methods)
library(base)

library(RosettaFeatures)

compare_sample_sources(
	\"analysis_configuration.json\",
	verbose=T)
"


	dir.create(run_base, recursive=T)
	cat(submit.sh, file=paste(run_base, "submit.sh", sep="/"))
	cat(run_job.sh, file=paste(run_base, "run_job.sh", sep="/"))
	dir.create(paste(run_base, "jobs", sep="/"))
	for(analysis_script in configuration$analysis_scripts){
		job_base_dir <- paste(
			run_base, "jobs", analysis_script %>% stringr::str_replace_all("/", "__"), sep="/")
		dir.create(job_base_dir)
		cat(features_analysis.R, file=paste(job_base_dir, "analyze_features.R", sep="/"))
		job_configuration <- initial_configuration
		job_configuration$analysis_scripts <- c(analysis_script)
		job_configuration %>%
			jsonlite::toJSON(pretty=T) %>%
			cat(file=paste(job_base_dir, "analysis_configuration.json", sep="/"))
	}

	cat("\nTo run benchmark on Sun Grid Engine\n\n  cd ", run_base, " && bash submit.sh\n", sep="")
}

