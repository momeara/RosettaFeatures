#!/usr/bin/env Rscript
# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

set_db_cache_size <- function(con, cache_size){
	if(is.null(cache_size)){
		stop("ERROR: unable to set database cache size because the cache_size is null.")
	}
	dbGetQuery(con,
		paste("PRAGMA cache_size=",as.integer(cache_size),";",sep=""))
}


query_sample_source <- function(
	sample_source,
	sele,
	sele_args_frame = NULL,
	cache_size=database_configuration$db_cache_size,
	char_as_factor=T
	){
	tryCatch(sele,error=function(e){
		cat("ERROR: The select statement is not defined.\n")
	})

	#Get the features from the sample_source
	tryCatch(c(sample_source),error=function(e){
		cat("ERROR: The specified sample source is not defined.\n")
	})

	cat("loading:", as.character(sample_source$sample_source), "... ")
	if( is.na(sample_source$sample_source[1]) ){
		stop("Specified sample source is not defined")
	}
	con <- dbConnect(engine, as.character(sample_source$fname))
	set_db_cache_size(con, cache_size);

	#Allow select statements to be prefaced with arbitrary statements.
	#This allows the creation of temporary tables, indices, etc.
	sele_split <- paste(strsplit(sele, ";\\W*", perl=TRUE)[[1]], ";", sep="")
	l_ply(sele_split[-length(sele_split)], function(sele){
			dbGetQuery(con, sele)
	})

	last_stmt <- sele_split[length(sele_split)] #Real statement for data
	if (! is.null(sele_args_frame)){
		features <- dbGetPreparedQuery(con, sele, bind.data=sele_args_frame)
	}
	else {
		features <- dbGetQuery(con, last_stmt)
	}

	dbDisconnect(con)
	#cat(as.character(round(timing[3], 2)),"s\n")

	if(nrow(features)==0){
		cat("WARNING: The following query returned no rows:\n")
		cat(sele)
	}
	if(char_as_factor){
		for(col in names(features)){
			if(is.character(features[,col])){
				features[,col] <- factor(features[,col])
			}
		}
	}
	features
}

query_sample_sources <- function(
	sample_sources,
	sele,
	sele_args_frame = NULL,
	cache_size=database_configuration$db_cache_size,
	char_as_factor=T
	){
	tryCatch(sele,error=function(e){
		cat("ERROR: The select statement is not defined.\n")
	})
	features <- ddply(sample_sources, c("sample_source"), function(ss){
		tryCatch(c(ss),error=function(e){
			cat("ERROR: The specified sample source is not defined.\n")
		})

		cat("loading:", as.character(ss$sample_source), "... ")
		if( is.na(ss$sample_source[1]) ){
			stop("Specified sample source is not defined")
		}
		timing <- system.time({
			con <- dbConnect(engine, as.character(ss$fname))
			set_db_cache_size(con, cache_size);

			#Allow select statements to be prefaced with arbitrary statements.
			#This allows the creation of temporary tables, indices, etc.
			sele_split <- paste(strsplit(sele, ";\\W*", perl=TRUE)[[1]], ";", sep="")
			l_ply(sele_split[-length(sele_split)], function(sele){
					dbGetQuery(con, sele)
			})

			last_stmt <- sele_split[length(sele_split)] #Real statement for data
			if (! is.null(sele_args_frame)){
				df <- dbGetPreparedQuery(con, sele, bind.data=sele_args_frame)
			}
			else {
				df <- dbGetQuery(con, last_stmt)
			}

			dbDisconnect(con)
		})
		cat(as.character(round(timing[3], 2)),"s\n")
		df
	})
	if(nrow(features)==0){
		cat("WARNING: The following query returned no rows:\n")
		cat(sele)
	}
	if(char_as_factor){
		for(col in names(features)){
			if(is.character(features[,col])){
				features[,col] <- factor(features[,col])
			}
		}
	}
	features
}

query_sample_sources_against_ref <- function(
	sample_sources,
	sele,
	sele_args_frame = NULL,
	cache_size=database_configuration$db_cache_size,
	char_as_factor=T
	){
	tryCatch(sele,error=function(e){
		cat("ERROR: The select statement ", sele, " is not defined.\n")
	})

	if(nrow(sample_sources) < 2) {
		stop(paste("Please provide 2 or more sample sources.

sample_sources provided:
	'", paste(sample_sources$sample_source_ids, collapse= "',\n\t'"), "'

This function will execute the query for all but the first sample source
	The first sample source will be available as a database named 'ref'
	The second sample source will be available as a database named 'new'.

In the returned data.frame the there will be the following columns:
	'ref_sample_source' -> id of ref sample_source
	'new_sample_source' -> id of new sample_source
", sep=""))
	}

	ref_ss <- sample_sources[1,]
	con <- dbConnect(engine)
	set_db_cache_size(con, cache_size);
	dbGetQuery(con, paste("ATTACH DATABASE '", ref_ss$fname, "' AS ref;", sep=""))

	features <- ddply(sample_sources[seq(2,nrow(sample_sources)),], c("sample_source"), function(ss){
		tryCatch(c(ss),error=function(e){
			cat("ERROR: The specified sample source is not defined.\n")
		})
		cat("loading: ref:", as.character(ref_ss$sample_source),
				" new:", as.character(ss$sample_source), " ... ", sep="")
		if( is.na(ss$sample_source[1]) ){
			stop("Specified sample source is not defined")
		}
		timing <- system.time({
			dbGetQuery(con, paste("ATTACH DATABASE '", ss$fname, "' AS new;", sep=""))

			#Allow select statements to be prefaced with arbitrary statements.
			#This allows the creation of temporary tables, indices, etc.
			sele_split <- paste(strsplit(sele, ";\\W*", perl=TRUE)[[1]], ";", sep="")
			l_ply(sele_split[-length(sele_split)], function(sele){
				dbGetQuery(con, sele)
			})
			last_stmt <- sele_split[length(sele_split)] #Real statement for data
			if(! is.null(sele_args_frame)){
				df <- dbGetPreparedQuery(con, last_stmt, bind.data=sele_args_frame)
			}
			else{
				df <- dbGetQuery(con, last_stmt)
			}
		})
		dbGetQuery(con, "DETACH DATABASE new;")
		cat(as.character(round(timing[3],2)),"s\n")
		df
	})
	dbDisconnect(con)

	if(nrow(features)==0){
		cat("WARNING: The following query returned no rows:\n")
		cat(sele)
		return(features)
	}
	if(char_as_factor){
		for(col in names(features)){
			if(is.character(features[,col])){
				features[,col] <- factor(features[,col])
			}
		}
	}
	data.frame(
		ref_sample_source = factor(ref_ss$sample_source[1]),
		new_sample_source = factor(features$sample_source),
		sample_source = factor(features$sample_source),
		subset(features, select= -sample_source))
}
