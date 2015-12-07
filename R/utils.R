#!/usr/bin/env Rscript
# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


#' create date code of the form YYMMDD
#' inputs:
#'    d (optional): date code in the format of Sys.Date() for which to generate the date code, defaulting to 'today'
#' @export
date_code <- function(d=NA){
    # reference http://www.r-cookbook.com/node/17
    if(is.na(d)) d <- Sys.Date()
    pattern <- '20([[:digit:]]{2})-([[:digit:]]{2})-([[:digit:]]{2})'
    paste(
        sub(pattern, '\\1', d), sub(pattern, '\\2', d), sub(pattern, '\\3', d),
        sep="")
}


#' @importFrom magrittr %>%
NULL
