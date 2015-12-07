# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#Extending ggplot2
#http://rstudio-pubs-static.s3.amazonaws.com/108934_8537676801dd4548a96f6451bae01e94.html
GeomIndicator <- ggplot2::ggproto(
	"GeomIndicator", ggplot2::Geom,
	required_aes = c("indicator"),
	default_aes = ggplot2::aes(
		colour = "black",
		xpos="right",
		ypos="top",
		xjust=NULL,
		yjust=NULL,
		size = 5,
		group=1,
		angle = 0,
		alpha = NA,
		family = "",
		fontface = 1,
		lineheight = 1.2),

	draw_group = function(
		data, scales, coordinates,
		parse = FALSE, na.rm = FALSE, check_overlap = FALSE){

		indicator <- data$indicator[1]
		if(is.na(indicator) || is.null(indicator)){
			return(grid::nullGrob())
		}

		if("xpos" %in% names(data)){
			if(data$xpos[1] == "left"){
				xpos <- .07
			} else if( data$xpos[1] == "center"){
				xpos <- .5
			} else if( data$xpos[1] == "right"){
				xpos <- .97
			} else if( is.numeric(data$xpos[1]) && 0 <= data$xpos[1] && data$xpos[1] <= 1){
				xpos <- data$xpos[1]
			} else {
				stop(paste("In geom_indicator(), unrecognized value xpos=\"", data$xpos[1],"\". Please use 'left', 'right' or 'center', or a value from 0 to 1.", sep=""))
			}
		} else {
			xpos <- .97
		}

		if("ypos" %in% names(data)){
			if(data$ypos[1] == "top"){
				ypos <- .97
			} else if( data$ypos[1] == "center"){
				ypos <- .5
			} else if( data$ypos[1] == "bottom"){
				ypos <- .03
			} else if( is.numeric(data$ypos[1]) && 0L <= data$ypos[1] && data$ypos[1] <= 1L){
				ypos <- data$ypos[1]
			} else {
				stop(paste("In geom_indicator(), unrecognized value ypos=\"", data$ypos[1],"\". Please use 'top', 'bottom' or 'center', or a value from 0 to 1.", sep=""))
			}
		} else {
			ypos <- .97
		}

		if(!is.null(data$xjust[1])){
			if(data$xjust[1] %in% c("left", "center", "right")){
				xjust <- data$xjust[1]
			} else {
				stop(paste("In geom_indicator(), unrecognized value xjust=\"", data$xjust[1],"\". Please use 'left', 'right' or 'center'.", sep=""))
			}
		} else {
			if(xpos < 1/3){
				xjust <- "left"
			} else if(xpos >= 1/3 && xpos < 2/3){
				xjust <- "center"
			} else {
				xjust <- "right"
			}
		}

		if(!is.null(data$yjust[1])){
			if(data$yjust[1] %in% c("top", "center", "bottom")){
				yjust <- data$yjust[1]
			} else {
				stop(paste("In geom_indicator(), unrecognized value yjust=\"", data$yjust[1],"\". Please use 'top', 'bottom' or 'center'.", sep=""))
			}
		} else {
			if(ypos < 1/3){
				yjust <- "bottom"
			} else if(ypos >= 1/3 && ypos < 2/3){
				yjust <- "center"
			} else {
				yjust <- "top"
			}
		}


		if (parse) {
			# this is what the new geom_text does
			indicator_display_value <- parse(text = as.character(data$indicator[1]))
		} else {
			if(is.character(indicator)){
				indicator_display_value <- indicator
			} else {
				indicator_display_value <- prettyNum(data$indicator[1], big.mark=",")
			}
		}

		size <- data$size[1]
		level <- data$group[1] - 1

		textGrob(
			indicator_display_value,
			unit(xpos, "npc"),
			unit(ypos, "npc") - unit(level, "line"),
			just=c(xjust, yjust),
			gp=gpar(
				col=alpha(data$colour[1], data$alpha[1]),
				fontsize=size*12/5,
				fontfamily=data$family[1],
				fontface=data$fontface[1],
				lineheight=data$lineheight[1]),
			check.overlap = check_overlap)})

#' @export
geom_indicator <- function (
	mapping = NULL,
	data = NULL,
	stat = "identity",
	position = "identity",
	parse = FALSE,
	...,
	nudge_x = 0,
	nudge_y = 0,
	check_overlap = FALSE,
	na.rm = FALSE,
	show.legend = NA,
	inherit.aes = TRUE
) {

	if (!missing(nudge_x) || !missing(nudge_y)) {
		if (!missing(position)) {
			stop("Specify either `position` or `nudge_x`/`nudge_y`", call. = FALSE)
		}

		position <- position_nudge(nudge_x, nudge_y)
	}

	ggplot2::layer(
		data = data,
		mapping = mapping,
		stat = stat,
		geom = GeomIndicator,
		position = position,
		inherit.aes = inherit.aes,
		show.legend=F,
		params=list(
			parse = parse,
			check_overlap = check_overlap,
			na.rm = na.rm,
			...
		)
	)
}
