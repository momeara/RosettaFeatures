


check_basic_input <- function(
	data_sources,
	dry_run,
	config_filename){
	if(dry_run & length(data_sources) == 0 & is.null(config_filename)){
		cat(
			"ERROR: No configuration file or sample_source databases was supplied",
			"",
			"##################################################",
			"NAME",
			"    compare_sample_source.R",
			"",
			"SYNOPSIS",
			"    ./compare_sample_sources.R [OPTIONS] --script <analysis_script> features_<ss_id1>.db3 [features_<ss_id2>.db3 ...]",
			"    ./compare_sample_sources.R [OPTIONS] --analysis_dir <analysis_dir> features_<ss_id1>.db3 [features_<ss_id2>.db3 ...]",
			"",
			"    From within an R session:",
			"       source(\"compare_sample_sources_iscript.R\") # This will re-run the last './compare_sample_sources.R' interactively",
			"",
			"DESCRIPTION",
			"    Compare structural features coming from different sample sources.",
			"",
			"    See ./compare_sample_sources.R --help for available options.",
			"    See https://wiki.rosettacommons.org/index.php/FeaturesScientificBenchmark for application documentation.",
			"    To cite: contact mattjomeara@gmail.com as the work is in progress.",
			"",
			"EXAMPLE",
			"    To compare how Rosetta distorts native structures, extract feature databases for a set of structures from the pdb",
			"    and the same set of structures optimized with your favorite prediction protocol. Assume the resulting feature databases",
			"    are stored in 'features_natives.db3' and 'features_rosetta.db3'. To compare the length of hydrogen bonds",
			"    conditional on the donor and acceptor chemical types between the two sample sources, run:",
			"",
			"    ./compare_sample_sources.R --script scripts/analysis/plots/hbonds/AHdist_chem_type.R features_rosetta.db3 features_rosetta.db3",
			"",
			"AUTHOR",
			"    Matthew O'Meara (mattjomeara@gmail.com)\n", sep="\n")
	  quit()
	}
}


#' Compare sample sources
#' @param script path to a single analysis script
#' @param analysis_dir directory where the analysis scripts are located. The supplied directory is searched recursively for files of the form \"*.R\"
#' @param output_dir directory where the output plots and statistics will be generated
#' @param config_filename path to json file with additional configuration information
#' @param verbose print extra output
#' @param dry_run prepare the analyses but do not run them
#' @param db_cache_size number of 1k pages of cache to use for database queries
#' @param add_footer add footer to plots saying the analysis script and run date
#' @param ncores run in parallel with specified number of cores
#' @param general_adjust_kernel multiplicative factor in the kernel bandwith for density estimation
#' @param out_disable_ss_names disable concatonating the sample source names for output.  Be sure to combine this with out_dir or all plots will go to build
#' @param generate_website add footer to plots saying the analysis script and run date

compare_sample_sources(
	script=NULL,
	analysis_dir=NULL,
	output_dir="build",
	config_filename=NULL,
	verbose=F,
	dry_run=F,
	db_cache_size=10000,
	add_footer=T,
	ncores=1,
	general_adjust_kernel=1`
	generate_website=F,
	out_disable_ss_names
){		


}