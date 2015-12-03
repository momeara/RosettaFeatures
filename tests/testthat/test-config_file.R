context("Configuration File")


test_that("A simple configuration file is parsed correctly", {
	verbose = T
	config_filename <- "test-config_file__analysis_configuration.json"
	configuration <- load_config_file(config_filename, verbose=verobse)
	configuration <- initialize_output_dir(configuration, verbose)
	
	db_engine <- initialize_db_engine(NULL, configuration, verbose=verbose)
	configuration <- initialize_sample_sources(configuration, db_engine, verbose=verbose)
	configuration <- initialize_analysis_scripts(configuration, verbose)
	configuration <- initialize_output_formats(configuration, verbose)
	
	
	expect_equal(configuration$output_dir, "build/native_vs_talaris2014")
})
