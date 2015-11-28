context("Configuration File")


test_that("A simple configuration file is parsed correctly", {
	configuration <- initialize_config_file("test-config_file__analysis_configuration.json")
        expect_equal(configuration$output_dir, "build/native_vs_talaris2014")
})
