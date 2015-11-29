
library(readr)

output_formats <- read_tsv("output_formats.tsv")

devtools::use_data(output_formats, internal=TRUE)
