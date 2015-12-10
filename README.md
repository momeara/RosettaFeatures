# Features Analysis Benchmark Protocol Capture

## Overview
High-resolution crystal structures are assumed to be
thermodynamically stable, which means that each structure should
occupy a minimum of a high quality energy function. In particular, in
aggregate, Rosetta should reproduce local structure features observed
in the high quality crystal structures. Therefore if we relax native
crystal structures and observe systematic displacement of local
geometric features (specific types of bonds are lost, polar contacts
become too close, etc.) this indicates that when new structures are
predicted, they will likely not be stable in the lab.

In our 2015 H-bond paper, we observed that the interaction between
related terms in the energy function create areas of double counting
for certain types of interactions or contexts that may be maybe
difficult to observe in full structure benchmarks. By plotting and
looking closely at distributions of local features, it is possible to
gain a deeper understanding of the full energy that arises from the
complex interaction between the range of active terms and the implicit
kinematic constraints.


## Literature references

Combined Covalent-Electrostatic Model of Hydrogen Bonding Improves Structure Prediction with Rosetta Matthew J. O’Meara, Andrew Leaver-Fay, Michael D. Tyka, Amelie Stein, Kevin Houlihan, Frank DiMaio, Philip Bradley, Tanja Kortemme, David Baker, Jack Snoeyink, and Brian Kuhlman J. Chem. Theory Comput., 2015, 11 (2), pp 609–622 DOI: 10.1021/ct500864r

Scientific benchmarks for guiding macromolecular energy function improvement Andrew Leaver-Fay, Matthew J O’Meara, Mike Tyka, Ron Jacak, Yifan Song, Elizabeth H Kellogg, James Thompson, Ian W Davis, Roland A Pache, Sergey Lyskov, Jeffrey J Gray, Tanja Kortemme, Jane S Richardson, James J Havranek, Jack Snoeyink, David Baker, Brian Kuhlman Methods in enzymology 523 2013


## Usage:

To install this package, in R:

    if (packageVersion("devtools") < 1.6) {
      install.packages("devtools")
    }
    devtools::install_github("momeara/RosettaFeatures")

Generate features databases following the features_benchmark protocol capture

    https://github.com/RosettaCommons/demos/tree/master/protocol_capture/features_benchmark/README.md


Then to report features, in R:

    library(RosettaFeatures)
    compare_sample_sources(
      config_filename="analysis_configuration.json")

Where the `analysis_configuration.json` looks like:

    {
      "output_dir" : "native_vs_relax_native",

      "sample_sources" : [{
        "database_path" : "native_features/features.db3",
        "id" : "Native",
        "reference" : true
      }, {
        "database_path" : "relax_native_features/features.db3",
        "id" : "talaris2014",
        "reference" : false
      }],
    
      "analysis_scripts" : [
        "scripts/analysis/plots/EXAMPLE_PLOT.R"
      ],

      "output_formats" : [
        "output_print_pdf"
      ]
    }


 

       




   
