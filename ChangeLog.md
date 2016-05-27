# Change Log
All notable changes to `libpll` will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.1.0] - 2016-05-25
### Added
 - AVX implementation for updating partials
 - Core functions
 - Selection between character states and cond. probability vectors for tips
 - Compatibility with MS Windows
 - Added basic PHYLIP parsing
 - Functions for obtaining first/second derivatives to approximate LH function
 - Trimming of line-feed/carriage-return when reading FASTA files
 - Ascertainment bias correction (Lewis, Felsenstein and Stamatakis algorithms) for edge and root likelihood computation

### Changed
 - Selection of a specific rate matrix selection for each rate category

### Fixed
 - Newick parsers handle arbitrary string literals as taxon names
