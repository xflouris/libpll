# Change Log
All notable changes to `libpll` will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.2.0] - 2016-09-09
### Added
 - Methods for ascertainment bias correction (Lewis, Felsenstein, Stamatakis)
 - Functions for performing NNI and SPR moves on unrooted trees
 - Iterator function (pll_utree_every)
 - Vectorized (AVX) computation of derivatives and sumtable
 - restructure of source code such that higher-level functions call lower level
   (core) functions
 - Vectorized (AVX and SSE) functions for log-likelihood computation on rooted
   and unrooted trees
 - Implementation of Sankoff's dynamic programming parsimony algorithm, which
   allows for arbitrary scoring matrices and facilitates ancestral state
   reconstruction
 - Faster vectorization scheme for computing partials under an arbitrary number
   of states
 - SSE vectorization for computing partials
 - Faster lookup table constructionf or tip-tip case (PLL_ATTRIB_PATTERN_TIP)
 - More error codes
 - AVX and SSE vectorized functions for updating transition probability matrices
 - SVG visualization for unrooted trees
 - man pages

### Fixed
 - When the number of states is not a power of 2, the log-likelihood results to
   -inf if PLL_ATTRIB_PATTERN_TIP is used
 - Missing memory initialization when using site pattern compression
 - Erroneous charmap usage in pattern compression
 - Removed likelihood computation from the function computing derivatives
 - bug in unrooted tree clone function

## [0.1.0] - 2016-05-25
### Added
 - AVX implementation for updating partials
 - Core functions
 - Selection between character states and cond. probability vectors for tips
 - Compatibility with MS Windows
 - Added basic PHYLIP parsing
 - Functions for obtaining first/second derivatives to approximate LH function
 - Trimming of line-feed/carriage-return when reading FASTA files

### Changed
 - Selection of a specific rate matrix selection for each rate category

### Fixed
 - Newick parsers handle arbitrary string literals as taxon names
