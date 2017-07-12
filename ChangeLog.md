# Change Log
All notable changes to `libpll` will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.3.2] - 2017-07-12
### Added
 - Optional per-rate category scalers for protein and generic kernels
 - Hardware detection for APPLE builds defaults to assembly code
### Fixed
 - Improved fix for negative p-matrix values.
 - Derivatives computation for Lewis/Felsenstein ascertainment bias correction
 - set_tipchars() for ascertainment bias correction with non-DNA sequences
 - Excessive memory allocation when compressing site patterns
 - Issue with PHYLIP parsing when header ends with CRLF

## [0.3.1] - 2017-05-17
### Added
 - Checks for older versions of clang and gcc to use assembly instructions
   for cpu features detection
 - Include guards for pll.h
### Fixed
 - Correct updating of padded eigen-decomposition arrays for models with a
   number of states not being a power of two
 - Changed to the usage of builtin functions for cpu features detection
 - Check for x86intrin.h

## [0.3.0] - 2017-05-15
### Added
 - Run-time detection of cpu features
 - Vectorized (AVX) computation of 20-state transition probability matrices
 - Faster tip-inner kernels for 20-state models
 - Improved AVX vectorization of derivatives
 - Faster PHYLIP parser
 - vectorized scaling for 20-state and arbitrary-state models
 - AVX2 vectorizations for partials, likelihood and derivatives
 - Unweighted parsimony functions including SSE, AVX and AVX2 vectorizations
 - Randomized stepwise addition
 - Portable functions for parsing trees from a C-string
 - Optional per-rate category scalers to prevent numerical underflows on large
   trees
 - Setting of identity matrix if all exponentiations of eigenvalues multiplied
   by branch length and rate are approximately equal to one
 - Re-entrant cross-platform pseudo-random number generator
 - Wrapper tree structures
 - Custom exporting of tree structures using a callback function
 - Support for median category rates in discrete gamma model
 
### Fixed
 - Derivatives computation
 - Parsing of branch lengths in newick trees
 - Invariant sites computation
 - Multiplication of log-likelihood with pattern weight after scaling term
 - Added destructors for eliminating memory leaks when tree parsing fails
 - Sumtable computation when having multiple substitution matrices
 - Ascertainment bias computation
 - Per-site log-likelihood computation
 - Uninitialized values in testing framework



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
