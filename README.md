# libpll

[![Build Status](https://travis-ci.org/xflouris/libpll.svg?branch=dev)](https://travis-ci.org/xflouris/libpll)
[![License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

## Introduction

The aim of this project is to implement a versatile high-performance software
library for phylogenetic analysis. The library should serve as a lower-level
interface of PLL (Flouri et al. 2015) and should have the following
properties:

* open source code with an appropriate open source license.
* 64-bit multi-threaded design that handles very large datasets.
* easy to use and well-documented.
* SIMD implementations of time-consuming parts.
* as fast or faster likelihood computations than RAxML (Stamatakis 2014).
* fast implementation of the site repeats algorithm (Kobert 2017).
* functions for tree visualization.
* bindings for Python.
* generic and clean design.
* Linux, Mac, and Microsoft Windows compatibility.

## Compilation instructions

Currently, `libpll` requires that [GNU Bison](http://www.gnu.org/software/bison/)
and [Flex](http://flex.sourceforge.net/) are installed on the target system. On
a Debian-based Linux system, the two packages can be installed using the command

`apt-get install flex bison`

The library also requires that a GNU system is available as it uses several
functions (e.g. `asprintf`) which are not present in the POSIX standard.
This, however will change in the future in order to have a more portable
and cross-platform library.

The library can be compiled using either of the following two ways.

**Cloning the repo** Clone the repo and bild the executable and documentation
using the following commands.

```bash
git clone https://github.com/xflouris/libpll.git
cd libpll
./autogen.sh
./configure
make
make install    # as root, otherwise run: sudo make install
```

When using the cloned repository version, you will also need
[autoconf](https://www.gnu.org/software/autoconf/autoconf.html),
[automake](https://www.gnu.org/software/automake/) and
[libtool](https://www.gnu.org/software/libtool/) installed. On a Debian-based
Linux system, the packages can be installed using the command

```bash
sudo apt-get install autotools-dev autoconf libtool
```

The library will be installed on the operating system's standard paths.  For
some GNU/Linux distributions it might be necessary to add that standard path
(typically `/usr/local/lib`) to `/etc/ld.so.conf` and run `ldconfig`.

Microsoft Windows compatibility was tested with a cross-compiler and seems to
work out-of-the-box using [MingW](http://www.mingw.org/).

## Available functionality

libpll currently implements the General Time Reversible (GTR) model (Tavare
1986) which can be used for nucleotide and amino acid data. It supports models
of variable rates among sites, the Inv+&Gamma; (Gu et al. 1995) and has
functions for computing the discretized rate categories for the gamma model
(Yang 1994). Furthermore, it supports several methods for
[ascertainment bias correction](https://github.com/xflouris/libpll/wiki/Ascertainment-bias-correction)
(Kuhner et al. 2000, McGill et al. 2013, Lewis 2011, Leach&eacute; et al.
2015). Additional functionality includes tree visualization, functions for
parsimony (minimum mutation cost) calculation and ancestral state
reconstruction using Sankoff's method (Sankoff 1975, Sankof and Rousseau 1975).
The functions for computing partials, evaluating the log-likelihood and
updating transition probability matrices are vectorized using both SSE3, AVX and
AVX2 instruction sets.

## Documentation

Please refer to the [wiki page](https://github.com/xflouris/libpll/wiki).

## Usage examples

Please refer to the [wiki page](https://github.com/xflouris/libpll/wiki) and/or
the [examples directory](https://github.com/xflouris/libpll/tree/master/examples).

## libpll license and third party licenses

The libpll code is currently licensed under the
[GNU Affero General Public License version 3](http://www.gnu.org/licenses/agpl-3.0.en.html).
Please see LICENSE.txt for details.

libpll includes code from several other projects. We would like to thank the
authors for making their source code available.

libpll includes code from GNU Compiler Collection distributed under the GNU
General Public License.

## Code

The code is written in C with some parts written using in-line assembler and intrinsic functions.

|     File                    | Description                                                                                  |
| --------------------------- | -------------------------------------------------------------------------------------------- |
| **compress.c**              | Functions for compressing alignment into site patterns.                                      |
| **core_derivatives_avx2.c** | AVX2 vectorized core functions for computing derivatives of the likelihood function.         |
| **core_derivatives_avx.c**  | AVX vectorized core functions for computing derivatives of the likelihood function.          |
| **core_derivatives.c**      | Core functions for computing derivatives of the likelihood function.                         |
| **core_derivatives_sse.c**  | SSE vectorized core functions for computing derivatives of the likelihood function.          |
| **core_likelihood_avx2.c**  | AVX2 vectorized core functions for computing the log-likelihood.                             |
| **core_likelihood_avx.c**   | AVX vectorized core functions for computing the log-likelihood.                              |
| **core_likelihood.c**       | Core functions for computing the log-likelihood, that do not require partition instances.    |
| **core_likelihood_sse.c**   | SSE vectorized core functions for computing the log-likelihood.                              |
| **core_partials_avx2.c**    | AVX2 vectorized core functions for updating vectors of conditional probabilities (partials). |
| **core_partials_avx.c**     | AVX vectorized core functions for updating vectors of conditional probabilities (partials).  |
| **core_partials.c**         | Core functions for updating vectors of conditional probabilities (partials).                 |
| **core_partials_sse.c**     | SSE vectorized core functions for updating vectors of conditional probabilities (partials).  |
| **core_pmatrix_avx2.c**     | AVX2 vectorized core functions for updating transition probability matrices.                 |
| **core_pmatrix_avx.c**      | AVX vectorized core functions for updating transition probability matrices.                  |
| **core_pmatrix.c**          | Core functions for updating transition probability matrices.                                 |
| **core_pmatrix_sse.c**      | SSE vectorized core functions for updating transition probability matrices.                  |
| **derivatives.c**           | Functions for computing derivatives of the likelihood function.                              |
| **fasta.c**                 | Functions for parsing FASTA files.                                                           |
| **fast_parsimony_avx2.c**   | AVX2 fast unweighted parsimony functions.                                                    |
| **fast_parsimony_avx.c**    | AVX fast unweighted parsimony functions.                                                     |
| **fast_parsimony.c**        | Non-vectorized fast unweighted parsimony functions.                                          |
| **fast_parsimony_sse.c**    | SSE fast unweighted parsimony functions.                                                     |
| **gamma.c**                 | Functions related to Gamma (&Gamma;) function and distribution.                              |
| **hardware.c**              | Hardware detection functions.                                                                |
| **lex_rtree.l**             | Lexical analyzer for parsing newick rooted trees.                                            |
| **lex_utree.l**             | Lexical analyzer for parsing newick unrooted trees.                                          |
| **likelihood.c**            | Functions ofr computing the log-likelihood of a tree given a partition instance.             |
| **list.c**                  | (Doubly) Linked-list implementations.                                                        |
| **maps.c**                  | Character mapping arrays for converting sequences to the internal representation.            |
| **models.c**                | Model parameters related functions.                                                          |
| **output.c**                | Functions for output in terminal (i.e. conditional likelihood arrays, probability matrices). |
| **parse_rtree.y**           | Functions for parsing rooted trees in newick format.                                         |
| **parse_utree.y**           | Functions for parsing unrooted trees in newick format.                                       |
| **parsimony.c**             | Parsimony functions.                                                                         |
| **partials.c**              | Functions for updating vectors of conditional probabilities (partials).                      |
| **phylip.c**                | Functions for parsing phylip files.                                                          |
| **pll.c**                   | Functions for setting PLL partitions (instances).                                            |
| **random.c**                | Re-entrant multi-platform pseudo-random number generator.                                    |
| **rtree.c**                 | Rooted tree manipulation functions.                                                          |
| **utree.c**                 | Unrooted tree manipulation functions.                                                        |
| **utree_moves.c**           | Functions for topological rearrangements on unrooted trees.                                  |
| **utree_svg.c**             | Functions for SVG visualization of unrooted trees.                                           |

## Bugs

The source code in the master branch is thoroughly tested before commits.
However, mistakes may happen. All bug reports are highly appreciated. You may
submit a bug report here on GitHub as an issue, or you could send an email to
[t.flouris@ucl.ac.uk](mailto:t.flouris@ucl.ac.uk).

## libpll core team

* Tom&aacute;&scaron; Flouri
* Diego Darriba
* Kassian Kobert
* Mark T. Holder
* Alexey Kozlov
* Alexandros Stamatakis

## Acknowledgements

Special thanks to the following people for patches and suggestions:

* Frederick Matsen
* Ben Redelings
* Andreas Tille
* Ziheng Yang

## Contributing to libpll

Please read the section [Contributing to `libpll`](https://github.com/xflouris/libpll/wiki#contributing-to-libpll)
of the [wiki](https://github.com/xflouris/libpll/wiki).

## References

* Flouri T., Izquierdo-Carrasco F., Darriba D., Aberer AJ, Nguyen LT, Minh BQ, von Haeseler A., Stamatakis A. (2015)
**The Phylogenetic Likelihood Library.**
*Systematic Biology*, 64(2): 356-362.
doi:[10.1093/sysbio/syu084](http://dx.doi.org/10.1093/sysbio/syu084)

* Gu X., Fu YX, Li WH. (1995)
**Maximum Likelihood Estimation of the Heterogeneity of Substitution Rate among Nucleotide Sites.**
*Molecular Biology and Evolution*, 12(4): 546-557.

* Kobert K., Stamatakis A., Flouri T. (2017)
**Efficient detection of repeating sites to accelerate phylogenetic likelihood calculations.**
*Systematic Biology*, 66(2): 205-217.
doi:[10.1093/sysbio/syw075](http://dx.doi.org/10.1093/sysbio/syw075)

* Leach&eacute; AL, Banbury LB, Felsenstein J., de Oca ANM, Stamatakis A. (2015)
**Short Tree, Long Tree, Right Tree, Wrong Tree: New Acquisition Bias Corrections for Inferring SNP Phylogenies.**
*Systematic Biology*, 64(6): 1032-1047.
doi:[10.1093/sysbio/syv053](http://dx.doi.org/10.1093/sysbio/syv053)

* Lewis LO. (2001)
**A Likelihood Approach to Estimating Phylogeny from Discrete Morphological Character Data.**
*Systematic Biology*, 50(6): 913-925.
doi:[10.1080/106351501753462876](http://dx.doi.org/10.1080/106351501753462876)

* Sankoff D. (1975)
**Minimal Mutation Trees of Sequences.**
*SIAM Journal on Applied Mathematics*, 28(1): 35-42.
doi:[10.1137/0128004](http://dx.doi.org/10.1137/0128004)

* Sankoff D, Rousseau P. (1975)
**Locating the Vertices of a Steiner Tree in Arbitrary Metric Space.**
*Mathematical Programming*, 9: 240-246.
doi:[10.1007/BF01681346](http://dx.doi.org/10.1007/BF01681346)

* Stamatakis A. (2014)
**RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies.**
*Bioinformatics*, 30(9): 1312-1313.
doi:[10.1093/bioinformatics/btu033](http://dx.doi.org/10.1093/bioinformatics/btu033)

* Tavar&eacute; S. (1986)
**Some probabilistic and statistical problems in the analysis of DNA sequences.**
*American Mathematical Sciety: Lectures on Mathematics in the Life Sciences*, 17: 57-86.

* Yang Z. (2014)
**Maximum likelihood phylogenetic estimation from dna sequences with variable rates over sites: Approximate methods.**
*Journal of Molecular Evolution*, 39(3): 306-314.
doi:[10.1007/BF00160154](http://dx.doi.org/10.1007/BF00160154)
