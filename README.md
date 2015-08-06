# libpll

[![Build Status](https://magnum.travis-ci.com/xflouris/libpll.svg?token=rjft2y6GBHow4SDyjuoy&branch=master)](https://magnum.travis-ci.com/xflouris/libpll)
[![License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

## Introduction

The aim of this project is to implement a versatile high-performance software
library for phylogenetic analysis. The library should serve as a lower-level
interface of the PLL (Flouri et al. 2014) and should have the following
properties:

* have an open source code with an appropriate open source license.
* 64-bit multi-threaded design that handles very large datasets.
* easy to use and well-documented.
* SIMD implementations of time-consuming parts.
* Comparably fast likelihood computations as RAxML (Stamatakis 2014).
* generic and clean design.

## Compilation instructions


Currently, `libpll` requires that [GNU Bison](http://www.gnu.org/software/bison/)
and [Flex](http://flex.sourceforge.net/) are installed on the target system. On
a Debian-based Linux system, the two packages can be installed using the command

`apt-get install flex bison`

The library also requires that a GNU system is available as it uses several
functions (e.g. `asprintf`) which are not present in the POSIX standard.
This, however might change in the future in order to have a more portable
and cross-platform library.

The library can be compiled using the included Makefile:

`make`

## Available functionality

libpll currently implements the General Time Reversible (GTR) model (Tavare
1986) which can be used for nucleotide and amino acid data. It supports models
of variable rates among sites, the Inv+&Gamma; (Gu et al. 1995) and has
functions for computing the discretized rate categories for the gamma model
(Yang 1994).

Below is a list of available functions in the current version.

### Partition (instance) manipulation

* `pll_partition_t * pll_partition_create(unsigned int tips, unsigned int clv_buffers, unsigned int states, unsigned int sites, unsigned int rate_matrices, unsigned int prob_matrices, unsigned int rate_cats, unsigned int scale_buffers, int attributes);`
* `void pll_partition_destroy(pll_partition_t * partition);`

### Linked lists

* `int pll_dlist_append(pll_dlist_t ** dlist, void * data);`
* `int pll_dlist_prepend(pll_dlist_t ** dlist, void * data);`
* `int pll_dlist_remove(pll_dlist_t ** dlist, void * data);`

### Models setup

* `void pll_set_subst_params(pll_partition_t * partition, unsigned int params_index, const double * params);`
* `void pll_set_frequencies(pll_partition_t * partition, unsigned int params_index, const double * frequencies);`
* `void pll_set_category_rates(pll_partition_t * partition, const double * rates);`
* `void pll_update_prob_matrices(pll_partition_t * partition, unsigned int params_index, unsigned int * matrix_indices, double * branch_lenghts, unsigned int count);`
* `int pll_set_tip_states(pll_partition_t * partition, unsigned int tip_index, const unsigned int * map, const char * sequence);`
* `void pll_set_tip_clv(pll_partition_t * partition, unsigned int tip_index, const double * clv);`
* `void pll_set_pattern_weights(pll_partition_t * partition, const unsigned int * pattern_weights);`
* `int pll_update_invariant_sites_proportion(pll_partition_t * partition, unsigned int params_index, double prop_invar);`
* `int pll_update_invariant_sites(pll_partition_t * partition);`

### Likelihood computation

* `void pll_update_partials(pll_partition_t * partition, const pll_operation_t * operations, unsigned int count);`
* `double pll_compute_root_loglikelihood(pll_partition_t * partition, unsigned int clv_index, int scaler_index, unsigned int freqs_index);`
* `double pll_compute_edge_loglikelihood(pll_partition_t * partition, unsigned int parent_clv_index, int parent_scaler_index, unsigned int child_clv_index, int child_scaler_index, unsigned int matrix_index, unsigned int freqs_index);`

### Output functions

* `void pll_show_pmatrix(pll_partition_t * partition, unsigned int index, unsigned int float_precision);`
* `void pll_show_clv(pll_partition_t * partition, unsigned int index, int scaler_index, unsigned int float_precision);`

### Functions for parsing files

* `pll_fasta_t * pll_fasta_open(const char * filename, const unsigned int * map);`
* `int pll_fasta_getnext(pll_fasta_t * fd, char ** head, long * head_len, char ** seq, long * seq_len, long * seqno);`
* `void pll_fasta_close(pll_fasta_t * fd);`
* `long pll_fasta_getfilesize(pll_fasta_t * fd);`
* `long pll_fasta_getfilepos(pll_fasta_t * fd);`
* `int pll_fasta_rewind(pll_fasta_t * fd);`
* `pll_utree_t * pll_utree_parse_newick(const char * filename, unsigned int * tip_count);`
* `pll_rtree_t * pll_rtree_parse_newick(const char * filename, unsigned int * tip_count);`
* `void pll_utree_destroy(pll_utree_t * root);`
* `void pll_rtree_destroy(pll_rtree_t * root);`

### Tree manipulation functions

* `void pll_utree_show_ascii(pll_utree_t * tree);`
* `void pll_rtree_show_ascii(pll_rtree_t * tree);`
* `char * pll_utree_export_newick(pll_utree_t * root);`
* `char * pll_rtree_exprot_newick(pll_rtree_t * root);`
* `int pll_utree_query_tipnodes(pll_utree_t * root, pll_utree_t ** node_list);`
* `int pll_utree_query_innernodes(pll_utree_t * root, pll_utree_t ** node_list);`
* `void pll_utree_create_operations(pll_utree_t ** trav_buffer, unsigned int trav_buffer_size, double * branches, unsigned int * pmatrix_indices, pll_operation_t * ops, unsigned int * matrix_count, unsigned int * ops_count);`
* `int pll_rtree_query_tipnodes(pll_rtree_t * root, pll_rtree_t ** node_list);`
* `int pll_rtree_query_innernodes(pll_rtree_t * root, pll_rtree_t ** node_list);`
* `void pll_rtree_create_operations(pll_utree_t ** trav_buffer, unsigned int trav_buffer_size, double * branches, unsigned int * pmatrix_indices, pll_operation_t * ops, unsigned int * matrix_count, unsigned int * ops_count);`
* `int pll_utree_traverse(pll_utree_t * tree, int (*cbtrav)(pll_utree_t *), pll_utree_t ** outbuffer);`
* `int pll_rtree_traverse(pll_rtree_t * tree, int (*cbtrav)(pll_rtree_t *), pll_rtree_t ** outbuffer);`

### Core functions

* `void pll_core_update_partial(unsigned int states, unsigned int sites, unsigned int rate_cats, double * parent_clv, unsigned int * parent_scaler, const double * left_clv, const double * right_clv, const double * left_matrix, const double * right_matrix, const unsigned int * left_scaler, const unsigned int * right_scaler, unsigned int attrib);`

### Auxiliary functions

* `int pll_compute_gamma_cats(double alpha, unsigned int categories, double * output_rates);`

## Usage examples

Please refer to the [wiki page](https://github.com/xflouris/libpll/wiki) and/or the [examples directory](https://github.com/xflouris/libpll/tree/master/examples).

## Documentation

Please refer to the [wiki page](https://github.com/xflouris/libpll/wiki).

## libpll license and third party licenses

The code is currently licensed under the [GNU Affero General Public License version 3](http://www.gnu.org/licenses/agpl-3.0.en.html).

## Code

The code is written in C.

    File              | Description
----------------------|----------------
**fasta.c**           | Functions for parsing FASTA files.
**gamma.c**           | Functions related to Gamma (&Gamma;) function.
**lex_rtree.l**       | Lexical analyzer parsing newick rooted trees.
**lex_utree.l**       | Lexical analyzer parsing newick unrooted trees.
**likelihood.c**      | Likelihood computation functions.
**core_likelihood.c** | Core likelihood functions that do not require partition instances.
**list.c**            | (Doubly) Linked-list implementations.
**Makefile**          | Makefile.
**maps.c**            | Character mapping arrays for converting sequences to the internal representation.
**models.c**          | Model parameters related functions.
**output.c**          | Functions for output in terminal (i.e. conditional likelihood arrays, probability matrices).
**pll.c**             | Functions for setting PLL partitions (instances).
**rtree.c**           | Rooted tree manipulation functions.
**utree.c**           | Unrooted tree manipulation functions.
**parse_rtree.y**     | Functions for parsing rooted trees in newick format.
**parse_utree.y**     | Functions for parsing unrooted trees in newick format.

## Bugs

The source code in the master branch is thoroughly tested before commits.
However, mistakes may happen. All bug reports are highly appreciated.

## libpll core team

* Tom&aacute;&scaron; Flouri
* Diego Darriba
* Kassian Kobert
* Mark T. Holder
* Alexandros Stamatakis

## Contributing to libpll

Please read the section [Contributing to `libpll`](https://github.com/xflouris/libpll/wiki#contributing-to-libpll) 
of the [wiki](https://github.com/xflouris/libpll/wiki).

## References

* Flouri T., Izquierdo-Carrasco F., Darriba D., Aberer AJ, Nguyen LT, Minh BQ, von Haeseler A., Stamatakis A. (2014)
**The Phylogenetic Likelihood Library.**
*Systematic Biology*, 64(2): 356-362.
doi:[10.1093/sysbio/syu084](http://dx.doi.org/10.1093/sysbio/syu084)

* Gu X., Fu YX, Li WH. (1995)
**Maximum Likelihood Estimation of the Heterogeneity of Substitution Rate among Nucleotide Sites.**
*Molecular Biology and Evolution*, 12(4): 546-557.

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
