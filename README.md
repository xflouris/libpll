# libpll

## Introduction

The aim of this project is to implement a versatile high-performance software library
for phylogenetic analysis. The library should serve as a lower-level interface of the
PLL (Flouri et al. 2014) and should have the following properties:

* have an open source code with an appropriate open source license.
* 64-bit multi-threaded design that handles very large datasets.
* easy to use and well-documented.
* SIMD implementations of time-consuming parts.
* provide comparably fast likelihood computations as RAxML (Stamatakis 2014)

## Compilation instructions

Currently, libpll can be compiled using the included Makefile:

`make`

## Available functionality

libpll implements the General Time Reversible (GTR) model (Tavare 1986) which
can be used for nucleotide and amino acid data. It also supports models of
variable rates among sites, and has functions for computing the discretized
rate categories for the gamma model (Yang 1994).

Below is list of available functions in the current version.

### Partition (instance) manipulation

* `pll_partition_t * pll_create_partition(int tips, int clv_buffers, int states, int sites, int rate_matrices, int prob_matrices, int rate_cats, int scale_buffers, int attributes);`
* `int pll_destroy_partition(pll_partition_t * partition);`

### Linked lists

* `int pll_dlist_append(pll_dlist_t ** dlist, void * data);`
* `int pll_dlist_prepend(pll_dlist_t ** dlist, void * data);`
* `int pll_dlist_remove(pll_dlist_t ** dlist, void * data);

### Models setup

* `void pll_set_subst_params(pll_partition_t * partition, int params_index, double * params, int count);`
* `void pll_set_frequencies(pll_partition_t * partition, pll_partition_t * partition, int params_index, double * frequencies);`
* `void pll_set_category_rates(pll_partition_t * partition, double * rates);`
* `void pll_update_prob_matrices(pll_partition_t * partition, int params_index, int * matrix_indices, double * branch_lenghts);`
* `void pll_set_tip_states(pll_partition_t * partition, int tip_index, const char * sequence);`

### Likelihood computation

* `void pll_update_partials(pll_partition_t * partition, pll_operation_t * operations, int count);`
* `double pll_compute_root_loglikelihood(pll_partition_t * partition, int clv_index, int freqs_index);`
* `double pll_compute_edge_loglikelihood(pll_partition_t * partition, int parent_clv_index, int child_clv_index, int matrix_index, int freqs_index);`

### Output functions

* `void pll_show_pmatrix(pll_partition_t * partition, int index);`
* `void pll_show_clv(pll_partition_t * partition, int index);`

## Usage examples

Please refer to our wiki page.

## libpll license and third party licenses

The code is currently licensed under the GNU Affero General Public License version 3.

## Code

The code is written in C.

    File         | Description
-----------------|----------------
**likelihood.c** | Likelihood computation functions.
**list.c**       | (Doubly) Linked-list implementations.
**Makefile**     | Makefile
**maps.c**       | Character mapping arrays for converting sequences to the internal representation.
**models.c**     | Model parameters related functions
**output.c**     | Functions for output in terminal (i.e. conditional likelihood arrays, probability matrices)
**pll.c**        | Functions for setting PLL partitions (instances).

## Bugs

The source code in the master branch is thoroughly tested before commits.
However, mistakes may happen. All bug reports are highly appreciated.

## libpll core team

* Tom&aacute;&scaron; Flouri
* Diego Darriba
* Alexandros Stamatakis

## References

* Flouri T., Izquierdo-Carrasco F., Darriba D., Aberer AJ, Nguyen LT, Minh BQ, von Haeseler A., Stamatakis A. (2014)
**The Phylogenetic Likelihood Library.**
*Systematic Biology*, 64(2): 356-362.
doi:[10.1093/sysbio/syu084](http://dx.doi.org/10.1093/sysbio/syu084)

* Stamatakis A. (2014)
**RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenie.s**
*Bioinformatics*, 30(9): 1312-1313.
doi:[10.1093/bioinformatics/btu033](http://dx.doi.org/10.1093/bioinformatics/btu033)

* Tavar&eacute;  S. (1986)
**Some probabilistic and statistical problems in the analysis of DNA sequences**
*American Mathematical Sciety: Lectures on Mathematics in the Life Sciences*, 17: 57-86.

* Yang Z. (2014)
**Maximum likelihood phylogenetic estimation from dna sequences with variable rates over sites: Approximate methods.**
*Journal of Molecular Evolution*, 39(3): 306-314.
doi:[10.1007/BF00160154](http://dx.doi.org/10.1007/BF00160154)
