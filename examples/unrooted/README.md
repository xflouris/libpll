# Unrooted tree example

This examples evaluates the log-likelihood of the unrooted tree presented in the
figure below, by creating a custom post-order traversal that drives the
likelihood computation.

![unrooted tree](https://github.com/xflouris/assets/raw/master/libpll/images/unrooted.png)

## Explanation of the example

The program first instantiates a partition using the function call

```C
partition = pll_create_partition(4,
                                 2,
                                 4,
                                 6,
                                 1,
                                 5,
                                 4,
                                 2,
                                 PLL_ATTRIB_ARCH_AVX);
```

The parameters of the function (in the order passed) indicate
* the number of tip sequences that will be used for this partition (4),
* the extra number of Conditional Likelihood Vectors (CLVs) that should be allocated apart from those created for the tip sequences (typically, number of tips minus two for unrooted trees) (2),
* number of states in the dataset (for instance 4 for nucleotide datasets, 20 for aminoacid datasets) (4),
* the length of the alignment, i.e. the number of sites at the tip sequences (6),
* how many different substitution models we want to have at one time (1),
* the number of probability matrices that should be allocated (typically 2 times the number of tip sequences minus 3 for unrooted trees) (5) ,
* number of discrete rate categories (rate heterogeneity) (4),
* number of scale buffers to be allocated (2),
* attributes that specify what kind of hardware acceleration should be used (PLL\_ATTRIB\_ARCH\_AVX).

For a more detailed explanation of the function arguments refer to the [API Reference](https://github.com/xflouris/libpll/wiki/API-Reference#pll_create_partition).

Model parameters are set with the function calls

[`pll_set_frequencies(partition, 0, frequencies);`](https://github.com/xflouris/libpll/wiki/API-Reference#void-pll_set_frequencies)

[`pll_set_subst_params(partition, 0, subst_params, 6);`](https://github.com/xflouris/libpll/wiki/API-Reference#void-pll_set_subst_params)

`pll_set_category_rates(partition, rate_cats);`


The CLVs at tips are set by calling, for example,

[`pll_set_tip_states(partition, 0, pll_map_nt, "WAAAAB");`](https://github.com/xflouris/libpll/wiki/API-Reference#void-pll_set_tip_states)

This function sets the sequence for tip 0 (first tip). Note that it is
necessary to specify a map for converting the sequence to partials. The map is
essentially a 256-long array of  `unsigned int` values that map each ASCII
character to a number. The set bits of the number indicate which CLV entries
will be set to 1.  For nucleotide data, the predefined map `pll_map_nt` may be
used, which accounts for degenerate (ambiguous) base-calls as well. Also note
that CLVs for tips are numbered from 0 to _n-1_, where _n_ is the number of tips.
For inner node CLVs you may use any number starting from _n_ upto _n+x-1_, where
_x_ is the second parameter of [`pll_create_partition(...)`](https://github.com/xflouris/libpll/wiki/API-Reference#pll_create_partition).

Next, we compute the transition probability matrices for each unique branch
length in the tree by calling the function

```C
unsigned int params_indices[4] = {0,0,0,0};

pll_update_prob_matrices(partition,
                         params_indices,
                         matrix_indices,
                         branch_lengths,
                         5);
```
which computes the probability matrices at indices `matrix_indices` from the
corresponding branch lengths specified in `branch_lengths` and the
corresponding rate matrices whose index is specified with `params_indices`
(second parameter). The last argument indicates the size of the last two
arrays. Note that the function will compute probability matrices for all
available rate categories. For more information on this function check the
[documentation](https://github.com/xflouris/libpll/wiki/Updating-transition-probability-matrices).

The next step is to create a traversal descriptor for driving the likelihood
computation. This is done by allocating a `pll_operation_t` structure which we
fill with information on how to compute each CLV corresponding to an internal
node.

```C
operations[0].parent_clv_index    = 4;
operations[0].child1_clv_index    = 0;
operations[0].child2_clv_index    = 1;
operations[0].child1_matrix_index = 0;
operations[0].child2_matrix_index = 1;
operations[0].parent_scaler_index = 0;
operations[0].child1_scaler_index = PLL_SCALE_BUFFER_NONE;
operations[0].child2_scaler_index = PLL_SCALE_BUFFER_NONE;
```

`parent_clv_index` indicates which CLV we want to update/compute, using the
CLVs specified by `child1_clv_index` and `child2_clv_index`. The probability
matrices used for each child CLV is specified by `child1_matrix_index` and
`child2_matrix_index`, respectively.  Note that the node indices in the tree
illustration directly correspond to the indices used in the example program.
For the two tip nodes (child1 and child2) we use no scalers for maintaining
numerical stability, and we use one scaler for the inner node with index 4.
Note that, a scaler is an array of integer values which holds how many times
each partial entry was scaled.

Now we can use the created `pll_operation_t` structure to compute the CLVs by
calling

`pll_update_partials(partition, operations, 2);`

where the last parameter (2) sets the number of operations to be carried in
sequential order (0 to 1).

Finally, we obtain the log-likelihood of the dataset by a call to

```C
logl = pll_compute_edge_loglikelihood(partition,
                                      4,
                                      0,
                                      5,
                                      1,
                                      4,
                                      params_indices,
                                      NULL);
```

where the parameters (in order of appearance) specify the

* partition
* index of one of the end-point nodes of the edge on which the log-likelihood will be computed
* index of the scaler for the particular CLV
* index of the other end-point node
* index of the scaler for the other end-point node CLV
* index of the probability matrix for that particular edge,
* array of indices indicating the frequency array to use for each rate category.
* array to store per-site log-likelihood values, or `NULL` otherwise.

Before exiting, we dispose of the allocated memory by calling

`pll_destroy_partition(partition);`

## Instructions to compile

Before proceeding with the compilation, make sure that `pll.h` is accessible,
by copying it to the current example directory, or modifying the line

`#include "pll.h"`

You will also need to make the shared object `libpll.so` accessible by either
placing it in your system's library directory (typically `/lib` or
`/usr/local/lib`), or by copying it to the current example directory.

Now you may compile the example by running the included Makefile using

`make`

## Instructions to run

Make absolutely sure that the shared object `libpll.so` is accessible by either
placing it in your system's library directory (typically `/lib` or
`/usr/local/lib`). If it is not, you may copy it to a directory (for example
$PWD), and run the command

`export LD_LIBRARY_PATH=$PWD`

Now, run the example by changing to the directory of the compiled file and
typing:

`./unrooted`
