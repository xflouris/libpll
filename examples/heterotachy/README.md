# Heterotachous models example

This examples evaluates the log-likelihood of the unrooted tree presented in
the figure below, by creating a custom post-order traversal that drives the
likelihood computation. The model used accounts for heterotachy, using 3
different rate matrices, which are applied on different branches of the tree.

![unrooted tree](https://github.com/xflouris/assets/raw/master/libpll/images/unrooted.png)

## Explanation of the example

The example is based on the unrooted tree example. Here are explained the
features related to heterotachous models.

```C
partition = pll_partition_create(4,
                                 2,
                                 4,
                                 6,
                                 rmatrix_count,
                                 5,
                                 4,
                                 1,
                                 PLL_ATTRIB_ARCH_SSE);
```

The number of different rate matrices (composed from substitution rates and
frequencies) has to be defined when calling the function
`pll_partition_create`. It might happen that these number changes through time,
for example if we want to optimize the number of per-branch models. In that
case, here specify the **maximum** number of different models we would like to
store.

For a more detailed explanation of the function arguments refer to the [API Reference](https://github.com/xflouris/libpll/wiki/API-Reference#pll_partition_create).

We set the frequencies and substitution rates for each of the `rmatrix_count`
rate matrices.

[`pll_set_frequencies(partition, i, 0, frequencies[i]);`](https://github.com/xflouris/libpll/wiki/API-Reference#void-pll_set_frequencies)

[`pll_set_subst_params(partition, i, 0, subst_params[i], 6);`](https://github.com/xflouris/libpll/wiki/API-Reference#void-pll_set_subst_params)

Next, we compute the transition probability matrices for each unique branch
length in the tree. The plan here is to create the first two transition
matrices (p-matrices) using the first rate matrix, p-matrices 3 and 4 using the
second rate matrix, and p-matrix 5 using the third rate matrix. We accomplish
this assignment with the following three arrays

```C
  unsigned int matrix_indices[5] = {0,1,2,3,4};
  unsigned int matrix_start[3]   = {0,2,4};
  unsigned int matrix_count[3]   = {2,2,1};
```

and the custom function `update_pmatrices`

```C
static void update_pmatrix_het(pll_partition_t * partition,
                               unsigned int * matrix_indices,
                               unsigned int * matrix_start,
                               unsigned int * matrix_count,
                               double * branch_lengths)
{
  unsigned int i;
  unsigned int params_indices[4];

  for (i=0; i<num_models; i++)
    for (j=0; j < 4; ++j)
      params_indices[j] = i;
    pll_update_prob_matrices(partition,
                             params_indices,
                             matrix_indices + matrix_start[i],
                             branch_lengths + matrix_start[i],
                             matrix_count[i]);
}
```
Note that, each p-matrix in fact comprises of four p-matrices, one for each
rate category. Hence, the second parameter to `pll_update_prob_matrices`
specifies which rate matrix should be used for each rate category. Since we
want to use the same rate matrix, we duplicate the index in the array
`params_indices`.


The operation structures and the computation of the conditional probability
vectors is done exactly as in the single rate matrix case. The only difference
is in computing the log-likelihood, where the corresponding frequencies must be
specified. We use an auxiliary array `freqs_indices` filled with 2's, to
indicate that the frequencies array with index 2 should be used for all 4
p-matrix (remember we have four rates) when evaluating the log-likelihood.

```C
unsigned int freqs_indices[4] = {2,2,2,2};

logl = pll_compute_edge_loglikelihood(partition,
                                      4,
                                      5,
                                      4,
                                      freqs_indices,
                                      NULL);
```

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

`./heterotachy`
