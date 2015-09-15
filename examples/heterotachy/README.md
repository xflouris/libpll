# Heterotachous models example

This examples evaluates the log-likelihood of the unrooted tree presented in the
figure below, by creating a custom post-order traversal that drives the
likelihood computation. The model used is accounts for heterotachy, using 3
sets of parameters applied on different branches.

![unrooted tree](https://github.com/xflouris/assets/raw/master/libpll/images/unrooted.png)

## Explanation of the example

The example is based on the unrooted tree example. Here are explained the
features related to heterotachous models. 

```C
partition = pll_create_partition(4, 
                                 2, 
                                 4, 
                                 6,
                                 0, 
                                 num_models, 
                                 5, 
                                 4, 
                                 1, 
                                 PLL_ATTRIB_ARCH_SSE);
```

The number of different model parameters (substitution rates and frequencies)
has to be defined when calling the pll_create_partition. It might happen that
these number changes through time, for example if we want to optimize the number
of per-branch models. In that case, here specify the **maximum** number of
different models we would like to store. 

For a more detailed explanation of the function arguments refer to the [API Reference](https://github.com/xflouris/libpll/wiki/API-Reference#pll_create_partition).

Model parameters are set for each of the different models. The frequencies and 
substitution parameters arrays are now two-dimensional, with `num_models` sets
of values.

[`pll_set_frequencies(partition, i, 0, frequencies[i]);`](https://github.com/xflouris/libpll/wiki/API-Reference#void-pll_set_frequencies)

[`pll_set_subst_params(partition, i, 0, subst_params[i], 6);`](https://github.com/xflouris/libpll/wiki/API-Reference#void-pll_set_subst_params)

Next, we compute the transition probability matrices for each unique branch
length in the tree. The matrix indices are divided in *num_models* groups, one
for each different model. The way to define which model correspond to which
branch is up to the user. Here we propose a simple method where matrices
between `matrix_start[i]` and `matrix_start[i]+matrix_count[i]` are computed
with the *i-th* model parameters. 

```C
  unsigned int matrix_indices[5] = {0,1,2,3,4};
  unsigned int matrix_start[3]   = {0,2,4};
  unsigned int matrix_count[3]   = {2,2,1};
```

And the assignment is done in the custom function `update_pmatrix_het`
 
```C
static void update_pmatrix_het(pll_partition_t *partition,
                               unsigned int * matrix_indices,
                               unsigned int * matrix_start,
                               unsigned int * matrix_count,
                               double * branch_lengths)
{
  unsigned int i;
  for (i=0; i<num_models; i++)
      pll_update_prob_matrices(partition,
                           i,
                           &(matrix_indices[matrix_start[i]]),
                           &(branch_lengths[matrix_start[i]]),
                           matrix_count[i]);
}
```

The traversal descriptor, the CLVs are computed exactly as in the single 
model case. The only difference would be when computing the edge log-likelihood,
where the corresponding index of the frequencies in the root branch must be
specified.

```C
logl = pll_compute_edge_loglikelihood(partition,
                                      4,
                                      5,
                                      4,
                                      frequencies_index);
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
