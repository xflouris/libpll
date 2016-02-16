# Model optimization example

This examples optimizes the model parameters and branch lengths for an alignment 
and tree read from command file.

The syntax is as follows:

`model-optimize/model-optimize [newick] [fasta] [model]`

[newick] is a file containing the tree in NEWICK format
[fasta]  is a file containing the MSA in FASTA format
[model]  is the DNA substitution scheme. This a six-integer code that defines 
         how the substitution rate parameters are linked to each other. 
         For example, JC/F81 model is 000000, HKY/K80 = 010010 and 
         GTR/SYM = 012345. 
         There is no need to start in 0. For example, 174673 would work exactly as 010010.
         
## Explanation of the example

The program reads a MSA, a tree and DNA substitution rates scheme, 
and optimizes parameter values and branch lengths.

Following constants define which parameters are optimized:

```C
/* parameters to optimize */
#define OPTIMIZE_BRANCHES     1
#define OPTIMIZE_SUBST_PARAMS 1
#define OPTIMIZE_ALPHA        1
#define OPTIMIZE_FREQS        1
#define OPTIMIZE_PINV         1
```

The thoroughness of the optimization is defined with 2 constants:

```C
/* tolerances */
#define OPT_EPSILON       1e-1
#define OPT_PARAM_EPSILON 1e-2
```

An iteration finishes after optimizing all the selected parameters.
OPT_EPSILON defines the tolerance for the global optimization. The optimization
loop finishes when the difference in the log-Likelihood between starting and 
finishing the iteration is lower than OPT_EPSILON.
OPT_PARAM_EPSILON defines the tolerance for the optimization of each single
parameter. 

Initially, the parameters are set to empirical or default values:

* `pll_compute_empirical_frequencies (partition);` for the frequencies.
* `pll_compute_empirical_subst_rates (partition);` for the substitution rates.
* `1.0` for Gamma shape parameter (alpha)
* `0.5 * pll_compute_empirical_invariant_sites (partition);` for the invariant sites

Model parameters can be optimized using high-level or low-level functions
provided by the library. Low level functions are more flexible, and they are
usually used for achieving the best performance.

After creating the partition as described in previous examples, we set the
parameters and compute the initial logLikelihood score: `pll_update_prob_matrices`,
`pll_update_partials` and `pll_compute_edge_loglikelihood`.

Prior to starting the optimization loop, we need to set up the structures that
we need for calling the optimization functions. We use here 2 approaches. 

The first one uses the high level optimization functions provided by `pll_optimize`.
These functions have the prefix `pll_optimize_parameters_*`, with the suffix
`onedim`, for one single dimention, or `multidim`, for one or more dimensions.
`onedim` uses Brent algorithm and `multidim` L-BFGS-B. Nevertheless, the user
just need to set up the `pll_optimize_options_t` structure. The parameters to
optimize are set in the field `which_parameter` in a mask format. So, for example
if we optimize the Gamma shape parameter only, we set `which_parameter = PLL_OPT_PARAMETER_ALPHA`. 
If we want to optimize the Gamma shape and the invariant sites together, we would
call `multidim` with `which_parameter = PLL_OPT_PARAMETER_ALPHA | PLL_OPT_PARAMETER_PINV`
 
The second one uses the low level optimization functions also provided by `pll_optimize`.
These functions have the prefix `pll_minimize_`, and the suffix is the algorithm to use.
For example, for optimizing several variables together we would use `pll_minimize_lbfgsb`.
Low level optimization functions require a callback *score function* that computes the score 
for a given state. The functions will minimize the score function. 
For the low level optimization, we use a custom score function, `target_*_opt`, with
a custom structure as argument, that includes all the fields we need.

### Optimization loop

The parameters are optimized in the following order:

1.- Branch lengths (Newton-Raphson)
2.- Substitution rates (L-BFGS-B) 
3.- Gamma shape (Brent) 
4.- Proportion of invariant sites (Brent)
5.- State frequencies (L-BFGS-B)

A whole iteration involves the optimization of each parameter. At the end, we
compare the final and initial logLikelihood scores, and we finalize the optimization
when the difference is lower than the optimization threshold (i.e., the process
converged to optimal values).

## Instructions to compile

Before proceeding with the compilation, make sure that `pll_optimize.h`
and `pll.h` are accessible, by copying it to the current example directory, 
or modifying the line

`#include "pll_optimize.h"`

You will also need to make the shared object `libpll_optimize.so` accessible 
by either placing it in your system's library directory (typically `/lib` or
`/usr/local/lib`), or by copying it to the current example directory.

Now you may compile the example by running the included Makefile using

`make`

## Instructions to run

Make absolutely sure that the shared object `libpll_optimize.so` and 
`libpll.so`are accessible by either placing it in your system's library 
directory (typically `/lib` or `/usr/local/lib`). If it is not, you may 
copy it to a directory (for example `$PWD`), and run the command

`export LD_LIBRARY_PATH=$PWD`

Now, run the example by changing to the directory of the compiled file and
typing:

`./model-optimize [tree_file] [dna_msa_file] [model] `

Where `model` is the 6-character representation of the substitution rate
matrix symmetries. This a 6-integer code that defines how the substitution rate
parameters are linked to each other. For example, JC/F81 model is 000000,
HKY/K80 = 010010 and GTR/SYM = 012345. There is no need to start in 0.
For example, 174673 would work exactly as 010010.
