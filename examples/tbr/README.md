# TBR moves example

This examples applies a TBR move to an unrooted topology, computes the marginal
likelihoods and locally optimizes the branch lengths on the reconnection point.

![TBR move](https://github.com/xflouris/assets/raw/master/libpll/images/tbr-move.png)

## Explanation of the example

After reading the tree, we select a random bisection point among the inner nodes.
The bisection point must *not* be a tip branch:

```C
bisect_edge = innernodes[rand () % inner_nodes_count];
```

For the reconnection points, we select two edges among those located at a certain
distance from the bisection point. 
The reconnection point is a structure of type `pll_edge_t`, and it contains 2
nodes, `parent` and `child`, that actually represent branches.

```C
typedef struct
{
  union
  {
    struct
    {
      pll_utree_t * parent;
      pll_utree_t * child;
    } utree;
    struct
    {
      pll_rtree_t * parent;
      pll_rtree_t * child;
    } rtree;
  } edge;
  int additional_pmatrix_index;
  double length;
} pll_edge_t;
```

The following function outputs an array of nodes, `nodes_at_dist`, and its size,
`n_nodes_at_dist`, that are located at `distance` from the bisection point, 
`bisect_edge`.

```C
/* from one side of the bisection */
pll_utree_nodes_at_node_dist (bisect_edge, 
                              nodes_at_dist, 
                              &n_nodes_at_dist,
                              distance, 
                              fixed);
/* select random edge */
reconnect.edge.utree.parent = nodes_at_dist[rand () % n_nodes_at_dist];

/* and from the other side */
pll_utree_nodes_at_node_dist (bisect_edge->back, 
                              nodes_at_dist,
                              &n_nodes_at_dist, 
                              distance, 
                              fixed);
/* select random edge */
reconnect.edge.utree.child = nodes_at_dist[rand () % n_nodes_at_dist];
                                  
``` 

The parameter `fixed` (boolean) means wheter the selected nodes are located at
exactly the `distance` provided. If `fixed` is false, then `distance` represents
the search radius (i.e., will select *all* nodes within a distance of `distance`
or less.

We have set the `parent` reconnection edge first, and then the `child`. 
There is no difference if you do it the other way around.

After choosing the bisection and reconnection points, we apply the move:

```C
pll_utree_TBR (bisect_edge, &reconnect)
```

We check the integrity of the new tree:

```C
pll_utree_check_integrity (tree)
```

With the new tree topology, we need to recompute the tree traversal and update
the conditional likelihood vectors. Afterwards, we can compute the marginal 
likelihoods at the reconnection points:

```C
  /* compute marginal likelihoods */
  printf("\nMarginal likelihoods:\n");
  logl = pll_compute_root_loglikelihood (partition,
                                         tree->clv_index,
                                         tree->scaler_index,
                                         frequencies_index);
  printf ("  Log-L Partial at %s: %f\n", tree->label, logl);

  logl = pll_compute_root_loglikelihood (partition,
                                         tree->back->clv_index,
                                         tree->back->scaler_index,
                                         frequencies_index);
  printf ("  Log-L Partial at %s: %f\n", tree->back->label, logl);
```

And the new likelihood of the tree:

```C
  /* compute global likelihood */
  logl = pll_compute_edge_loglikelihood (partition,
                                         tree->clv_index,
                                         tree->scaler_index,
                                         tree->back->clv_index,
                                         tree->back->scaler_index,
                                         tree->pmatrix_index,
                                         frequencies_index);
```

Finally, we perform a local branch length optimization (4 `smoothings`) in the 
new branches around the reconnection points within a `radius` of 2 nodes. 
The last parameter, `keep_update` (boolean), determines whether the optimal
length of each branch is updated before evaluating the next branch.

```C
  pll_optimize_branch_lengths_local (partition,
                                     tree,
                                     parameters_index,
                                     frequencies_index,
                                     1e-2, /* tolerance */
                                     4,    /* smoothings */
                                     2,    /* radius */
                                     1);   /* keep update */
```

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

`./tbr [tree_file] [dna_msa_file] [model]`
