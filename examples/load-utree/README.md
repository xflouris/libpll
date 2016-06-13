# Rooted tree example

This example demonstrates how to write a function that loads an unrooted binary
tree from a newick file. In case the newick file contains a binary rooted tree,
it is automatically unrooted by removing the root node and summing up the
root's two edges.

## Explanation of the example

The example program composes of the following function

```C
pll_utree_t * load_tree_unrooted(const char * filename,
                                 unsigned int * tip_count)
{
  pll_utree_t * utree;
  pll_rtree_t * rtree;

  if (!(rtree = pll_rtree_parse_newick(filename, tip_count)))
  {
    if (!(utree = pll_utree_parse_newick(filename, tip_count)))
    {
      fprintf(stderr, "%s\n", pll_errmsg);
      return NULL;
    }
  }
  else
  {
    utree = pll_rtree_unroot(rtree);

    /* optional step if using default PLL clv/pmatrix index assignments */
    pll_utree_reset_template_indices(utree, *tip_count);
  }

  return utree;
}
```

which accepts two parameters, the name of the newick file and a variable where
the number of tips will be stored. The function first attempts to load the tree
as a rooted one, and if it fails, it attempts to read it as unrooted.

In case the file contained a rooted binary tree, then it is converted to unrooted
using the call

```C
utree = pll_rtree_unroot(rtree);
```

and, optionally, calls the function

```C
pll_utree_reset_template_indices(utree, *tip_count);
```

to reset CLV, p-matrix and scaler indices to the defaults used by the PLL functions.
That means, tips will be assigned CLV and p-matrix indices ranging from 0 to `tip_count - 1`, and
inner nodes indices from `tip_count` to `2*tip_count - 3`. Tips are assigned `PLL_SCALE_BUFFER_NONE`
for scaler indices (since no scaling is required for tips) and inner nodes get scaler indices from
0 to `tip_count - 3`.

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

`./rooted`
