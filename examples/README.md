# List of examples

## [`rooted`](https://github.com/xflouris/libpll/tree/master/examples/rooted)

This example shows how to evaluate the log-likelihood of a rooted tree by
creating a custom post-order traversal that drives the likelihood computation.
It also demonstrates how to account for invariant sites when computing the
likelihood.

The example covers the following:

* partition instance creation,
* setting variable rates among sites (four discrete rates from &Gamma;),
* setting the parameters for the GTR model,
* computation of probability matrices,
* setting frequencies,
* setting tip CLVs by passing the sequence and a map,
* creating a tree traversal for driving the likelihood function,
* computing the CLVs at inner nodes,
* evaluating the log-likelihood at the root,
* considering the Inv+&Gamma; model,
* displaying the CLVs and probability matrices on the console,
* deallocating the partition instance.

## [`newick-fasta-unrooted`](https://github.com/xflouris/libpll/tree/master/examples/newick-fasta-unrooted)

This example shows how to evaluate the log-likelihood of an
[unrooted binary tree](http://en.wikipedia.org/wiki/Unrooted_binary_tree) tree
which is loadded from a [newick](http://en.wikipedia.org/wiki/Newick_format)
file using the `libpll` newick parser. The alignment is loaded from a FASTA
file using the `libpll` fasta parser, and each sequence is associated to the
corresponding tree taxon via the use of the standard
[GNU C Library (glibc)](http://www.gnu.org/software/libc/) hash table. The
parsed newick tree is stored in an intermediate structure which may be used to
automagically and conveniently create the `libpll` operations structure.

The example covers the following:

* loading a newick tree file using the `libpll` parser into an intermediate unrooted binary tree structure,
* loading sequences from a FASTA file,
* querying the intermediate tree structure for a list of tip names
* associating sequences to taxa using the `glibc` hash table functions,
* getting a filled operations structure by using library functions to iterate over the parsed tree structure,
* setting default branch lengths to branches which did not have an associated length in the newick format,

## [`unrooted`](https://github.com/xflouris/libpll/tree/master/examples/unrooted)

This example shows how to evaluate the log-likelihood at an edge of an
[unrooted binary tree](http://en.wikipedia.org/wiki/Unrooted_binary_tree) tree.

This examples covers the following:

* obtaining the rate categories from a discretized gamma distribution by providing the alpha shape parameter,
* evaluating the log-likelihood at an edge of the unroted tree.

## [`newton`](https://github.com/xflouris/libpll/tree/master/examples/newton)

This example shows how to evaluate the log-likelihood at an edge of an
[unrooted binary tree](http://en.wikipedia.org/wiki/Unrooted_binary_tree) tree,
and how to optimize the length of a given branch using [Newton's method](https://en.wikipedia.org/wiki/Newton's_method).

This examples covers the following:

* obtaining the rate categories from a discretized gamma distribution by providing the alpha shape parameter,
* evaluating the log-likelihood at an edge of the unroted tree.
* optimizing the length of a given branch using Newton's method.
