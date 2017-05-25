# Test cases

In this directory go the source of the test cases. For each test case here,
there must be a matching expected output in out/ directory.

When implementing a test case it is necessary to ensure that there is no
errors in the output, and of course it must be deterministic. Afterwards,
place the expected output in the out/ directory and re-run the validations
again just to make sure that everything is working properly.

Each test must focus on evaluating one or a reduced set of features. Also
testing cases that should fail and it is determined how the library would
behave is interesting. For example, attempting to read an unexistent file.

* 00010_NMPU_lkcalc compute likelihood score for a simple unrooted tree
* 00020_NMPR_lkcalc compute likelihood score for a simple rooted tree

## alpha-cats

Evaluate the likelihood for different alpha shape parameters and number of
categories.

## blopt-minimal

(optimize module) Optimize branch lengths for a minimal tree with 3 tips and
3 branches.

## derivatives

Evaluate the computation of the likelihood derivatives at different branch 
lengths on a small tree and msa.

The derivatives are computed twice at an inner edge and at a tip edge using 3 
different alphas, 4 proportion of invariant sites, 3 sets of rate categories 
and 9 branches ranging from 0.1 to 90.

## derivatives-oddstates

Analogous to `derivatives` but using an odd number of states.

## fasta-dna

Read a DNA MSA in FASTA format, load the sequences into the PLL partition 
and evaluate the likelihood.

## fasta-prot

Read an amino acid MSA in FASTA format. Checks also that 'pll_fasta_get_next' 
would fail if the states map is wrong (dna instead of protein map). It does
not evaluate the likelihood.

## hky

Evaluate the likelihood for different transition-transversion ratios in
HKY models.

## odd-states

Evaluate the likelihood for a data set with 7 states. This is specially
important where vector intrinsics are used and the states are padded to fit
the alignment.

## partial-traversal

Perform partial traversals on the tree.

## protein-models

Evaluate the likelihood of a short sequence under all the available empirical 
amino acid replacement models

## rooted

Evaluate the likelihood of a tree rooted at an inner-inner node with 
4 categories and 4 different proportions of invariant sites, from 0.0 to 0.9

    /\
   /  \
  /    \
 /\    /\
/  \  /  \

## rooted-tipinner

Evaluate the likelihood of a tree rooted at a tip-inner node with 4 categories 
and 4 different proportions of invariant sites, from 0.0 to 0.9

    /\
   /  \
  /    \
 /     /\
/     /  \

## treemove-nni

Validate Nearest Neighbor Interchange moves.

## scaling

Validate CLV scaling on large trees (per-site and per-rate scaling modes)

## pmatrix

Validate pmatrix computation, and specifically check for negative values which could
appear due to numerical issues

NOTE: expected output file for this test is intentionally missing, since negative 
values problem has not been fully fixed yet 

## treemove-spr

Validate Subtree Prunning and Regrafting moves.

## treemove-tbr

Perform bisection, reconnection and local branch length optimization.
