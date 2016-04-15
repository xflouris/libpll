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

## alpha-cats

Evaluate the likelihood for different alpha shape parameters and number of
categories.

## blopt-minimal

(optimize module) Optimize branch lengths for a minimal tree with 3 tips and
3 branches.

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

## treemove-nni

Validate Nearest Neighbor Interchange moves.

## treemove-spr

Validate Subtree Prunning and Regrafting moves.

## treemove-tbr

Perform bisection, reconnection and local branch length optimization.