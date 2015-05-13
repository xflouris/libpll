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

## fasta-parser

Read a MSA in FASTA format, load the sequences into the PLL partition and
evaluate the likelihood.

## hky

Evaluate the likelihood for different transition-transversion ratios in
HKY models.

## protein-models

Evaluate the likelihood of a short sequence under all the available empirical 
amino acid replacement models
