# List of examples

## [`rooted`](https://github.com/xflouris/libpll/tree/master/examples/rooted)

This examples shows how to evaluate the log-likelihood of a rooted tree by
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
* considering the Inv+&Gamma; model.
* displaying the CLVs and probability matrices on the console,
* deallocating the partition instance.
