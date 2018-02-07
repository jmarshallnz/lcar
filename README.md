# lcaR

An implementation of inference for prevalence, sensitivity and specificity where we have
multiple (independent) diagnostic tests in the absence of a gold standard across multiple (independent)
populations.

This code can handle arbitrary numbers of tests and populations. Note that for 2 tests and
2 (independent) populations the model is just identifiable (6 independent parameters, 6
independent data), however for other cases (e.g. 2 tests, 1 population) this will no longer
be the case.

## Model description

A description of the statistical model, including likelihood specification, priors, and
conditional posteriors (via a Gibbs MCMC sampler) can be found in the [`lcar.Rmd`(https://github.com/jmarshallnz/lcar/blob/master/lcar.Rmd) vignette.

A pre-built version of this is available in [`mlcar.pdf`](https://github.com/jmarshallnz/lcar/blob/master/lcar.pdf).

## Running

The [`lca.R`](https://github.com/jmarshallnz/lcar/blob/master/lcar.R) script contains all needed R code. Data are specified at the top of the file,
as are priors on parameters, which are particularly important when the model is non-identifiable.

Then the number of iterations from the MCMC Gibbs sampler are specified, and the MCMC chain is
run. Output of the posterior is then plot towards the end of the file.



