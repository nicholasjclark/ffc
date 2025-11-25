# Define functions with dynamic coefficients in ffc formulae

Set up smooth terms with time-varying (dynamic) coefficients for use in
ffc models

## Usage

``` r
fts(
  ...,
  mean_only = FALSE,
  k = NA,
  time_k = 10,
  bs = "cr",
  time_bs = "ts",
  m = NA,
  time_m = 2,
  share_penalty = TRUE,
  d = NA,
  by = NA,
  xt = NULL,
  pc = NULL
)
```

## Arguments

- ...:

  a list of variables that are the covariates that this smooth is a
  function of. Transformations whose form depends on the values of the
  data are best avoided here: e.g. `fts(log(x), z)` is fine, but
  `fts(I(x / sd(x)), z)` is not.

- mean_only:

  `Logical` indicating whether to only include a single basis function
  to be modelled as time-varying. This can be helpful if you wish to
  include a temporal trend that can then be forecasted ahead using
  appropriate time series models. Default is `FALSE`

- k:

  the dimension(s) of the bases used to represent the smooth term. If
  not supplied then set either to `10` (if only a single covariate is
  supplied) or to `5 ^ d`, where `d` is the number of covariates
  supplied. If supplied as a single number then this basis dimension is
  used for each basis. If supplied as an array then the elements are the
  dimensions of the component (marginal) bases of the tensor product.

- time_k:

  the dimension of the bases to be used in the time-varying coefficient
  smooths (see **Details** below). Arbitrarily set to `10` by default.

- bs:

  array (or single character string) specifying the type for each
  marginal basis. `"cr"` for cubic regression spline; `"cs"` for cubic
  regression spline with shrinkage; `"cc"` for periodic/cyclic cubic
  regression spline; `"tp"` for thin plate regression spline; `"ts"` for
  t.p.r.s. with extra shrinkage. See
  [`smooth.terms`](https://rdrr.io/pkg/mgcv/man/smooth.terms.html) for
  details and full list. User defined bases can also be used here (see
  [`smooth.construct`](https://rdrr.io/pkg/mgcv/man/smooth.construct.html)
  for an example). If only one basis code is given then this is used for
  all bases.

- time_bs:

  a two letter `character` string indicating the (penalized) smoothing
  basis to use for the time-varying basis coefficients (eg `"tp"` for
  thin plate regression spline, `"cr"` for cubic regression spline). see
  [`smooth.terms`](https://rdrr.io/pkg/mgcv/man/smooth.terms.html) for
  an over view of what is available. It is generally recommended that
  you stick with one of the doubly-penalized bases (i.e. `"cs"` or
  `"ts"`) as this helps to ensure the resulting basis function
  coefficient time series are estimated on an appropriate scale for
  later forecasting

- m:

  The order of the spline and its penalty (for smooth classes that use
  this) for each term. If a single number is given then it is used for
  all terms. A vector can be used to supply a different `m` for each
  margin. For marginals that take vector `m` (e.g.
  [`p.spline`](https://rdrr.io/pkg/mgcv/man/smooth.construct.ps.smooth.spec.html)
  and
  [`Duchon.spline`](https://rdrr.io/pkg/mgcv/man/smooth.construct.ds.smooth.spec.html)),
  then a list can be supplied, with a vector element for each margin.
  `NA` autoinitializes. `m` is ignored by some bases (e.g. `"cr"`).

- time_m:

  the order of the penalty for the time-varying coefficient smooths
  (e.g. `2` for normal cubic spline penalty with 2nd derivatives). Only
  some smooth classes use this. The `"ps"` class can use a `2` item
  `array` giving the basis and penalty order separately.

- share_penalty:

  `logical` specifying whether the time-varying coefficient smooths for
  this term should all share a smoothing penalty. Defaults to `TRUE` for
  single-parameter families, but automatically set to `FALSE` for
  distributional families (list formulae) to prevent mgcv fitting
  issues. Changing to `FALSE` gives more flexibility to capture
  time-varying functions.

- d:

  array of marginal basis dimensions. For example if you want a smooth
  for 3 covariates made up of a tensor product of a 2 dimensional
  t.p.r.s. basis and a 1-dimensional basis, then set `d=c(2,1)`.
  Incompatibilities between built in basis types and dimension will be
  resolved by resetting the basis type.

- by:

  a `factor` variable of the same dimension as each covariate, used to
  create a replicate of the smooth for each factor level.

- xt:

  Either a single object, providing any extra information to be passed
  to each marginal basis constructor, or a list of such objects, one for
  each marginal basis.

- pc:

  If not `NULL`, signals a point constraint: the smooth should pass
  through zero at the point given here (as a vector or list with names
  corresponding to the smooth names). Never ignored if supplied. See
  [`identifiability`](https://rdrr.io/pkg/mgcv/man/identifiability.html).

## Value

A class `fts.spec` object defining the smooth function to be evaluated
and configured for estimating time-varying basis function coefficients

## Details

ffc will evaluate the basis from these smooths and add the basis
functions to the internal model design matrix. Once these basis function
predictors are added to the data, a model will then be estimated that
allows their coefficients to change through time using terms such as
`s(time, by = bfun_1, id = 1, bs = 'cr') + s(time, by = bfun_2, id = 1, bs = 'cr') + ...`.
By linking the smoothing parameters using the `id` argument, the
time-varying function will be efficiently regularised.

## Author

Nicholas J Clark

## Examples

``` r
# Define a time-varying smooth function of age
fts(age, k = 8, bs = "cr", time_k = 10)
#> $call
#> [1] "s(age, k = 8, bs = \"cr\")"
#> 
#> $term
#> [1] "age"
#> 
#> $by
#> [1] "NA"
#> 
#> $time_bs
#> [1] "ts"
#> 
#> $time_k
#> [1] 10
#> 
#> $time_m
#> [1] 2
#> 
#> $label
#> [1] "fts_age"
#> 
#> $mean_only
#> [1] FALSE
#> 
#> $share_penalty
#> [1] TRUE
#> 
```
