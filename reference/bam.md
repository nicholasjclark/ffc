# Large dataset GAM fitting

Fits GAMs for large datasets using mgcv::bam. In ffc, specify
`engine = "bam"` in
[`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md)
for large datasets.

## Usage

``` r
bam(
  formula,
  family = gaussian(),
  data = list(),
  weights = NULL,
  subset = NULL,
  na.action = na.omit,
  offset = NULL,
  method = "fREML",
  control = list(),
  select = FALSE,
  scale = 0,
  gamma = 1,
  knots = NULL,
  sp = NULL,
  min.sp = NULL,
  paraPen = NULL,
  chunk.size = 10000,
  rho = 0,
  AR.start = NULL,
  discrete = FALSE,
  cluster = NULL,
  nthreads = 1,
  gc.level = 0,
  use.chol = FALSE,
  samfrac = 1,
  coef = NULL,
  drop.unused.levels = TRUE,
  G = NULL,
  fit = TRUE,
  drop.intercept = NULL,
  in.out = NULL,
  ...
)
```

## See also

[`bam`](https://rdrr.io/pkg/mgcv/man/bam.html),
[`ffc_gam()`](https://nicholasjclark.github.io/ffc/reference/ffc_gam.md)

## Examples

``` r
# Use bam engine for large datasets
# mod <- ffc_gam(y ~ fts(x), data = large_data, time = "time", engine = "bam")
```
