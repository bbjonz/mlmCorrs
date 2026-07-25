# mlmCorrs
R package for computing multilevel correlation matrices and intraclass correlations.

Requires a dataframe with numeric/integer columns and, for the multilevel functions, the cluster variable must be included in the dataframe.

- `corstars()` — APA-style correlation table with significance stars, optionally split by a grouping variable
- `icc.corrs()` — correlation table combined with intraclass correlations (ICC1/ICC2) from a random-intercept model per variable
- `lgm()` — two-level (within/between) correlation matrix fit via a saturated multilevel SEM

# Installation
```r
remotes::install_github("bbjonz/mlmCorrs")
```
