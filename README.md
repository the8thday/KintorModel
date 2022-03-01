
<!-- README.md is generated from README.Rmd. Please edit that file -->

# KintorModel

<!-- badges: start -->
<!-- badges: end -->

The goal of KintorModel is a backup of scripts for model build.

it uses tidymodels as the main package, and introduce detail steps in
clinical research.

``` r
library(tidymodels)
#> Registered S3 method overwritten by 'tune':
#>   method                   from   
#>   required_pkgs.model_spec parsnip
#> ── Attaching packages ────────────────────────────────────── tidymodels 0.1.4 ──
#> ✓ broom        0.7.9      ✓ recipes      0.1.17
#> ✓ dials        0.0.10     ✓ rsample      0.1.0 
#> ✓ dplyr        1.0.7      ✓ tibble       3.1.5 
#> ✓ ggplot2      3.3.5      ✓ tidyr        1.1.4 
#> ✓ infer        1.0.0      ✓ tune         0.1.6 
#> ✓ modeldata    0.1.1      ✓ workflows    0.2.4 
#> ✓ parsnip      0.1.7      ✓ workflowsets 0.1.0 
#> ✓ purrr        0.3.4      ✓ yardstick    0.0.8
#> ── Conflicts ───────────────────────────────────────── tidymodels_conflicts() ──
#> x purrr::discard() masks scales::discard()
#> x dplyr::filter()  masks stats::filter()
#> x dplyr::lag()     masks stats::lag()
#> x recipes::step()  masks stats::step()
#> • Learn how to get started at https://www.tidymodels.org/start/
```

Clinial research always use linear regression, logistic regression, Cox
regression. Easystats make this much easier!

``` r
library(easystats)
#> # Attaching packages: easystats 0.4.2
#> ✔ insight     0.14.4     ✔ datawizard  0.2.1   
#> ✔ bayestestR  0.11.0     ✔ performance 0.8.0   
#> ✔ parameters  0.14.0.1   ✔ effectsize  0.5     
#> ✔ modelbased  0.7.0      ✔ correlation 0.7.1   
#> ✔ see         0.6.8      ✔ report      0.4.0
```
