Luciernaga
================
David Rach
<h4>  
README updated: <i>Aug-04-2024</i>  
</h4>

<!-- To modify Package/Title/Description/Authors fields, edit the DESCRIPTION file -->
<!-- badges: start -->

[![R build
status](https://github.com/DavidRach/Luciernaga/workflows/rworkflows/badge.svg)](https://github.com/DavidRach/Luciernaga/actions)
[![License: AGPL (\>=
3)](https://img.shields.io/badge/license-AGPL%20(%3E=%203)-blue.svg)](https://cran.r-project.org/web/licenses/AGPL%20(%3E=%203))
[![](https://img.shields.io/badge/devel%20version-0.1.0-black.svg)](https://github.com/DavidRach/Luciernaga)
[![](https://img.shields.io/github/languages/code-size/DavidRach/Luciernaga.svg)](https://github.com/DavidRach/Luciernaga)
[![](https://img.shields.io/github/last-commit/DavidRach/Luciernaga.svg)](https://github.com/DavidRach/Luciernaga/commits/master)
[![codecov](https://codecov.io/gh/DavidRach/Luciernaga/graph/badge.svg?token=GHWZ3NJ7IK)](https://codecov.io/gh/DavidRach/Luciernaga)
<br> <!-- badges: end -->

<img src="inst/hex/hex.png" width="50%" style="display: block; margin: auto;" />

## `Luciernaga`: Quality checks, signature interrogation, and data visualization for Spectral Flow Cytometry (SFC) unmixing controls.

Having good quality unmixing controls is critical to accurate unmixing
of full-stained samples in Spectral Flow Cytometry, but few tools exist
to check whether this is the case. In this package we provide tools for
quality control checks, signature interrogation and data visualization
to allow individuals responsible for the pre-processing to make informed
decisions. Tools are implemented for individual .fcs files at a single
timepoint or across an experimental run.

If you use `Luciernaga`, please cite:

<!-- Modify this by editing the file: inst/CITATION  -->

## Installation

``` r
if(!require("remotes")) install.packages("remotes")

remotes::install_github("https://github.com/DavidRach/Luciernaga")
library(Luciernaga)
```

## Documentation

### [Website](https://davidrach.github.io/Luciernaga)

### [Get started](https://davidrach.github.io/Luciernaga/articles/Luciernaga)

<br>
