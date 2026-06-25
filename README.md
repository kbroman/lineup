### R/lineup: detecting and correcting sample mix-ups <a href="https://github.com/kbroman/lineup"><img src="https://kbroman.org/lineup/lineup_logo.png" align="right" height="138" alt="R/lineup logo"/></a>

[![R-CMD-check](https://github.com/kbroman/lineup/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kbroman/lineup/actions/workflows/R-CMD-check.yaml)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/lineup)](https://cran.r-project.org/package=lineup)
[![r-universe badge](https://kbroman.r-universe.dev/lineup/badges/version)](https://kbroman.r-universe.dev/lineup)
[![zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4277411.svg)](https://doi.org/10.5281/zenodo.4277411)

[R/lineup](https://github.com/kbroman/lineup) is an
[R](https://www.r-project.org) package with tools for detecting and
correcting sample mix-ups between two sets of measurements, such as
between gene expression data on two tissues.

This is particularly aimed at eQTL data for an experimental cross.

---

### Installation

Install the R/lineup package from [CRAN](https://cran.r-project.org):

```r
install.packages("lineup")
```

Alternatively, install it from [R
universe](https://kbroman.r-universe.dev):

```r
install.packages("lineup", repos=c("https://kbroman.r-universe.dev",
                                   "https://cloud.r-project.org"))
```

Or use [remotes](https://remotes.r-lib.org) to install it from its GitHub source:

```r
install.packages("remotes")
remotes::install_github("kbroman/lineup")
```

---

### Vignette

A vignette describing the use of the package is available
[on the web](https://kbroman.org/lineup/lineup.html).
Or view it from within R by loading the package and then using the
`vignette()` function.

```r
library(lineup2)
vignette("lineup", package="lineup")
```


---

### Citation

To cite R/lineup in publications use

- Broman KW, Keller MP, Broman AT, Kendziorski C, Yandell BS, Sen Ś,
  Attie AD (2015) Identification and correction of
  sample mix-ups in expression genetic data: A case study.
  G3 5:2177-2186
  [doi:10.1534/g3.115.019778](https://doi.org/10.1534/g3.115.019778)


---

### License

This package is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License, version 3, as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the GNU
General Public License for more details.

A copy of the GNU General Public License, version 3, is available at
<https://www.r-project.org/Licenses/GPL-3>
