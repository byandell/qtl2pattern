### R/qtl2pattern

[Karl Broman](http://kbroman.org) & [Brian Yandell](http://www.stat.wisc.edu/~yandell)

R/qtl2pattern is a support library for extending [R/qtl2](http://kbroman.org/qtl2), primarily associated with capabilities.reimplementation of [qtl2ggplot](https://github.com/byandell/qtl2ggplot) (for data visualization). See
[R/qtl2](http://kbroman.org/qtl2) (aka qtl2) for the bigger story of the qtl2 suite of routines.

---

### Installation

R/qtl2 is early in development and so is not yet available on
[CRAN](http://cran.r-project.org).

You can install R/qtl2 from [GitHub](https://github.com/rqtl).

You first need to install the
[devtools](https://github.com/hadley/devtools) package, plus a set of
package dependencies: [yaml](https://cran.r-project.org/package=yaml),
[jsonlite](https://cran.r-project.org/package=jsonlite),
[data.table](https://cran.r-project.org/package=data.table),
and [RcppEigen](https://github.com/RcppCore/RcppEigen).
(Additional, secondary dependencies will also be installed)

    install.packages(c("devtools","yaml","jsonlite","data.table","RcppEigen","fst"))

You will also need the following packages for qtl2ggplot:

    install.packages(c("tidyverse", "RColorBrewer"))

Then, install R/qtl2 using `devtools::install_github()`.

    library(devtools)
    install_github(c("rqtl/qtl2", "rqtl/qtl2fst"))

Once you have installed these, install qtl2pattern as

    install_github("byandell/qtl2pattern")
    
To install vignettes:

    install_github("byandell/qtl2pattern", build_vignettes = TRUE)

---

### Vignettes

- [qtl2pattern](https://github.com/byandell/qtl2pattern/blob/master/vignettes/qtl2pattern.Rmd)
- [effectplot](https://github.com/byandell/qtl2pattern/blob/master/vignettes/effectplot.Rmd)

---

#### License

[Licensed](License.md) under [GPL-3](http://www.r-project.org/Licenses/GPL-3).
