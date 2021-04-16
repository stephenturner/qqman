# qqman: An R package for creating Q-Q and manhattan plots from GWAS results.

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/qqman)](https://cran.r-project.org/package=qqman)

[![DOI](https://joss.theoj.org/papers/10.21105/joss.00731/status.svg)](https://doi.org/10.21105/joss.00731)

![qqman.gif](tools/qqman.gif)

## Citation

If you'd like to cite qqman (appreciated but not required), please cite the publication below:

Turner, (2018). qqman: an R package for visualizing GWAS results using Q-Q and manhattan plots. _Journal of Open Source Software_, 3(25), 731, https://doi.org/10.21105/joss.00731

More details are found in the preprint:

Turner, S.D. qqman: an R package for visualizing GWAS results using Q-Q and manhattan plots. [*biorXiv* DOI: 10.1101/005165](https://biorxiv.org/content/early/2014/05/14/005165).

## Installation

Install the stable release from CRAN:

```coffee
install.packages("qqman")
```

Or install directly from github using devtools

```coffee
library(devtools)
install_github("stephenturner/qqman")
```

Or install the most recent development release with devtools (note, there be dragons here):

```coffee
library(devtools)
install_github("stephenturner/qqman", ref="dev")
```

Load the package each time you use it:

```coffee
library(qqman)
```

## Usage

See the [online package vignette](https://cran.r-project.org/web/packages/qqman/vignettes/qqman.html) for more examples:

```coffee
vignette("manhattan")
```

Take a look at the built-in data:

```coffee
head(gwasResults)
```

Basic manhattan plot using built-in data:

```coffee
manhattan(gwasResults)
```

Basic Q-Q plot using built-in data:

```coffee
qq(gwasResults$P)
```

Get help:

```coffee
?manhattan
?qq
```

## Notes

* This release is substantially simplified for the sake of maintainability and creating an R package. The old code that allows confidence intervals on the Q-Q plot and allows more flexible annotation and highlighting is still available at the [version 0.0.0 tag](https://github.com/stephenturner/qqman/tree/v0.0.0).
* Special thanks to Dan Capurso and Tim Knutsen for useful contributions and bugfixes.
* Thanks to all the blog commenters for pointing out bugs and other issues.
