# qqman: An R package for creating Q-Q and manhattan plots from GWAS results.

## Installation

Install the stable release from CRAN:

```coffee
# Install once
# Note, not yet available on CRAN
# install.packages("qqman")
# library(qqman)
```

Install the most recent development release with devtools (note, there be dragons here):

```coffee
# Install once with devtools
library(devtools)
install_github("stephenturner/qqman", ref="pkg")

# Load
library(qqman)
```

## Usage

See the package vignette for more examples:

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

* Submission to CRAN pending.
* This release is substantially simplified for the sake of maintainability and creating an R package. The old code that allows confidence intervals on the Q-Q plot and allows more flexible annotation and highlighting is still available at the [version 0.0.0 tag](https://github.com/stephenturner/qqman/tree/v0.0.0).
* Special thanks to Dan Capurso and Tim Knutsen for useful contributions and bugfixes.
* Thanks to all the blog commenters for pointing out bugs and other issues.
