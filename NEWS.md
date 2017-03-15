# qqman 0.1.4

* Minor fix to location to referenced image for pandoc self-contained README.html generation. 

# qqman 0.1.3

* Annotate SNPs below a p-value threshold with the `annotatePval=` option. See vignette for details.
* Annotate the top SNP on each chromosome with the `annotateTop=` option. See vignette for details.

# qqman 0.1.2

* Does not assume that SNPs are evenly distributed across chromosomes when deciding where to place the tick in the center of the chromosome.
* Changed single chromosome x-axis notation to use Mb instead of raw pos
* `qq()` accepts graphical parameters the same way as `manhattan()`
* Removed default `xlim`
* Citation details on package load
* Added axis label options
* Removed `ymax` argument in favor of allowing user to set `ylim` in `...`
* Option to *not* take log of p-value

# qqman 0.1.1

* Fixed a bunch of typos in the vignette
