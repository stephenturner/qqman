#' Creates a Q-Q plot
#' 
#' Creates a quantile-quantile plot from p-values from a GWAS study.
#' 
#' @param x A numeric vector of p-values OR a data.frame with columns "BP," "CHR," "P," and optionally, "SNP."
#' @param p A string denoting the column name for the p-value. Defaults to 
#'   PLINK's "P." Said column must be numeric.
#' @param snp A string denoting the column name for the SNP name (rs number). 
#'   Defaults to PLINK's "SNP." Said column should be a character.
#' @param highlight A character vector of SNPs in your dataset to highlight. 
#'   These SNPs should all be in your dataset.
#' @param col A character vector indicating which colors to alternate.
#' @param ... Other arguments passed to \code{plot()}
#' 
#' @return A Q-Q plot.
#' 
#' @keywords visualization qq qqplot
#' 
#' @examples
#' qq(gwasResults$P)
#' 
#' @export

qq = function(x, p="P", snp="SNP", highlight=NULL, col="black", ...) {
    
    # Check for sensible input
    if (is.numeric(x)) {
	# convert to data.frame
	x=data.frame(snp=NA,p=x)
	names(x)=c(snp,p)
    } else if (is.data.frame(x)) {
        if (!(p %in% names(x))) stop(paste("Column", p, "not found!"))
    	if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
    } else {
	stop("Input must be numeric or a gwasResults data.frame.")
    }

    # Create a new data.frame with columns called SNP (if available) and P,
    # limited to not missing, not nan, not null, not infinite, between 0 and 1
    d <- x[!is.na(x[[p]]) & !is.nan(x[[p]]) & !is.null(x[[p]]) & is.finite(x[[p]]) & x[[p]]<1 & x[[p]]>0,which(names(x) %in% c(snp,p))]
    # Order ascending by p-value
    d <- d[order(d[[p]]),]
    d$o = -log10(d[[p]])
    d$e = -log10( ppoints(nrow(d) ))

    # Highlight snps from a character vector
    if (!is.null(highlight)) {
        if (any(!(highlight %in% d$SNP))) warning("You're trying to highlight SNPs that don't exist in your results.")
        d.highlight=d$SNP %in% highlight
    } else {
        d.highlight = rep(F, nrow(d))
    }

    attach(d)

#     # The old way
#     plot(e, o, pch=20, 
#          xlab=expression(Expected~~-log[10](italic(p))), 
#          ylab=expression(Observed~~-log[10](italic(p))), 
#          ...)
    
    # The new way to initialize the plot.
    ## See http://stackoverflow.com/q/23922130/654296
    ## First, define your default arguments
    def_args <- list(pch=20, xlim=c(0, max(e)), ylim=c(0, max(o)), 
		     col=ifelse(d.highlight,"green3", col),
                     xlab=expression(Expected~~-log[10](italic(p))), 
                     ylab=expression(Observed~~-log[10](italic(p)))
    )
    ## Next, get a list of ... arguments
    #dotargs <- as.list(match.call())[-1L]
    dotargs <- list(...)
    ## And call the plot function passing NA, your ... arguments, and the default
    ## arguments that were not defined in the ... arguments.
    tryCatch(do.call("plot", c(list(x=e, y=o), def_args[!names(def_args) %in% names(dotargs)], dotargs)), warn=stop)

    detach(d)

    # Add diagonal
    abline(0,1,col="red")
}
