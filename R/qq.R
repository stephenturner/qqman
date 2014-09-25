#' Creates a Q-Q plot
#' 
#' Creates a quantile-quantile plot from p-values from a GWAS study.
#' 
#' @param pvector A numeric vector of p-values.
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

qq = function(pvector, ...) {
    
    # Check for sensible input
    if (!is.numeric(pvector)) stop("Input must be numeric.")
    
    # limit to not missing, not nan, not null, not infinite, between 0 and 1
    pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector<1 & pvector>0]
    
    # Observed and expected
    o = -log10(sort(pvector,decreasing=FALSE))
    e = -log10( ppoints(length(pvector) ))
    
    
#     # The old way
#     plot(e, o, pch=20, 
#          xlab=expression(Expected~~-log[10](italic(p))), 
#          ylab=expression(Observed~~-log[10](italic(p))), 
#          ...)
    
    # The new way to initialize the plot.
    ## See http://stackoverflow.com/q/23922130/654296
    ## First, define your default arguments
    def_args <- list(pch=20, xlim=c(0, max(e)), ylim=c(0, max(o)), 
                     xlab=expression(Expected~~-log[10](italic(p))), 
                     ylab=expression(Observed~~-log[10](italic(p)))
    )
    ## Next, get a list of ... arguments
    #dotargs <- as.list(match.call())[-1L]
    dotargs <- list(...)
    ## And call the plot function passing NA, your ... arguments, and the default
    ## arguments that were not defined in the ... arguments.
    tryCatch(do.call("plot", c(list(x=e, y=o), def_args[!names(def_args) %in% names(dotargs)], dotargs)), warn=stop)

    # Add diagonal
    abline(0,1,col="red")
    
}
