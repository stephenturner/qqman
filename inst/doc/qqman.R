
## ----, include=FALSE-----------------------------------------------------
library(qqman)
library(knitr)
opts_chunk$set(comment=NA, fig.width=12, fig.height=9, message=FALSE)


## ----generatedata, eval=FALSE, echo=FALSE--------------------------------
## # This code used to generate the test data. Runs slow, but does the job.
## chrstats <- data.frame(chr=1:22, nsnps=1500)
## chrstats$nsnps <- with(chrstats, round(nsnps/chr^(1/3)))
## chrstats
## 
## d <- data.frame(SNP=rep(NA, sum(chrstats$nsnps)),
##                 CHR=rep(NA, sum(chrstats$nsnps)),
##                 BP=rep(NA, sum(chrstats$nsnps)),
##                 P=rep(NA, sum(chrstats$nsnps)))
## snpi <- 1
## set.seed(42)
## for (i in chrstats$chr) {
##     for (j in 1:chrstats[i, 2]) {
##         d[snpi, ]$SNP=paste0("rs", snpi)
##         d[snpi, ]$CHR=i
##         d[snpi, ]$BP=j
##         d[snpi, ]$P=runif(1)
##         snpi <- snpi+1
##     }
## }
## 
## divisor <- c(seq(2,50,2), seq(50,2,-2))
## divisor <- divisor^4
## length(divisor)
## d[3026:3075, ]$P <- d[3026:3075, ]$P/divisor
## snpsOfInterest <- paste0("rs", 3001:3100)
## qq(d$P)
## manhattan(d, highlight=snpsOfInterest)
## gwasResults <- d
## save(gwasResults, file="data/gwasResults.RData")


## ------------------------------------------------------------------------
str(gwasResults)
head(gwasResults)
tail(gwasResults)


## ------------------------------------------------------------------------
as.data.frame(table(gwasResults$CHR))


## ------------------------------------------------------------------------
manhattan(gwasResults)


## ------------------------------------------------------------------------
manhattan(gwasResults, main="Manhattan Plot", cex=0.5, cex.axis=0.8)


## ------------------------------------------------------------------------
manhattan(gwasResults, col=c("blue4", "orange3"), ymax=12)


## ------------------------------------------------------------------------
manhattan(gwasResults, suggestiveline=F, genomewideline=F)


## ------------------------------------------------------------------------
manhattan(subset(gwasResults, CHR==1))


## ------------------------------------------------------------------------
str(snpsOfInterest)
manhattan(gwasResults, highlight=snpsOfInterest)


## ------------------------------------------------------------------------
manhattan(subset(gwasResults, CHR==3), highlight=snpsOfInterest, main="Chr 3")


## ------------------------------------------------------------------------
qq(gwasResults$P, main="Q-Q plot of GWAS p-values")


