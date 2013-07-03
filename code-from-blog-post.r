## Source the qqman functions
download.file("https://raw.github.com/stephenturner/qqman/master/qqman.r", destfile="./qqman.r", method="curl")
source("./qqman.r")

## Load some example GWAS data, have a look at relevant columns
download.file("https://raw.github.com/stephenturner/qqman/master/plink.assoc.txt.gz", destfile="./plink.assoc.txt.gz", method="curl")
results <- read.table(gzfile("./plink.assoc.txt.gz"), header=T)
head(subset(results, select=c(SNP, CHR, BP, P)))

## Draw a manhattan plot
manhattan(results)

## Draw a manhattan plot changing the colors, and remove the suggestive and genomewide lines
manhattan(results, pt.col=c("black","grey50","darkorange1"), pch=20, genomewideline=F, suggestiveline=F)

## Download some SNPs to annotate, then create a manhattan plot, highlighting the SNPs of interest
download.file("https://raw.github.com/stephenturner/qqman/master/snps.txt", destfile="./snps.txt", method="curl")
snps_to_highlight <- scan("./snps.txt", character())
manhattan(results, highlight=snps_to_highlight, pch=20, main="Manhattan Plot")
manhattan(subset(results, CHR==11), pch=20, annotate=snps_to_highlight, main="Chromosome 11")

## Make a q-q plot of the p-values
qq(results$P)