# Stephen Turner
# http://StephenTurner.us/
# http://GettingGeneticsDone.blogspot.com/

# Daniel Capurso
# UCSF
# http://www.linkedin.com/in/dcapurso

# Last major update: June 10, 2013
# R code for making manhattan plots and QQ plots from plink output files. 

### This is for testing purposes. ######################################
#   set.seed(42)
#   nchr=22
#   nsnps=1000
#   test_data = data.frame(
#       SNP=sapply(1:(nchr*nsnps), function(x) paste("rs",x,sep='')),
#       CHR=rep(1:nchr,each=nsnps), 
#       BP=rep(1:nsnps,nchr), 
#       P=runif(nchr*nsnps)
#   )
#   ### d[d$SNP=='rs20762',]$P = 1e-29
#   top_snps = c('rs13895','rs20762')
#   surrounding_snps = list(as.character(test_data$SNP[13795:13995]),as.character(test_data$SNP[20662:20862]))
#   
#   pvector = test_data$P
#   names(pvector) = test_data$SNP
#   
#   #CALL:
#   #manhattan(test_data,annotate=top_snps,highlight=surrounding_snps)
#   #qq(pvector,annotate=top_snps,highlight=top_snps)
#########################################################################

# manhattan plot using base graphics
manhattan <- function(dataframe, limitchromosomes=NULL,pt.col=c('gray10','gray50'),pt.bg=c('gray10','gray50'),
    pt.cex=0.45,pch=21,cex.axis=0.7,gridlines=F,gridlines.col='gray83',gridlines.lty=1,gridlines.lwd=1,ymax=8, ymax.soft=T, annotate=NULL,annotate.cex=0.7,annotate.font=3,
    suggestiveline=-log10(1e-5), suggestiveline.col='blue', suggestiveline.lwd=1.5, suggestiveline.lty=1, 
    genomewideline=-log10(5e-8), genomewideline.col='red', genomewideline.lwd=1.5, genomewideline.lty=1, 
    highlight=NULL,highlight.col=c('green3','magenta'),highlight.bg=c('green3','magenta'),  ...) {
    #============================================================================================
    ######## Check data and arguments
    d = dataframe
    if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
    d = d[( !is.na(d$CHR) & !is.na(d$BP) & !is.na(d$P) ), ]
    
    if (TRUE %in% is.na(suppressWarnings(as.numeric(d$CHR)))) warning('non-numeric, non-NA entries in CHR column of dataframe. attempting to remove..')
    if (TRUE %in% is.na(suppressWarnings(as.numeric(d$BP)))) warning('non-numeric, non-NA entries in BP column of dataframe. attempting to remove..')
    if (TRUE %in% is.na(suppressWarnings(as.numeric(d$P)))) warning('non-numeric, non-NA entries in P column of dataframe. attempting to remove..')
    
    d = d[!is.na(suppressWarnings(as.numeric(d$CHR))),] # remove rows with non-numeric, non-NA entries
    d = d[!is.na(suppressWarnings(as.numeric(d$BP))),]
    d = d[!is.na(suppressWarnings(as.numeric(d$P))),]
    
    
    if (!is.null(annotate)){
        if ('SNP' %in% names(d)){
        	missing_annotate = annotate[!(annotate %in% d$SNP)]
        	annotate = annotate[annotate %in% d$SNP]
        	
        	if (length(missing_annotate)>0){
        		print('These SNPs were not annotated because of missing data:')
        		print(missing_annotate)
        	}
        } else {
            stop("D'oh! Dataframe must have a column $SNP with rs_ids to use annotate feature.")
        }
    }
    if (!is.numeric(annotate.cex) | annotate.cex<0) annotate.cex=0.7
    if (!is.numeric(annotate.font)) annotate.font=3
    
    if (is.character(gridlines.col[1]) & !(gridlines.col[1] %in% colors())) gridlines.col = 'gray83'
    if (!is.numeric(pt.cex) | pt.cex<0) pt.cex=0.45
    if (is.character(pt.col) & (FALSE %in% (pt.col %in% colors()))) pt.col = c('gray10','gray50')
    if (is.character(pt.bg) & (FALSE %in% (pt.bg %in% colors()))) pt.bg = F
    if (is.character(highlight.col) & (FALSE %in% (highlight.col %in% colors()))) highlight.col = c('green3','magenta')
    if (is.character(highlight.bg) & (FALSE %in% (highlight.bg %in% colors()))) highlight.bg = F
    if (is.character(suggestiveline.col[1]) & !(suggestiveline.col[1] %in% colors())) suggestiveline.col = 'blue'
    if (is.character(genomewideline.col[1]) & !(genomewideline.col[1] %in% colors())) genomewideline.col = 'red'
        
    if(!is.null(limitchromosomes)){
        if (TRUE %in% is.na(suppressWarnings(as.numeric(limitchromosomes)))){
            stop('limitchromosomes argument is not numeric') 
        } else {  
            d = d[d$CHR %in% as.numeric(limitchromosomes), ]
        }
    }
    

    ######################
    
    # Set positions, ticks, and labels for plotting
    d=subset(d[order(d$CHR, d$BP), ], (P>0 & P<=1)) # sort, and keep only 0<P<=1
    d$logp = -log10(d$P)
    d$pos=NA
    
    
    # Ymax
    if(is.na(suppressWarnings(as.numeric(ymax)))){  # not numeric
        ymax = ceiling(max(-log10(d$P)))
        warning('non-numeric ymax argument.')
    } else if (as.numeric(ymax) < 0){           # negative
        ymax = ceiling(max(-log10(d$P)))
        warning('negative ymax argument.')
    }
    if (ymax.soft==T){ #if soft, ymax is just the lower limit for ymax
        ymax = max(ymax, ceiling(max(-log10(d$P))))
        
        # make ymax larger if top annotate SNP is very high
        if (!is.null(annotate)){
            annotate.max = max(d[which(d$SNP %in% annotate),]$logp)
            if ((ymax - annotate.max) < 0.18*ymax){
                ymax = annotate.max + 0.18*ymax
            }
        }
    } #else, ymax = ymax
    
    ## Fix for the bug where one chromosome is missing. Adds index column #####
    d$index=NA
    ind = 0
    for (i in unique(d$CHR)){
        ind = ind + 1
        d[d$CHR==i,]$index = ind
    }
    ########
    
    nchr=length(unique(d$CHR))
    if (nchr==1) {
        d$pos=d$BP
        ticks=floor(length(d$pos))/2+1
        xlabel = paste('Chromosome',unique(d$CHR),'position')
        labs = ticks
    } else {
        ticks = rep(NA,length(unique(d$CHR))+1)
        ticks[1] = 0
        for (i in 1:max(d$index)) {
            d[d$index==i, ]$pos   =    (d[d$index==i, ]$BP - d[d$index==i,]$BP[1]) +1 +ticks[i]
            ticks[i+1] = max(d[d$index==i,]$pos)
        }
        xlabel = 'Chromosome'
        labs = append(unique(d$CHR),'')
    }
    
    # Initialize plot
    xmax = max(d$pos) * 1.03
    xmin = max(d$pos) * -0.03
    #ymax = ceiling(ymax * 1.03)
    ymin = -ymax*0.03
    plot(0,col=F,xaxt='n',bty='n',xaxs='i',yaxs='i',xlim=c(xmin,xmax), ylim=c(ymin,ymax),
            xlab=xlabel,ylab=expression(-log[10](italic(p))),las=1,cex.axis=cex.axis, ...)
    
    # stagger labels
    blank = rep('',length(labs))
    lowerlabs = rep('',length(labs))
    upperlabs = rep('',length(labs))
    
    for (i in 1:length(labs)){
        if (i %% 2 == 0){
            lowerlabs[i] = labs[i]
        } else{
            upperlabs[i] = labs[i]
        }
    }
    
    axis(1,at=ticks,labels=blank,lwd=0,lwd.ticks=1,cex.axis=cex.axis)
    axis(1,at=ticks,labels=upperlabs,lwd=0,lwd.ticks=0,cex.axis=cex.axis,line=-0.25)
    axis(1,at=ticks,labels=lowerlabs,lwd=0,lwd.ticks=0,cex.axis=cex.axis,line=0.25)
    
    yvals = par('yaxp')
    yinterval = par('yaxp')[2] / par('yaxp')[3]
    axis(2,at= (seq(0,(ymax+yinterval/2),yinterval) - yinterval/2),labels=F,lwd=0,lwd.ticks=1,cex.axis=cex.axis)
    
    # Gridlines
    if (isTRUE(gridlines)){
        
        abline(v=ticks,col=gridlines.col[1],lwd=gridlines.lwd,lty=gridlines.lty) #at ticks
        abline(h=seq(0,ymax,yinterval),col=gridlines.col[1],lwd=gridlines.lwd,lty=gridlines.lty) # at labeled ticks
        #abline(h=(seq(0,ymax,yinterval) - yinterval/2),col=gridlines.col[1],lwd=1.0) # at unlabeled ticks
    }
    
    # Points, with optional highlighting
    pt.col = rep(pt.col,max(d$CHR))[1:max(d$CHR)]
    pt.bg = rep(pt.bg,max(d$CHR))[1:max(d$CHR)]
    d.plain = d
    if (!is.null(highlight)) {
        if(class(highlight)!='character' & class(highlight)!='list'){
            stop('"highlight" must be a char vector (for 1 color) or list (for multi color).')
        }
        
        
        if ('SNP' %in% names(d)){
            missing_highlight = highlight[!(highlight %in% d$SNP)]
        	highlight = highlight[highlight %in% d$SNP]
        	
        	if (length(missing_highlight)>0){
        		print('These SNPs were not highlightd because of missing data:')
        		print(missing_highlight)
        	}
        
        } else {
            stop("D'oh! Dataframe must have a column $SNP with rs_ids to use highlight feature.")
        }
        
        if (class(highlight)=='character'){ #if char vector, make list for consistency in plotting below
            highlight = list(highlight)
        }
        
        highlight.col = rep(highlight.col,length(highlight))[1:length(highlight)]
        highlight.bg = rep(highlight.bg,length(highlight))[1:length(highlight)]
        
        for (i in 1:length(highlight)){
            d.plain = d.plain[which(!(d.plain$SNP %in% highlight[[i]])), ]
        }
    }
    
    icol=1
    for (i in unique(d.plain$CHR)) {
        with(d.plain[d.plain$CHR==i, ],points(pos, logp, col=pt.col[icol],bg=pt.bg[icol],cex=pt.cex,pch=pch,...))
        icol=icol+1
    }
    
    if (!is.null(highlight)){   
        for (i in 1:length(highlight)){
            d.highlight=d[which(d$SNP %in% highlight[[i]]), ]
            with(d.highlight, points(pos, logp, col=highlight.col[i],bg=highlight.bg[i],cex=pt.cex,pch=pch,...)) 
        }
    }
    
    # Significance lines
    if (is.numeric(suggestiveline)) abline(h=suggestiveline, col=suggestiveline.col[1],lwd=suggestiveline.lwd,lty=suggestiveline.lty)
    if (is.numeric(genomewideline)) abline(h=genomewideline, col=genomewideline.col[1],lwd=genomewideline.lwd,lty=genomewideline.lty)

    # Annotate
    if (!is.null(annotate)){
        d.annotate = d[which(d$SNP %in% annotate),]
        text(d.annotate$pos,(d.annotate$logp + 0.019*ymax),labels=d.annotate$SNP,srt=90,cex=annotate.cex,adj=c(0,0.48),font=annotate.font)      
    }

    # Box
    box()
}







## Make a pretty QQ plot of p-values ### Add ymax.soft
qq = function(pvector,gridlines=F,gridlines.col='gray83',gridlines.lwd=1,gridlines.lty=1,confidence=T,confidence.col='gray81',
    pt.cex=0.5,pt.col='black',pt.bg='black',pch=21,abline.col='red',abline.lwd=1.8,abline.lty=1,ymax=8,ymax.soft=T,
    highlight=NULL,highlight.col=c('green3','magenta'),highlight.bg=c('green3','magenta'),
    annotate=NULL,annotate.cex=0.7,annotate.font=3,cex.axis=0.95,...) {
    #======================================================================================================
    ######## Check data and arguments; create observed and expected distributions
    d = suppressWarnings(as.numeric(pvector))
    names(d) = names(pvector)
    d = d[!is.na(d)] # remove NA, and non-numeric [which were converted to NA during as.numeric()]
    d = d[d>0 & d<1] # only Ps between 0 and 1
    
    
    if (!is.null(highlight) | !is.null(annotate)){
        if (is.null(names(d))) stop("P-value vector must have names to use highlight or annotate features.")
        d = d[!is.na(names(d))]
        if (!is.null(highlight) & FALSE %in% (highlight %in% names(d))) stop ("D'oh! Highlight vector must be a subset of names(pvector).")
        if (!is.null(annotate) & FALSE %in% (annotate %in% names(d))) stop ("D'oh! Annotate vector must be a subset of names(pvector).")
    }
    
    d = d[order(d,decreasing=F)] # sort
    o = -log10(d)
    e = -log10( ppoints(length(d) ))
    if (!is.null(highlight) | !is.null(annotate)) names(e) = names(o) = names(d)
    
    if (!is.numeric(ymax) | ymax<max(o)) ymax <- max(o) 
    
    if (!is.numeric(pt.cex) | pt.cex<0) pt.cex=0.5
    if (!is.numeric(annotate.cex) | annotate.cex<0) annotate.cex=0.7
    if (!is.numeric(annotate.font)) annotate.font=3
    
    if (is.character(gridlines.col[1]) & !(gridlines.col[1] %in% colors())) gridlines.col = 'gray83'
    if (is.character(confidence.col[1]) & !(confidence.col[1] %in% colors())) confidence.col = 'gray81'
    if (is.character(abline.col[1]) & !(abline.col[1] %in% colors())) abline.col = 'red'
    
    if (FALSE %in% (pt.col %in% colors() | !is.na(suppressWarnings(as.numeric(pt.col))) )){
        pt.col = 'black'; warning("pt.col argument(s) not recognized. Setting to default: 'black'.")
    }

    if (FALSE %in% (pt.bg %in% colors() | !is.na(suppressWarnings(as.numeric(pt.bg))) )){
        pt.bg = 'black'; warning("pt.bg argument(s) not recognized. Setting to default: 'black'.")
    }
    
    if (FALSE %in% (highlight.col %in% colors() | !is.na(suppressWarnings(as.numeric(highlight.col))) )){
        highlight.col = 'blue'; warning("highlight.col argument(s) not recognized. Setting to default: 'blue'.")
    }

    if (FALSE %in% (highlight.bg %in% colors() | !is.na(suppressWarnings(as.numeric(highlight.bg))) )){
        highlight.bg = 'blue'; warning("highlight.bg argument(s) not recognized. Setting to default: 'blue'.")
    }
    
    # Ymax
    if(is.na(suppressWarnings(as.numeric(ymax)))){  # not numeric
        ymax = ceiling(max(o))
        warning('non-numeric ymax argument.')
    } else if (as.numeric(ymax) < 0){           # negative
        ymax = ceiling(max(o))
        warning('negative ymax argument.')
    }
    if (ymax.soft==T){ #if soft, ymax is just the lower limit for ymax
        ymax = max(ymax, ceiling(max(o)))
    } #else, ymax = ymax
            
    
    ################################
    
    # Initialize plot
    #print('Setting up plot.')
    #print(ymax)
    xspace = 0.078
    xmax = max(e) * 1.019
    xmin = max(e) * -0.035
    #ymax = ceiling(ymax * 1.03)
    ymin = -ymax*0.03
    plot(0,xlab=expression(Expected~~-log[10](italic(p))),ylab=expression(Observed~~-log[10](italic(p))),
            col=F,las=1,xaxt='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty='n',xaxs='i',yaxs='i',cex.axis=cex.axis, ...)
    axis(side=1,labels=seq(0,max(e),1),at=seq(0,max(e),1),cex.axis=cex.axis,lwd=0,lwd.ticks=1)
    
    # Grid lines
    if (isTRUE(gridlines)){
        yvals = par('yaxp')
        yticks = seq(yvals[1],yvals[2],yvals[2]/yvals[3])
        abline(v=seq(0,max(e),1),col=gridlines.col[1],lwd=gridlines.lwd,lty=gridlines.lty)
        abline(h=yticks,col=gridlines.col[1],lwd=gridlines.lwd,lty=gridlines.lty)
    }
    
     #Confidence intervals
     find_conf_intervals = function(row){
        i = row[1]
        len = row[2]
        if (i < 10000 | i %% 100 == 0){
            return(c(-log10(qbeta(0.95,i,len-i+1)), -log10(qbeta(0.05,i,len-i+1))))
        } else { # Speed up
            return(c(NA,NA))
        }
     }

     # Find approximate confidence intervals
    if (isTRUE(confidence)){
        #print('Plotting confidence intervals.')
        ci = apply(cbind( 1:length(e), rep(length(e),length(e))), MARGIN=1, FUN=find_conf_intervals)
        bks = append(seq(10000,length(e),100),length(e)+1)
        for (i in 1:(length(bks)-1)){
            ci[1, bks[i]:(bks[i+1]-1)] = ci[1, bks[i]]
            ci[2, bks[i]:(bks[i+1]-1)] = ci[2, bks[i]]
        }
        colnames(ci) = names(e)
        # Extrapolate to make plotting prettier (doesn't affect intepretation at data points)
        slopes = c((ci[1,1] - ci[1,2]) / (e[1] - e[2]), (ci[2,1] - ci[2,2]) / (e[1] - e[2]))
        extrap_x = append(e[1]+xspace,e) #extrapolate slightly for plotting purposes only
        extrap_y = cbind( c(ci[1,1] + slopes[1]*xspace, ci[2,1] + slopes[2]*xspace), ci)
        
        polygon(c(extrap_x, rev(extrap_x)), c(extrap_y[1,], rev(extrap_y[2,])),col = confidence.col[1], border = confidence.col[1]) 
    }
    
    # Points (with optional highlighting)
    #print('Plotting data points.')
    fills = rep(pt.bg,length(o))
    borders = rep(pt.col,length(o))
    names(fills) = names(borders) = names(o)
    if (!is.null(highlight)){   
        borders[highlight] = rep(NA,length(highlight))
        fills[highlight] = rep(NA,length(highlight))
    }
    points(e,o,pch=pch,cex=pt.cex,col=borders,bg=fills)
    
    if (!is.null(highlight)){
        points(e[highlight],o[highlight],pch=pch,cex=pt.cex,col=highlight.col,bg=highlight.bg)
    }
    
    #Abline
    abline(0,1,col=abline.col,lwd=abline.lwd,lty=abline.lty)
    
    # Annotate SNPs
    if (!is.null(annotate)){
        x = e[annotate] # x will definitely be the same
        y = -0.1 + apply(rbind(o[annotate],ci[1,annotate]),2,min)
        text(x,y,labels=annotate,srt=90,cex=annotate.cex,adj=c(1,0.48),font=annotate.font)      
    }
    # Box
    box()
}








# Old ggplot2 code --------------------------------------------------------
#
## manhattan plot using ggplot2
# gg.manhattan = function(dataframe, title=NULL, max.y="max", suggestiveline=0, genomewideline=-log10(5e-8), size.x.labels=9, size.y.labels=10, annotate=F, SNPlist=NULL) {
# library(ggplot2)
#     if (annotate & is.null(SNPlist)) stop("You requested annotation but provided no SNPlist!")
#   d=dataframe
#   #limit to only chrs 1-23?
#   d=d[d$CHR %in% 1:23, ]
#   if ("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d) ) {
#       d=na.omit(d)
#       d=d[d$P>0 & d$P<=1, ]
#       d$logp = -log10(d$P)
#       d$pos=NA
#       ticks=NULL
#       lastbase=0
#       #new 2010-05-10
#       numchroms=length(unique(d$CHR))
#       if (numchroms==1) {
#           d$pos=d$BP
#       } else {
#       
#           for (i in unique(d$CHR)) {
#               if (i==1) {
#                   d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
#               }   else {
#                   lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
#                   d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
#               }
#               ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
#           }
#           ticklim=c(min(d$pos),max(d$pos))
# 
#       }
#       mycols=rep(c("gray10","gray60"),max(d$CHR))
#       if (max.y=="max") maxy=ceiling(max(d$logp)) else maxy=max.y
#       if (maxy<8) maxy=8
#       if (annotate) d.annotate=d[as.numeric(substr(d$SNP,3,100)) %in% SNPlist, ]
#       if (numchroms==1) {
#           plot=qplot(pos,logp,data=d,ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"))
#       }   else {
#           plot=qplot(pos,logp,data=d, ylab=expression(-log[10](italic(p))) , colour=factor(CHR))
#           plot=plot+scale_x_continuous(name="Chromosome", breaks=ticks, labels=(unique(d$CHR)))
#           plot=plot+scale_y_continuous(limits=c(0,maxy), breaks=1:maxy, labels=1:maxy)
#           plot=plot+scale_colour_manual(value=mycols)
#       }
#       if (annotate)   plot=plot + geom_point(data=d.annotate, colour=I("green3")) 
#       plot=plot + opts(legend.position = "none") 
#       plot=plot + opts(title=title)
#       plot=plot+opts(
#           panel.background=theme_blank(), 
#           panel.grid.minor=theme_blank(),
#           axis.text.x=theme_text(size=size.x.labels, colour="grey50"), 
#           axis.text.y=theme_text(size=size.y.labels, colour="grey50"), 
#           axis.ticks=theme_segment(colour=NA)
#       )
#       if (suggestiveline) plot=plot+geom_hline(yintercept=suggestiveline,colour="blue", alpha=I(1/3))
#       if (genomewideline) plot=plot+geom_hline(yintercept=genomewideline,colour="red")
#       plot
#   }   else {
#       stop("Make sure your data frame contains columns CHR, BP, and P")
#   }
# }
# 
# ## QQ plot using ggplot2
# gg.qq = function(pvector, title=NULL, spartan=F) {
#   library(ggplot2)
#   o = -log10(sort(pvector,decreasing=F))
#   #e = -log10( 1:length(o)/length(o) )
#   e = -log10( ppoints(length(pvector) ))
#   plot=qplot(e,o, xlim=c(0,max(e)), ylim=c(0,max(o))) + stat_abline(intercept=0,slope=1, col="red")
#   plot=plot+opts(title=title)
#   plot=plot+scale_x_continuous(name=expression(Expected~~-log[10](italic(p))))
#   plot=plot+scale_y_continuous(name=expression(Observed~~-log[10](italic(p))))
#   if (spartan) plot=plot+opts(panel.background=theme_rect(col="grey50"), panel.grid.minor=theme_blank())
#   plot
# }
# 
# ## Make a qq and manhattan plot
# gg.qqman = function(data="plinkresults") {
#   myqqplot = ggqq(data$P)
#   mymanplot = ggmanhattan(data)
#   ggsave(file="qqplot.png",myqqplot,w=5,h=5,dpi=100)
#   ggsave(file="manhattan.png",mymanplot,width=12,height=9,dpi=100)
# }
# 
# ## make qq and manhattan plots for a list of files
# gg.qqmanall= function(command="ls *assoc") {
#   filelist=system(command,intern=T)
#   datalist=NULL
#   for (i in filelist) {datalist[[i]]=read.table(i,T)}
#   highestneglogp=ceiling(max(sapply(datalist, function(df) max(na.omit(-log10(df$P))))))
#   print(paste("Highest -log10(P) = ",highestneglogp),quote=F)
#   start=Sys.time()
#   for (i in names(datalist)) {
#       myqqplot=ggqq(datalist[[i]]$P, title=i)
#       ggsave(file=paste("qqplot-",    i, ".png", sep=""),myqqplot, width=5, height=5,dpi=100)
#       mymanplot=ggmanhattan(datalist[[i]], title=i, max.y=highestneglogp)
#       ggsave(file=paste("manhattan-", i, ".png", sep=""),mymanplot,width=12,height=9,dpi=100)
#   }
#   end=Sys.time()
#   print(elapsed<-end-start)
# }
