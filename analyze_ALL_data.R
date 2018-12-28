
# Most of these codes are copied from https://github.com/carlosproca/gene-expr-norm-paper
# The data is accessible through GEO Series accession number GSE13425.

setwd("C:/ALL-data-analysis")

source( "normalize_median_condec.r" )
source( "normalize_stdvec_condec.r" )
source( "select_nodup.r" )
source( "watson_u2.r" )

ma.data.dir <- "C:/ALL-data-analysis/GSE13425_RAW"   # raw data directory

stopifnot( ma.data.dir != "" )

# define experimental conditions and samples

ma.treatment <- c( "E2A.rearranged.EP" )
ma.control <- c( "BCR.ABL" )
ma.condition <- c( ma.treatment, unique( ma.control ) )

ma.condition.sample.num <- c( 8, 4 )
ma.sample.condition <- rep( ma.condition, ma.condition.sample.num )
ma.sample.postfix <- unlist( lapply( ma.condition.sample.num, 
                                     function( mcsn ) 1:mcsn ) )
ma.sample <- paste0( ma.sample.condition, '.', ma.sample.postfix )

ma.sample.num <- length( ma.sample )
ma.condition.num <- length( ma.condition )
ma.treatment.num <- length( ma.treatment )

library(GEOquery)
gpl96 <- getGEO(filename='GPL96.soft')
ma.db <- Table(gpl96) 

library(tidyr)
ma.annot <- Table(gpl96) %>% dplyr::select("ID","Gene Symbol", "Representative Public ID" )

names( ma.annot ) <- c("probe.id", "SYMBOL", "GENENAME")

ma.probe.platform <- ma.annot$probe.id
rownames(ma.annot) <- ma.probe.platform
# read raw expression data

ma.data.def.file.name <- "ALL_def.txt"

ma.data.def <- read.table( ma.data.def.file.name, header=T )

ma.sample.file.name.read <- ma.data.def$FileName
ma.sample.read <- ma.data.def$Sample
ma.sample.condition.read <- ma.data.def$Condition

library( affy )
ma.raw.data <- ReadAffy( filenames = file.path( ma.data.dir, ma.data.def$FileName ) )

# get background-corrected data

ma.expr.data <- log2( exprs( mas5( ma.raw.data, normalize=FALSE ) ) )

ma.expr.datafull <- ma.expr.data
# remove expression data below noise threshold

ma.mas5.call <- exprs( mas5calls( ma.raw.data ) )
apply(ma.mas5.call,2,table)
ma.expr.data[ !ma.mas5.call== "P"] <- NA

# remove control probes and probes with missing data
# reorder probes in expression data

ma.probe.control.bol <- grepl( "^AFFX-", ma.probe.platform )

ma.probe.missing.data.bol <- rowSums( is.na( ma.expr.data[,1:8] ) ) > 6 | rowSums( is.na( ma.expr.data[,9:12] ) ) > 2

ma.probe <- rownames(ma.expr.data)[ ! ma.probe.control.bol &
                                      ! ma.probe.missing.data.bol ]
present <- apply(ma.mas5.call, 1, function(x)(sum(x == "P")))
table(present)

#ma.expr.data <- ma.expr.data[ ma.probe, ]

ma.expr.data <- ma.expr.datafull[ ma.probe, ]
ma.annot <- ma.annot[ma.probe,]
ma.probe.num <- length( ma.probe )

# plotDensity(log2(exprs(ma.raw.data))[present == 12, ], main="Before normalization (log2-transformed) - Probe level", xlim=c(4,16)); grid()
# plotDensity(log2(exprs(ma.raw.data))[present == 0, ], main="Before normalization (log2-transformed) - Probe level", xlim=c(4,16)); grid()

#ma.expr.data <- ma.expr.data[ ma.probe, ]
# reorder and rename samples in expression data

colnames( ma.expr.data ) <- ma.sample



# identify samples and conditions for analysis

set.treatment <- ma.treatment
set.control <- ma.control
set.condition <- ma.condition

set.sample <- ma.sample
set.sample.condition <- ma.sample.condition

set.treatment.num <- length( set.treatment )
set.condition.num <- length( set.condition )
set.sample.num <- length( set.sample )


# define comparisons between experimental conditions

treatment.control.comparison <- cbind( set.treatment, set.control )

set.comparison <- rbind( treatment.control.comparison )

colnames( set.comparison ) <- NULL

set.comparison.name <- apply( set.comparison, 1, paste0, collapse=".vs." )
set.comparison.num <- length( set.comparison.name )


# select expression data for analysis

set.expr.data <- ma.expr.data[ , ma.sample.condition %in% set.condition ]


# normalize expression data

random.seed <- 5000000
set.seed( random.seed )

norm.method <- c("median.condec")

norm.method.pch <- c( 16, 17)
norm.method.lty <- c( 1, 2)
norm.method.col <- c( "black", "red3")

bg.col.offset <- c( 0.8, 0.8, 0.8, 0 )
norm.method.bg.col <- c( "orange3", "black")

names( norm.method.pch ) <- norm.method
names( norm.method.lty ) <- norm.method
names( norm.method.col ) <- norm.method
names( norm.method.bg.col ) <- norm.method

for ( nmeth in norm.method )
{
  if ( nmeth == "stdvec.condec" )
  {
    stdvec.condec.norm.result <- normalize.stdvec.condec( set.expr.data, 
                                                          set.sample.condition, verbose=TRUE )
    stdvec.condec.norm.data <- stdvec.condec.norm.result$data
  }
  else if ( nmeth == "median.condec" )
  {
    median.condec.norm.result <- normalize.median.condec( set.expr.data, 
                                                          set.sample.condition, convergence.threshold = c( 0.001, 0.3 ), 
                                                          verbose=TRUE )
    median.condec.norm.data <- median.condec.norm.result$data
  }
  else
    stop( "wrong normalization method ", nmeth )
}

#####################################################
x11(width=8, height=10) ## Open a graphical window with specific dimensions
par(mfrow=c(2,1)) ## Share this window between two plots
boxplot(log2(exprs(ma.raw.data)), pch=".", las=2, cex.axis=0.5, main="Before normalization (log2-transformed) - Probe level")
grid()
boxplot(median.condec.norm.data, pch=".", las=2, cex.axis=0.5, main="RMA-normalized - Probeset level")
grid()

## Plot the density distributions before and after normalization
plotDensity(log2(exprs(ma.raw.data)), main="Before normalization (log2-transformed) - Probe level", xlim=c(-5,16)); grid()
plotDensity(median.condec.norm.data, main="RMA-normalized - Probeset level", xlim=c(-5,16)); grid()

par(mfrow=c(1,1)) ## Restore single plot per page

## Compute a median for each row (probe)
m <- rowMedians(median.condec.norm.data)

## Centring: substract the median value of each probe
rle <- sweep(median.condec.norm.data,1, m, "-")

## Plot a box of the probe-wise centred values
x11(width=8, height=8)
boxplot(rle, pch=".", las=2, cex.axis=0.5)

###################################################
# x11(width=8, height=10) ## Open a graphical window with specific dimensions
# par(mfrow=c(2,1)) ## Share this window between two plots
# boxplot(log2(exprs(ma.raw.data)), pch=".", las=2, cex.axis=0.5, main="Before normalization (log2-transformed) - Probe level")
# grid()
# boxplot(stdvec.condec.norm.data, pch=".", las=2, cex.axis=0.5, main="RMA-normalized - Probeset level")
# grid()

## Plot the density distributions before and after normalization
# plotDensity(log2(exprs(ma.raw.data)), main="Before normalization (log2-transformed) - Probe level", xlim=c(-5,16)); grid()
# plotDensity(stdvec.condec.norm.data, main="RMA-normalized - Probeset level", xlim=c(-5,16)); grid()

#par(mfrow=c(1,1)) ## Restore single plot per page

## Compute a median for each row (probe)
#m <- rowMedians(stdvec.condec.norm.data)

## Centring: substract the median value of each probe
#rle <- sweep(stdvec.condec.norm.data,1, m, "-")

## Plot a box of the probe-wise centred values
# x11(width=8, height=8)
# boxplot(rle, pch=".", las=2, cex.axis=0.5)

###################################################
# function to conduct MCT test
MCT <- function(x,cl){
  n <- nrow(x)
  classes <- unique(cl)
  
  means.1 <- apply(x[,cl==classes[1]],1,mean,na.rm=TRUE)
  means.2 <- apply(x[,cl==classes[2]],1,mean,na.rm=TRUE)
  var.est.1 <- apply(x[,cl==classes[1]],1,var,na.rm=TRUE)
  var.est.2 <- apply(x[,cl==classes[2]],1,var,na.rm=TRUE)
  sd.est.1 <- sqrt(var.est.1)
  sd.est.2 <- sqrt(var.est.2)
  
  means.diff <- means.1 - means.2
  
  y <- !is.na(x)
  n.1 <- apply(y[,cl==classes[1]],1,sum,na.rm=TRUE)
  n.2 <- apply(y[,cl==classes[2]],1,sum,na.rm=TRUE)
  
  ## Calculate observed t value
  st.err.diff <- sqrt(var.est.1/n.1 + var.est.2/n.2)
  t.obs <- means.diff/st.err.diff
  
  h.lambda <- (var.est.1/n.1)/(var.est.1/n.1+var.est.2/n.2)
  
  
  p.value <- sapply(1:n, function(xx,h.lambda,t.obs){
    TNull <- rnorm(100000)/sqrt(h.lambda[xx]*rchisq(100000,n.1[xx]-1)/(n.1[xx]-1)+(1-h.lambda[xx])*rchisq(100000,n.2[xx]-1)/(n.2[xx]-1))
    mean(abs(TNull) >= abs(t.obs[xx]))
  }
  ,h.lambda,t.obs)
  
  list(dm=means.diff, p.value=p.value)
  
}

#stdvec.condec.norm.data <- set.expr.data

# analyze differential gene expression by iterating first over normalization 
# methods and then over differential gene expression methods


gene.diff.expr.method <- c( "MCT", "ttest")

set.comparison.gene.test <- as.vector( sapply( set.comparison.name, paste0, 
                                               ".", c( "fc", "pv", "fdr" ) ) )

for ( nmeth in norm.method )
{
  norm.data <- get( paste0( nmeth, ".norm.data" ) )
  
  stopifnot( rownames( norm.data ) == ma.probe )
  stopifnot( colnames( norm.data ) == set.sample )
  
  for ( gdemeth in gene.diff.expr.method )
  {
    if ( gdemeth == "MCT" )
    {
      gene.diff.expr <-  apply( set.comparison, 1, function( sc ) {
        tc.col <- set.sample.condition == sc[ 1 ] | 
          set.sample.condition == sc[ 2 ]
        tc.factor <- factor( ( 1*( set.sample.condition == sc[ 1 ] ) + 
                                 2*( set.sample.condition == sc[ 2 ] ) )[ tc.col ] )
        rtt <- MCT( norm.data[ , tc.col ], tc.factor )
        rtt$fdr <- p.adjust( rtt$p.value, method="fdr" )
        c( rtt$dm, rtt$p.value, rtt$fdr )
      } )
    }
    else if ( gdemeth == "ttest" )
    {
      gene.diff.expr <-  apply( set.comparison, 1, function( sc ) {
        tc.col <- set.sample.condition == sc[ 1 ] | 
          set.sample.condition == sc[ 2 ]
        tc.factor <- factor( ( 1*( set.sample.condition == sc[ 1 ] ) + 
                                 2*( set.sample.condition == sc[ 2 ] ) )[ tc.col ] )
        rtt=list()
        wlcht <- sapply(1:nrow(norm.data),function(z){ttp <- t.test(norm.data[z,]~tc.factor);c(ttp$estimate[1]-ttp$estimate[2],ttp$p.value)})
        rtt$dm <- wlcht[1,]
        rtt$p.value <- wlcht[2,]
        rtt$fdr <- p.adjust( rtt$p.value, method="fdr" )
        c(rtt$dm,rtt$p.value, rtt$fdr )
      } )
    }
    else
      stop( "wrong gene differential expression method ", gdemeth )
    
    dim( gene.diff.expr ) <- c( ma.probe.num, 
                                length( set.comparison.gene.test ) )
    rownames( gene.diff.expr ) <- ma.probe
    colnames( gene.diff.expr ) <- set.comparison.gene.test
    
    gene.diff.expr <- data.frame( probe=ma.probe, 
                                  gene.diff.expr[ ma.probe, ], stringsAsFactors=FALSE )
    
    assign( paste0( nmeth, ".norm.", gdemeth, ".gene.diff.expr" ), 
            gene.diff.expr )
  }
}

# plot roc curves

plot.norm.method <- c("median.condec")

plot.fdr.point <- c( 0.01, 0.05 )


for ( nmeth in norm.method )
{
  for ( gdemeth in gene.diff.expr.method )
  {
    
    
    gene.diff.expr <- get( paste0( nmeth, ".norm.", gdemeth, 
                                   ".gene.diff.expr" ) )
    
    gene.diff.expr.probe <- gene.diff.expr$probe
    gene.diff.expr.fdr <- gene.diff.expr$E2A.rearranged.EP.vs.BCR.ABL.pv
    hist(gene.diff.expr$E2A.rearranged.EP.vs.BCR.ABL.pv)
    gene.diff.expr.probe <- gene.diff.expr.probe[ 
      order( gene.diff.expr.fdr ) ]
    gene.diff.expr.fdr <- sort( gene.diff.expr.fdr )
    
    
    fdr.point <- sapply( plot.fdr.point, function( pfp ) 
      sum( gene.diff.expr.fdr < pfp ) )
    cat("gdemeth=",fdr.point,sep = "\n")
    assign( paste0( nmeth, ".norm.", gdemeth, ".gene.diff.expr",".pvalue" ), 
            gene.diff.expr$E2A.rearranged.EP.vs.BCR.ABL.pv )
    
  }
}

qqplot(log10(median.condec.norm.MCT.gene.diff.expr.pvalue),
       log10(median.condec.norm.ttest.gene.diff.expr.pvalue), xlab="MCT log10(p-value)", ylab="Welch's t log10(p-value)")
abline(0,1)

#################################################################

plot.norm.method <- c("median.condec")

plot.fdr.point <- c( 0.01 , 0.05)


for ( nmeth in norm.method )
{
  for ( gdemeth in gene.diff.expr.method )
  {
    
    
    gene.diff.expr <- get( paste0( nmeth, ".norm.", gdemeth, 
                                   ".gene.diff.expr" ) )
    
    gene.diff.expr.probe <- gene.diff.expr$probe
    gene.diff.expr.fdr <- gene.diff.expr$E2A.rearranged.EP.vs.BCR.ABL.fdr
    hist(gene.diff.expr$E2A.rearranged.EP.vs.BCR.ABL.fdr)
    gene.diff.expr.probe <- gene.diff.expr.probe[ 
      order( gene.diff.expr.fdr ) ]
    gene.diff.expr.fdr <- sort( gene.diff.expr.fdr )
    
    
    fdr.point <- sapply( plot.fdr.point, function( pfp ) 
      sum( gene.diff.expr.fdr < pfp ) )
    cat("gdemeth=",fdr.point,sep = "\n")
    assign( paste0( nmeth, ".norm.", gdemeth, ".gene.diff.expr",".fdr" ), 
            gene.diff.expr$E2A.rearranged.EP.vs.BCR.ABL.fdr )
    
  }
}

ma.annot[rownames(ma.annot)%in%rownames(ma.expr.data[median.condec.norm.MCT.gene.diff.expr.fdr<.01,]),]


qqplot(log10(median.condec.norm.MCT.gene.diff.expr.fdr),
       log10(median.condec.norm.ttest.gene.diff.expr.fdr), xlab="MCT log10(p-value)", ylab="Welch's t log10(p-value)")
abline(0,1)



####################################################################


postscript(file="PCAa.eps")
par(mar=c(8, 8, 4, 2))

X <- t(median.condec.norm.data)

svddec <- svd(X)
head(svddec$u)
expl <- svddec$d[1:2]^2/sum(svddec$d^2)*100

# PC <- prcomp(X,center=F,scale=F)
# head(scale(PC$x,center = F,scale=sqrt(PC$sdev^2*(12-1))))
# PC$sdev[1:2]^2/sum(PC$sdev^2)*100

plot(svddec$u[,1],svddec$u[,2],col=c(rep(1,8),rep(2,4)),pch=16,xlab="",ylab="",las=1, cex.axis=1.7)

title(main="(a)",cex.main=2)
title(xlab=paste0("PC-1 (",round(expl[1],2),"%)"), cex.lab=2,line = 4)
title(ylab=paste0("PC-2 (",round(expl[2],2),"%)"), cex.lab=2,line = 5)

dev.off()

postscript(file="PCAb.eps")
par(mar=c(8, 8, 4, 2))
#text(PC$x[,1],PC$x[,2], col=c(rep(1,8),rep(2,4)),labels = cancer.type)
significant.probesets.MCT <- median.condec.norm.MCT.gene.diff.expr.pvalue<0.01
significant.probesets.ttest <- median.condec.norm.ttest.gene.diff.expr.pvalue<0.01
XX <- t(median.condec.norm.data[which(significant.probesets.MCT!=significant.probesets.ttest),])
svddecsub <- svd(XX)
dim(svddecsub$u)
expl <- svddecsub$d[1:2]^2/sum(svddecsub$d^2)*100
plot(svddecsub$u[,1],svddecsub$u[,2],col=c(rep(1,8),rep(2,4)),pch=16,xlab="",ylab="",las=1, cex.axis=2)
title(main="(b)",cex.main=2)
title(xlab=paste0("PC-1 (",round(expl[1],2),"%)"), cex.lab=2,line = 4)
title(ylab=paste0("PC-2 (",round(expl[2],2),"%)"), cex.lab=2,line = 5)

dev.off()

####################################################################
library("mixOmics")
#################################################################################

X1 <- t(median.condec.norm.data)

nkeep <- c(6307,500,300,200,100)
spca1 <- spca(X1, center=TRUE,scale=TRUE, ncomp=2, keepX = c(nkeep[2],nkeep[1]))
spca2 <- spca(X1, center=TRUE,scale=TRUE, ncomp=2, keepX = c(nkeep[3],nkeep[1]))
spca3 <- spca(X1, center=TRUE,scale=TRUE, ncomp=2, keepX = c(nkeep[4],nkeep[1]))
spca4 <- spca(X1, center=TRUE,scale=TRUE, ncomp=2, keepX = c(nkeep[5],nkeep[1]))


postscript(file="SPCAa.eps")
par(mar=c(8, 8, 4, 2))
y1 <- scale(spca1$variates$X,scale=sqrt(apply(spca1$variates$X,2,var)*c(12-1)))
plot(y1[,1],y1[,2],xlab="", ylab="", pch=16, col=c(rep(1,8),rep(2,4)),cex=2, cex.axis=2.5,mgp = c(3, 1, 0),las=1)
title(main="(a)",cex.main=2)
title(xlab=paste0("PC-1 (",round(spca1$explained_variance[1]*100,2),"%)"), cex.lab=2,line = 5)
title(ylab=paste0("PC-2 (",round(spca1$explained_variance[2]*100,2),"%)"), cex.lab=2,line = 5)
dev.off()

postscript(file="SPCAb.eps")
par(mar=c(8, 8, 4, 2))
y2 <- scale(spca2$variates$X,scale=sqrt(apply(spca1$variates$X,2,var)*c(12-1)))
plot(y2[,1],y2[,2],xlab="", ylab="",pch=16, col=c(rep(1,8),rep(2,4)),cex=2, cex.axis=2.5,mgp = c(3, 1, 0),las=1)
title(main="(b)",cex.main=2)
title(xlab=paste0("PC-1 (",round(spca2$explained_variance[1]*100,2),"%)"), cex.lab=2,line = 5)
title(ylab=paste0("PC-2 (",round(spca2$explained_variance[2]*100,2),"%)"), cex.lab=2,line = 5)
dev.off()

postscript(file="SPCAc.eps")
par(mar=c(8, 8, 4, 2))
y3 <- scale(spca3$variates$X,scale=sqrt(apply(spca1$variates$X,2,var)*c(12-1)))
plot(y3[,1],y3[,2],xlab="", ylab="",pch=16, col=c(rep(1,8),rep(2,4)),cex=2, cex.axis=2.5,mgp = c(3, 1, 0),las=1)
title(main="(c)",cex.main=2)
title(xlab=paste0("PC-1 (",round(spca3$explained_variance[1]*100,2),"%)"), cex.lab=2,line = 5)
title(ylab=paste0("PC-2 (",round(spca3$explained_variance[2]*100,2),"%)"), cex.lab=2,line = 5)
dev.off()

postscript(file="SPCAd.eps")
par(mar=c(8, 8, 4, 2))
y4 <- scale(spca4$variates$X,scale=sqrt(apply(spca1$variates$X,2,var)*c(12-1)))
plot(y4[,1],y4[,2],xlab="", ylab="",pch=16, col=c(rep(1,8),rep(2,4)),cex=2, cex.axis=2.5,mgp = c(3, 1, 0),las=1)
title(main="(d)",cex.main=2)
title(xlab=paste0("PC-1 (",round(spca4$explained_variance[1]*100,2),"%)"), cex.lab=2,line = 5)
title(ylab=paste0("PC-2 (",round(spca4$explained_variance[2]*100,2),"%)"), cex.lab=2,line = 5)
dev.off()

library(xtable)

XXX <- t(median.condec.norm.data[which(significant.probesets.MCT!=significant.probesets.ttest),])
ctype <- as.factor(ma.sample.condition)
means <- aggregate(XXX, by=list(ctype), FUN=mean, na.rm=TRUE)
sds <- aggregate(XXX, by=list(ctype), FUN=sd, na.rm=TRUE)


# mns <- sapply(1:dim(XXX)[2],function(x){
#   y=XXX[,x]
#   z=data.frame(y,ctype) %>%  stats::setNames(c("y","ff")) 
#   ddply(z,~ff,summarise,means=mean(y), sds <- sd(y))
# },simplify = T)


spca1$explained_variance
spca2$explained_variance
spca3$explained_variance
spca4$explained_variance

table(significant.probesets.MCT) ## Count the number of significant probesets
table(significant.probesets.ttest)

MCT <- median.condec.norm.MCT.gene.diff.expr.pvalue[which(significant.probesets.MCT!=significant.probesets.ttest)]
welch <- median.condec.norm.ttest.gene.diff.expr.pvalue[which(significant.probesets.MCT!=significant.probesets.ttest)]

sv1 <- selectVar(spca1, comp = 1)
sv2 <- selectVar(spca2, comp = 1)
sv3 <- selectVar(spca3, comp = 1)
sv4 <- selectVar(spca4, comp = 1)

strars1 <- ifelse(colnames(means)[-1]%in%sv1$name,"*","")
strars2 <- ifelse(colnames(means)[-1]%in%sv2$name,"*","")
strars3 <- ifelse(colnames(means)[-1]%in%sv3$name,"*","")
strars4 <- ifelse(colnames(means)[-1]%in%sv4$name,"*","")

tab <- cbind(strars1,strars2,strars3,strars4,paste(round(means[1,-1],2)," (",round(sds[1,-1],2),")",sep=""),
             paste(round(means[2,-1],2)," (",round(sds[2,-1],2),")",sep=""),round(MCT,4),round(welch,4))
colnames(tab)=c("Probe Set ID","","","","mean (sd)", "mean (sd)", "MCT", "Welch")

gid <- unite(ma.annot[which(significant.probesets.MCT!=significant.probesets.ttest),],gid,sep="|",remove=F, c("SYMBOL","probe.id"))
rownames(tab) <- gid$gid

print(
  xtable(tab, label= "RRRRRRRRRRRR"),
  latex.environments=c("center"), 
  floating=FALSE, 
  include.rownames=TRUE
)

##################################################################################
#### q-values
significant.probesets.MCT.fdr <- median.condec.norm.MCT.gene.diff.expr.fdr<.01

XXX <- t(median.condec.norm.data[which(significant.probesets.MCT.fdr),])
ctype <- as.factor(ma.sample.condition)
means <- aggregate(XXX, by=list(ctype), FUN=mean, na.rm=TRUE)
sds <- aggregate(XXX, by=list(ctype), FUN=sd, na.rm=TRUE)

MCT <- median.condec.norm.MCT.gene.diff.expr.fdr[which(significant.probesets.MCT.fdr)]
welch <- median.condec.norm.ttest.gene.diff.expr.fdr[which(significant.probesets.MCT.fdr)]

strars1 <- ifelse(colnames(means)[-1]%in%sv1$name,"*","")
strars2 <- ifelse(colnames(means)[-1]%in%sv2$name,"*","")
strars3 <- ifelse(colnames(means)[-1]%in%sv3$name,"*","")
strars4 <- ifelse(colnames(means)[-1]%in%sv4$name,"*","")


tab <- cbind(strars1,strars2,strars3,strars4,paste(round(means[1,-1],2)," (",round(sds[1,-1],2),")",sep=""),
             paste(round(means[2,-1],2)," (",round(sds[2,-1],2),")",sep=""),round(MCT,4),round(welch,4))
colnames(tab)=c("Probe Set ID","","","","mean (sd)", "mean (sd)", "MCT", "Welch")

gid <- unite(ma.annot[which(significant.probesets.MCT.fdr),],gid,sep="|",remove=F, c("SYMBOL","probe.id"))
rownames(tab) <- gid$gid

print(
  xtable(tab, label= "RRRRRRRRRRRR"),
  latex.environments=c("center"), 
  floating=FALSE, 
  include.rownames=TRUE
)

#################################################################

nkeep <- c(6307,500,300,200,100)
X2 <- t(median.condec.norm.data)
plsda1 <- splsda(X2, ctype, scale=TRUE, ncomp=2, keepX = c(nkeep[2],nkeep[1]))
plsda2 <- splsda(X2, ctype, scale=TRUE, ncomp=2, keepX = c(nkeep[3],nkeep[1]))
plsda3 <- splsda(X2, ctype, scale=TRUE, ncomp=2, keepX = c(nkeep[4],nkeep[1]))
plsda4 <- splsda(X2, ctype, scale=TRUE, ncomp=2, keepX = c(nkeep[5],nkeep[1]))


plotIndiv(plsda1, comp = c(1,2),
          group = ctype, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE,
          title = 'SRBCT, sPLSDA comp 1 - 2')


postscript(file="SPLSDAa.eps")
par(mar=c(8, 8, 4, 2))
y1 <- scale(plsda1$variates$X,scale=sqrt(apply(plsda1$variates$X,2,var)*c(12-1)))
plot(y1[,1],y1[,2],xlab="", ylab="", pch=16, col=c(rep(1,8),rep(2,4)),cex=2, cex.axis=2.5,mgp = c(3, 1, 0),las=1)
title(main="(a)",cex.main=2)
title(xlab=paste0("component-1 (",round(plsda1$explained_variance$X[1]*100,2),"%)"), cex.lab=2,line = 5)
title(ylab=paste0("component-2 (",round(plsda1$explained_variance$X[2]*100,2),"%)"), cex.lab=2,line = 5)
dev.off()

postscript(file="SPLSDAb.eps")
par(mar=c(8, 8, 4, 2))
y2 <- scale(plsda2$variates$X,scale=sqrt(apply(plsda2$variates$X,2,var)*c(12-1)))
plot(y1[,1],y1[,2],xlab="", ylab="", pch=16, col=c(rep(1,8),rep(2,4)),cex=2, cex.axis=2.5,mgp = c(3, 1, 0),las=1)
title(main="(b)",cex.main=2)
title(xlab=paste0("component-1 (",round(plsda2$explained_variance$X[1]*100,2),"%)"), cex.lab=2,line = 5)
title(ylab=paste0("component-2 (",round(plsda2$explained_variance$X[2]*100,2),"%)"), cex.lab=2,line = 5)
dev.off()


postscript(file="SPLSDAc.eps")
par(mar=c(8, 8, 4, 2))
y3 <- scale(plsda3$variates$X,scale=sqrt(apply(plsda3$variates$X,2,var)*c(12-1)))
plot(y1[,1],y1[,2],xlab="", ylab="", pch=16, col=c(rep(1,8),rep(2,4)),cex=2, cex.axis=2.5,mgp = c(3, 1, 0),las=1)
title(main="(c)",cex.main=2)
title(xlab=paste0("component-1 (",round(plsda3$explained_variance$X[1]*100,2),"%)"), cex.lab=2,line = 5)
title(ylab=paste0("component-2 (",round(plsda3$explained_variance$X[2]*100,2),"%)"), cex.lab=2,line = 5)
dev.off()

postscript(file="SPLSDAd.eps")
par(mar=c(8, 8, 4, 2))
y4 <- scale(plsda4$variates$X,scale=sqrt(apply(plsda4$variates$X,2,var)*c(12-1)))
plot(y1[,1],y1[,2],xlab="", ylab="", pch=16, col=c(rep(1,8),rep(2,4)),cex=2, cex.axis=2.5,mgp = c(3, 1, 0),las=1)
title(main="(d)",cex.main=2)
title(xlab=paste0("component-1 (",round(plsda4$explained_variance$X[1]*100,2),"%)"), cex.lab=2,line = 5)
title(ylab=paste0("component-2 (",round(plsda4$explained_variance$X[2]*100,2),"%)"), cex.lab=2,line = 5)
dev.off()

plsda1$explained_variance
plsda2$explained_variance
plsda3$explained_variance
plsda4$explained_variance

XXX <- t(median.condec.norm.data[which(significant.probesets.MCT!=significant.probesets.ttest),])
ctype <- as.factor(ma.sample.condition)
means <- aggregate(XXX, by=list(ctype), FUN=mean, na.rm=TRUE)
sds <- aggregate(XXX, by=list(ctype), FUN=sd, na.rm=TRUE)


MCT <- median.condec.norm.MCT.gene.diff.expr.pvalue[which(significant.probesets.MCT!=significant.probesets.ttest)]
welch <- median.condec.norm.ttest.gene.diff.expr.pvalue[which(significant.probesets.MCT!=significant.probesets.ttest)]

sv1 <- selectVar(plsda1, comp = 1)
sv2 <- selectVar(plsda2, comp = 1)
sv3 <- selectVar(plsda3, comp = 1)
sv4 <- selectVar(plsda4, comp = 1)

strars1 <- ifelse(colnames(means)[-1]%in%sv1$name,"*","")
strars2 <- ifelse(colnames(means)[-1]%in%sv2$name,"*","")
strars3 <- ifelse(colnames(means)[-1]%in%sv3$name,"*","")
strars4 <- ifelse(colnames(means)[-1]%in%sv4$name,"*","")

tab <- cbind(strars1,strars2,strars3,strars4,paste(round(means[1,-1],2)," (",round(sds[1,-1],2),")",sep=""),
             paste(round(means[2,-1],2)," (",round(sds[2,-1],2),")",sep=""),round(MCT,4),round(welch,4))
colnames(tab)=c("Probe Set ID","","","","mean (sd)", "mean (sd)", "MCT", "Welch")

gid <- unite(ma.annot[which(significant.probesets.MCT!=significant.probesets.ttest),],gid,sep="|",remove=F, c("SYMBOL","probe.id"))
rownames(tab) <- gid$gid

print(
  xtable(tab, label= "RRRRRRRRRRRR"),
  latex.environments=c("center"), 
  floating=FALSE, 
  include.rownames=TRUE
)

##################################################################################
#### q-values
significant.probesets.MCT.fdr <- median.condec.norm.MCT.gene.diff.expr.fdr<.01

XXX <- t(median.condec.norm.data[which(significant.probesets.MCT.fdr),])
ctype <- as.factor(ma.sample.condition)
means <- aggregate(XXX, by=list(ctype), FUN=mean, na.rm=TRUE)
sds <- aggregate(XXX, by=list(ctype), FUN=sd, na.rm=TRUE)

MCT <- median.condec.norm.MCT.gene.diff.expr.fdr[which(significant.probesets.MCT.fdr)]
welch <- median.condec.norm.ttest.gene.diff.expr.fdr[which(significant.probesets.MCT.fdr)]

strars1 <- ifelse(colnames(means)[-1]%in%sv1$name,"*","")
strars2 <- ifelse(colnames(means)[-1]%in%sv2$name,"*","")
strars3 <- ifelse(colnames(means)[-1]%in%sv3$name,"*","")
strars4 <- ifelse(colnames(means)[-1]%in%sv4$name,"*","")


tab <- cbind(strars1,strars2,strars3,strars4,paste(round(means[1,-1],2)," (",round(sds[1,-1],2),")",sep=""),
             paste(round(means[2,-1],2)," (",round(sds[2,-1],2),")",sep=""),round(MCT,4),round(welch,4))
colnames(tab)=c("Probe Set ID","","","","mean (sd)", "mean (sd)", "MCT", "Welch")

gid <- unite(ma.annot[which(significant.probesets.MCT.fdr),],gid,sep="|",remove=F, c("SYMBOL","probe.id"))
rownames(tab) <- gid$gid

print(
  xtable(tab, label= "RRRRRRRRRRRR"),
  latex.environments=c("center"), 
  floating=FALSE, 
  include.rownames=TRUE
)

#################################################################
####################################################################
library("mixOmics")
#################################################################################

X1 <- t(median.condec.norm.data)

nkeep <- c(6307,500,300,200,100)
spca1 <- spca(X1, center=TRUE,scale=TRUE, ncomp=2, keepX = c(nkeep[2],nkeep[1]))
spca2 <- spca(X1, center=TRUE,scale=TRUE, ncomp=2, keepX = c(nkeep[3],nkeep[1]))
spca3 <- spca(X1, center=TRUE,scale=TRUE, ncomp=2, keepX = c(nkeep[4],nkeep[1]))
spca4 <- spca(X1, center=TRUE,scale=TRUE, ncomp=2, keepX = c(nkeep[5],nkeep[1]))


postscript(file="SPCAa.eps")
par(mar=c(8, 8, 4, 2))
y1 <- scale(spca1$variates$X,scale=sqrt(apply(spca1$variates$X,2,var)*c(12-1)))
plot(y1[,1],y1[,2],xlab="", ylab="", pch=16, col=c(rep(1,8),rep(2,4)),cex=2, cex.axis=2.5,mgp = c(3, 1, 0),las=1)
title(main="(a)",cex.main=2)
title(xlab=paste0("PC-1 (",round(spca1$explained_variance[1]*100,2),"%)"), cex.lab=2,line = 5)
title(ylab=paste0("PC-2 (",round(spca1$explained_variance[2]*100,2),"%)"), cex.lab=2,line = 5)
dev.off()

postscript(file="SPCAb.eps")
par(mar=c(8, 8, 4, 2))
y2 <- scale(spca2$variates$X,scale=sqrt(apply(spca1$variates$X,2,var)*c(12-1)))
plot(y2[,1],y2[,2],xlab="", ylab="",pch=16, col=c(rep(1,8),rep(2,4)),cex=2, cex.axis=2.5,mgp = c(3, 1, 0),las=1)
title(main="(b)",cex.main=2)
title(xlab=paste0("PC-1 (",round(spca2$explained_variance[1]*100,2),"%)"), cex.lab=2,line = 5)
title(ylab=paste0("PC-2 (",round(spca2$explained_variance[2]*100,2),"%)"), cex.lab=2,line = 5)
dev.off()

postscript(file="SPCAc.eps")
par(mar=c(8, 8, 4, 2))
y3 <- scale(spca3$variates$X,scale=sqrt(apply(spca1$variates$X,2,var)*c(12-1)))
plot(y3[,1],y3[,2],xlab="", ylab="",pch=16, col=c(rep(1,8),rep(2,4)),cex=2, cex.axis=2.5,mgp = c(3, 1, 0),las=1)
title(main="(c)",cex.main=2)
title(xlab=paste0("PC-1 (",round(spca3$explained_variance[1]*100,2),"%)"), cex.lab=2,line = 5)
title(ylab=paste0("PC-2 (",round(spca3$explained_variance[2]*100,2),"%)"), cex.lab=2,line = 5)
dev.off()

postscript(file="SPCAd.eps")
par(mar=c(8, 8, 4, 2))
y4 <- scale(spca4$variates$X,scale=sqrt(apply(spca1$variates$X,2,var)*c(12-1)))
plot(y4[,1],y4[,2],xlab="", ylab="",pch=16, col=c(rep(1,8),rep(2,4)),cex=2, cex.axis=2.5,mgp = c(3, 1, 0),las=1)
title(main="(d)",cex.main=2)
title(xlab=paste0("PC-1 (",round(spca4$explained_variance[1]*100,2),"%)"), cex.lab=2,line = 5)
title(ylab=paste0("PC-2 (",round(spca4$explained_variance[2]*100,2),"%)"), cex.lab=2,line = 5)
dev.off()


XXX <- t(median.condec.norm.data[which(significant.probesets.MCT!=significant.probesets.ttest),])
ctype <- as.factor(ma.sample.condition)
means <- aggregate(XXX, by=list(ctype), FUN=mean, na.rm=TRUE)
sds <- aggregate(XXX, by=list(ctype), FUN=sd, na.rm=TRUE)


# mns <- sapply(1:dim(XXX)[2],function(x){
#   y=XXX[,x]
#   z=data.frame(y,ctype) %>%  stats::setNames(c("y","ff")) 
#   ddply(z,~ff,summarise,means=mean(y), sds <- sd(y))
# },simplify = T)


spca1$explained_variance
spca2$explained_variance
spca3$explained_variance
spca4$explained_variance

table(significant.probesets.MCT) ## Count the number of significant probesets
table(significant.probesets.ttest)

MCT <- median.condec.norm.MCT.gene.diff.expr.pvalue[which(significant.probesets.MCT!=significant.probesets.ttest)]
welch <- median.condec.norm.ttest.gene.diff.expr.pvalue[which(significant.probesets.MCT!=significant.probesets.ttest)]

sv1 <- selectVar(spca1, comp = 1)
sv2 <- selectVar(spca2, comp = 1)
sv3 <- selectVar(spca3, comp = 1)
sv4 <- selectVar(spca4, comp = 1)

strars1 <- ifelse(colnames(means)[-1]%in%sv1$name,"*","")
strars2 <- ifelse(colnames(means)[-1]%in%sv2$name,"*","")
strars3 <- ifelse(colnames(means)[-1]%in%sv3$name,"*","")
strars4 <- ifelse(colnames(means)[-1]%in%sv4$name,"*","")

tab <- cbind(paste(round(means[1,-1],2)," (",round(sds[1,-1],2),")",sep=""),
             paste(round(means[2,-1],2)," (",round(sds[2,-1],2),")",sep=""),round(MCT,4),round(welch,4))
colnames(tab)=c("mean (sd)", "mean (sd)", "MCT", "Welch")

gid <- unite(ma.annot[which(significant.probesets.MCT!=significant.probesets.ttest),],gid,sep="|",remove=F, c("SYMBOL","probe.id"))
rownames(tab) <- gid$gid

print(
  xtable(tab, label= "RRRRRRRRRRRR"),
  latex.environments=c("center"), 
  floating=FALSE, 
  include.rownames=TRUE
)

##################################################################################
#### q-values
significant.probesets.MCT.fdr <- median.condec.norm.MCT.gene.diff.expr.fdr<.01

XXX <- t(median.condec.norm.data[which(significant.probesets.MCT.fdr),])
ctype <- as.factor(ma.sample.condition)
means <- aggregate(XXX, by=list(ctype), FUN=mean, na.rm=TRUE)
sds <- aggregate(XXX, by=list(ctype), FUN=sd, na.rm=TRUE)

MCT <- median.condec.norm.MCT.gene.diff.expr.fdr[which(significant.probesets.MCT.fdr)]
welch <- median.condec.norm.ttest.gene.diff.expr.fdr[which(significant.probesets.MCT.fdr)]

strars1 <- ifelse(colnames(means)[-1]%in%sv1$name,"*","")
strars2 <- ifelse(colnames(means)[-1]%in%sv2$name,"*","")
strars3 <- ifelse(colnames(means)[-1]%in%sv3$name,"*","")
strars4 <- ifelse(colnames(means)[-1]%in%sv4$name,"*","")


tab <- cbind(strars1,strars2,strars3,strars4,paste(round(means[1,-1],2)," (",round(sds[1,-1],2),")",sep=""),
             paste(round(means[2,-1],2)," (",round(sds[2,-1],2),")",sep=""),round(MCT,4),round(welch,4))
colnames(tab)=c("Probe Set ID","","","","mean (sd)", "mean (sd)", "MCT", "Welch")

gid <- unite(ma.annot[which(significant.probesets.MCT.fdr),],gid,sep="|",remove=F, c("SYMBOL","probe.id"))
rownames(tab) <- gid$gid

print(
  xtable(tab, label= "RRRRRRRRRRRR"),
  latex.environments=c("center"), 
  floating=FALSE, 
  include.rownames=TRUE
)

#################################################################

nkeep <- c(6307,500,300,200,100)
X2 <- t(median.condec.norm.data)
plsda1 <- splsda(X2, ctype, scale=TRUE, ncomp=2, keepX = c(nkeep[2],nkeep[1]))
plsda2 <- splsda(X2, ctype, scale=TRUE, ncomp=2, keepX = c(nkeep[3],nkeep[1]))
plsda3 <- splsda(X2, ctype, scale=TRUE, ncomp=2, keepX = c(nkeep[4],nkeep[1]))
plsda4 <- splsda(X2, ctype, scale=TRUE, ncomp=2, keepX = c(nkeep[5],nkeep[1]))


plotIndiv(plsda1, comp = c(1,2),
          group = ctype, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE,
          title = 'SRBCT, sPLSDA comp 1 - 2')


postscript(file="SPLSDAa.eps")
par(mar=c(8, 8, 4, 2))
y1 <- scale(plsda1$variates$X,scale=sqrt(apply(plsda1$variates$X,2,var)*c(12-1)))
plot(y1[,1],y1[,2],xlab="", ylab="", pch=16, col=c(rep(1,8),rep(2,4)),cex=2, cex.axis=2.5,mgp = c(3, 1, 0),las=1)
title(main="(a)",cex.main=2)
title(xlab=paste0("component-1 (",round(plsda1$explained_variance$X[1]*100,2),"%)"), cex.lab=2,line = 5)
title(ylab=paste0("component-2 (",round(plsda1$explained_variance$X[2]*100,2),"%)"), cex.lab=2,line = 5)
dev.off()

postscript(file="SPLSDAb.eps")
par(mar=c(8, 8, 4, 2))
y2 <- scale(plsda2$variates$X,scale=sqrt(apply(plsda2$variates$X,2,var)*c(12-1)))
plot(y1[,1],y1[,2],xlab="", ylab="", pch=16, col=c(rep(1,8),rep(2,4)),cex=2, cex.axis=2.5,mgp = c(3, 1, 0),las=1)
title(main="(b)",cex.main=2)
title(xlab=paste0("component-1 (",round(plsda2$explained_variance$X[1]*100,2),"%)"), cex.lab=2,line = 5)
title(ylab=paste0("component-2 (",round(plsda2$explained_variance$X[2]*100,2),"%)"), cex.lab=2,line = 5)
dev.off()


postscript(file="SPLSDAc.eps")
par(mar=c(8, 8, 4, 2))
y3 <- scale(plsda3$variates$X,scale=sqrt(apply(plsda3$variates$X,2,var)*c(12-1)))
plot(y1[,1],y1[,2],xlab="", ylab="", pch=16, col=c(rep(1,8),rep(2,4)),cex=2, cex.axis=2.5,mgp = c(3, 1, 0),las=1)
title(main="(c)",cex.main=2)
title(xlab=paste0("component-1 (",round(plsda3$explained_variance$X[1]*100,2),"%)"), cex.lab=2,line = 5)
title(ylab=paste0("component-2 (",round(plsda3$explained_variance$X[2]*100,2),"%)"), cex.lab=2,line = 5)
dev.off()

postscript(file="SPLSDAd.eps")
par(mar=c(8, 8, 4, 2))
y4 <- scale(plsda4$variates$X,scale=sqrt(apply(plsda4$variates$X,2,var)*c(12-1)))
plot(y1[,1],y1[,2],xlab="", ylab="", pch=16, col=c(rep(1,8),rep(2,4)),cex=2, cex.axis=2.5,mgp = c(3, 1, 0),las=1)
title(main="(d)",cex.main=2)
title(xlab=paste0("component-1 (",round(plsda4$explained_variance$X[1]*100,2),"%)"), cex.lab=2,line = 5)
title(ylab=paste0("component-2 (",round(plsda4$explained_variance$X[2]*100,2),"%)"), cex.lab=2,line = 5)
dev.off()

plsda1$explained_variance
plsda2$explained_variance
plsda3$explained_variance
plsda4$explained_variance

XXX <- t(median.condec.norm.data[which(significant.probesets.MCT!=significant.probesets.ttest),])
ctype <- as.factor(ma.sample.condition)
means <- aggregate(XXX, by=list(ctype), FUN=mean, na.rm=TRUE)
sds <- aggregate(XXX, by=list(ctype), FUN=sd, na.rm=TRUE)


MCT <- median.condec.norm.MCT.gene.diff.expr.pvalue[which(significant.probesets.MCT!=significant.probesets.ttest)]
welch <- median.condec.norm.ttest.gene.diff.expr.pvalue[which(significant.probesets.MCT!=significant.probesets.ttest)]

sv1 <- selectVar(plsda1, comp = 1)
sv2 <- selectVar(plsda2, comp = 1)
sv3 <- selectVar(plsda3, comp = 1)
sv4 <- selectVar(plsda4, comp = 1)

strars1 <- ifelse(colnames(means)[-1]%in%sv1$name,"*","")
strars2 <- ifelse(colnames(means)[-1]%in%sv2$name,"*","")
strars3 <- ifelse(colnames(means)[-1]%in%sv3$name,"*","")
strars4 <- ifelse(colnames(means)[-1]%in%sv4$name,"*","")

tab <- cbind(paste(round(means[1,-1],2)," (",round(sds[1,-1],2),")",sep=""),
             paste(round(means[2,-1],2)," (",round(sds[2,-1],2),")",sep=""),round(MCT,4),round(welch,4))
colnames(tab)=c("mean (sd)", "mean (sd)", "MCT", "Welch")

gid <- unite(ma.annot[which(significant.probesets.MCT!=significant.probesets.ttest),],gid,sep="|",remove=F, c("SYMBOL","probe.id"))
rownames(tab) <- gid$gid

print(
  xtable(tab, label= "RRRRRRRRRRRR"),
  latex.environments=c("center"), 
  floating=FALSE, 
  include.rownames=TRUE
)

##################################################################################
#### q-values
significant.probesets.MCT.fdr <- median.condec.norm.MCT.gene.diff.expr.fdr<.01

XXX <- t(median.condec.norm.data[which(significant.probesets.MCT.fdr),])
ctype <- as.factor(ma.sample.condition)
means <- aggregate(XXX, by=list(ctype), FUN=mean, na.rm=TRUE)
sds <- aggregate(XXX, by=list(ctype), FUN=sd, na.rm=TRUE)

MCT <- median.condec.norm.MCT.gene.diff.expr.fdr[which(significant.probesets.MCT.fdr)]
welch <- median.condec.norm.ttest.gene.diff.expr.fdr[which(significant.probesets.MCT.fdr)]

strars1 <- ifelse(colnames(means)[-1]%in%sv1$name,"*","")
strars2 <- ifelse(colnames(means)[-1]%in%sv2$name,"*","")
strars3 <- ifelse(colnames(means)[-1]%in%sv3$name,"*","")
strars4 <- ifelse(colnames(means)[-1]%in%sv4$name,"*","")


tab <- cbind(paste(round(means[1,-1],2)," (",round(sds[1,-1],2),")",sep=""),
             paste(round(means[2,-1],2)," (",round(sds[2,-1],2),")",sep=""),round(MCT,4),round(welch,4))
colnames(tab)=c("mean (sd)", "mean (sd)", "MCT", "Welch")

gid <- unite(ma.annot[which(significant.probesets.MCT.fdr),],gid,sep="|",remove=F, c("SYMBOL","probe.id"))
rownames(tab) <- gid$gid

print(
  xtable(tab, label= "RRRRRRRRRRRR"),
  latex.environments=c("center"), 
  floating=FALSE, 
  include.rownames=TRUE
)

#################################################################