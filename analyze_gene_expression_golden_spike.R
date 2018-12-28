# Most of these codes are copied from https://github.com/carlosproca/gene-expr-norm-paper
# Download the raw microarray data of the Golden Spike dataset, 
# and copy the *.CEL files into the rawdata subdirectory. 
# The dataset is available at the online version of the paper by 
# Choe et al., Genome Biology 6, R16, 2005, in the additional data files 6 and 7.


setwd("C:/Golden-Spike-data-analysis")
source( "normalize_median_condec.r" )
source( "normalize_stdvec_condec.r" )
source( "select_nodup.r" )
source( "watson_u2.r" )

# set directories

ma.data.dir <- "C:/Golden-Spike-data-analysis/rawdata"   # raw data directory

stopifnot( ma.data.dir != "" )

#gene.diff.expr.dir <- "gene_diff_expr_golden_spike"


# define experimental conditions and samples

ma.treatment <- c( "spike" )
ma.control <- c( "constant" )
ma.condition <- c( ma.treatment, unique( ma.control ) )

ma.condition.sample.num <- c( 3, 3 )
ma.sample.condition <- rep( ma.condition, ma.condition.sample.num )
ma.sample.postfix <- unlist( lapply( ma.condition.sample.num, 
                                     function( mcsn ) 1:mcsn ) )
ma.sample <- paste0( ma.sample.condition, '.', ma.sample.postfix )

ma.sample.num <- length( ma.sample )
ma.condition.num <- length( ma.condition )
ma.treatment.num <- length( ma.treatment )


# read platform annotation and identify probes
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("drosgenome1.db", version = "3.8")

ma.species <- "drosgenome1"
ma.db.name <- paste0( ma.species, ".db" )

library( ma.db.name, character.only=TRUE )
ma.db <- get( ma.db.name )

ma.probe.platform <- sort( keys( ma.db, "PROBEID" ) )

ma.annot <- select.nodup( ma.db, ma.probe.platform, 
                          c( "SYMBOL", "GENENAME", "FLYBASE" ), "PROBEID" )

names( ma.annot ) <- sub( "PROBEID", "probe.id", names( ma.annot ), fixed=TRUE )


# read probe set fold change from experiment design

exp.design.file.name <- "golden_spike_probe_set_fold_change.csv"

exp.design <- read.csv( exp.design.file.name, stringsAsFactors=FALSE, 
                        comment.char="#" )

names( exp.design ) <- sub( "probe.set.id", "probe", names( exp.design ), 
                            fixed=TRUE )
names( exp.design ) <- sub( "fold.change", "fc", names( exp.design ), 
                            fixed=TRUE )

stopifnot( exp.design$probe == ma.probe.platform )


# read raw expression data

ma.data.def.file.name <- "golden_spike_data_def.txt"

ma.data.def <- read.delim( ma.data.def.file.name, stringsAsFactors=FALSE )

ma.sample.file.name.read <- ma.data.def$FileName
ma.sample.read <- ma.data.def$Sample
ma.sample.condition.read <- ma.data.def$Condition

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("affy", version = "3.8")

require(affy)
ma.raw.data <- ReadAffy( filenames = 
                           file.path( ma.data.dir, ma.data.def$FileName ) )

# verify probes of raw data

stopifnot( sort( featureNames( ma.raw.data ) ) == ma.probe.platform )


# verify samples of raw data

stopifnot( sampleNames( ma.raw.data ) == ma.sample.file.name.read )

stopifnot( sort( ma.sample ) == sort( ma.sample.read ) )

stopifnot( ma.sample.condition[ order( ma.sample ) ] == 
             ma.sample.condition.read[ order( ma.sample.read ) ] )


# get background-corrected data

ma.expr.data <- log2( exprs( mas5( ma.raw.data, normalize=FALSE ) ) )

# remove expression data below noise threshold

ma.mas5.call <- exprs( mas5calls( ma.raw.data ) )
ma.expr.data.missing.bol <- ma.mas5.call != "P"

stopifnot( rownames( ma.expr.data ) == rownames( ma.expr.data.missing.bol ) )
stopifnot( colnames( ma.expr.data ) == colnames( ma.expr.data.missing.bol ) )

ma.expr.data[ ma.expr.data.missing.bol ] <- NA

present <- apply(ma.mas5.call, 1, function(x)(sum(x == "P")))
table(present)

# remove control probes and probes with missing data
# reorder probes in expression data

ma.probe.control.bol <- grepl( "^AFFX-", ma.probe.platform )
ma.probe.missing.data.bol <- rowSums( is.na( ma.expr.data ) ) > 0

ma.probe <- ma.probe.platform[ ! ma.probe.control.bol &
                                 ! ma.probe.missing.data.bol ]
ma.probe.num <- length( ma.probe )

ma.expr.data <- ma.expr.data[ ma.probe, ]


# identify probes from experiment design

ma.probe.known.pos <- intersect( ma.probe, 
                                 exp.design$probe[ exp.design$fc > 1 ] )
ma.probe.known.neg <- intersect( ma.probe, 
                                 exp.design$probe[ exp.design$fc == 1 ] )
ma.probe.known <- intersect( ma.probe, 
                             exp.design$probe[ exp.design$fc >= 1 ] )
ma.probe.unknown <- intersect( ma.probe, 
                               exp.design$probe[ exp.design$fc < 1 ] )


# reorder and rename samples in expression data

ma.expr.data <- ma.expr.data[ , match( ma.sample, ma.sample.read ) ]

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

random.seed <- 50000
set.seed( random.seed )

norm.method <- c( "stdvec.condec","median.condec")

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


#########################################################
x11(width=8, height=10) ## Open a graphical window with specific dimensions
par(mfrow=c(2,1)) ## Share this window between two plots
boxplot(log2(exprs(ma.raw.data)), pch=".", las=2, cex.axis=0.5, main="Before normalization (log2-transformed) - Probe level")
grid()
boxplot(stdvec.condec.norm.data, pch=".", las=2, cex.axis=0.5, main="RMA-normalized - Probeset level")
grid()

## Plot the density distributions before and after normalization
plotDensity(log2(exprs(ma.raw.data)), main="Before normalization (log2-transformed) - Probe level", xlim=c(0,16)); grid()
plotDensity(stdvec.condec.norm.data, main="RMA-normalized - Probeset level", xlim=c(0,16)); grid()

par(mfrow=c(1,1)) ## Restore single plot per page

## Compute a median for each row (probe)
m <- rowMedians(stdvec.condec.norm.data)

## Centring: substract the median value of each probe
rle <- sweep(stdvec.condec.norm.data,1, m, "-")

## Plot a box of the probe-wise centred values
x11(width=8, height=8)
boxplot(rle, pch=".", las=2, cex.axis=0.5)
###############################################################
x11(width=8, height=10) ## Open a graphical window with specific dimensions
par(mfrow=c(2,1)) ## Share this window between two plots
boxplot(log2(exprs(ma.raw.data)), pch=".", las=2, cex.axis=0.5, main="Before normalization (log2-transformed) - Probe level")
grid()
boxplot(median.condec.norm.data, pch=".", las=2, cex.axis=0.5, main="RMA-normalized - Probeset level")
grid()

## Plot the density distributions before and after normalization
plotDensity(log2(exprs(ma.raw.data)), main="Before normalization (log2-transformed) - Probe level", xlim=c(0,16)); grid()
plotDensity(median.condec.norm.data, main="RMA-normalized - Probeset level", xlim=c(0,16)); grid()

par(mfrow=c(1,1)) ## Restore single plot per page

## Compute a median for each row (probe)
m <- rowMedians(median.condec.norm.data)

## Centring: substract the median value of each probe
rle <- sweep(stdvec.condec.norm.data,1, m, "-")

## Plot a box of the probe-wise centred values
x11(width=8, height=8)
boxplot(rle, pch=".", las=2, cex.axis=0.5)
####################################################

# build model matrix for analysis of differential gene expression with limma
#stdvec.condec.norm.data <- set.expr.data
condition.factor <- factor( set.sample.condition, levels=set.condition )

model.design <- model.matrix( ~0 + condition.factor )
rownames( model.design ) <- set.sample
colnames( model.design ) <- set.condition

diff.expr.contrast <- apply( set.comparison, 1, function( sc ) {
  de.cont <- rep( 0, set.condition.num )
  de.cont[ set.condition == sc[ 1 ] ] <- 1
  de.cont[ set.condition == sc[ 2 ] ] <- -1
  de.cont
} )

rownames( diff.expr.contrast ) <- set.condition
colnames( diff.expr.contrast ) <- set.comparison.name

# function to conduct MCT test
MCT <- function(x,cl){
  n <- nrow(x)
  p <- ncol(x)
  classes <- unique(cl)
  
  means.1 <- apply(x[,cl==classes[1]],1,mean)
  means.2 <- apply(x[,cl==classes[2]],1,mean)
  var.est.1 <- apply(x[,cl==classes[1]],1,var)
  var.est.2 <- apply(x[,cl==classes[2]],1,var)
  sd.est.1 <- sqrt(var.est.1)
  sd.est.2 <- sqrt(var.est.2)
  
  means.diff <- means.1 - means.2
  
  n.1 <- sum(cl == classes[1])
  n.2 <- sum(cl == classes[2])
  
  ## Calculate observed t value
  st.err.diff <- sqrt(var.est.1/n.1 + var.est.2/n.2)
  t.obs <- means.diff/st.err.diff
  
  h.lambda <- (var.est.1/n.1)/(var.est.1/n.1+var.est.2/n.2)
  
  
  p.value <- sapply(1:n, function(xx,h.lambda,t.obs){
    TNull <- rnorm(100000)/sqrt(h.lambda[xx]*rchisq(100000,n.1-1)/(n.1-1)+(1-h.lambda[xx])*rchisq(100000,n.2-1)/(n.2-1))
    mean(abs(TNull) >= abs(t.obs[xx]))
  }
  ,h.lambda,t.obs)
  
  list(dm=means.diff, p.value=p.value)
  
}

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
gene.diff.expr.method <- c( "MCT", "ttest")

norm.method.pch <- c( 16, 17)
norm.method.lty <- c( 1, 3)
norm.method.col <- c( "black", "red")

norm.method.bg.col <- c("black", "red")

names( norm.method.pch ) <- gene.diff.expr.method
names( norm.method.lty ) <- gene.diff.expr.method
names( norm.method.col ) <- gene.diff.expr.method
names( norm.method.bg.col ) <- gene.diff.expr.method

plot.norm.method <- c( "stdvec.condec","median.condec")

plot.fdr.point <- c( 0.01)


par( mar = c( 1.90, 1.70, 0.15, 0.15 ), mgp = c( 0.90, 0, 0 ), tcl=-0.10 )

plot( 1, type="n", cex.lab=0.67, cex.axis=0.67, 
      xlim = c( 0, .1 ), ylim = c( 0, 1 ), 
      xlab="False positive rate", 
      ylab="True positive rate" )

for ( gdemeth in gene.diff.expr.method )
{
  for ( nmeth in plot.norm.method )
  {
    gene.diff.expr <- get( paste0( nmeth, ".norm.", gdemeth, 
                                   ".gene.diff.expr" ) )
    
    gene.diff.expr.probe <- gene.diff.expr$probe
    gene.diff.expr.fdr <- gene.diff.expr$spike.vs.constant.fdr
    
    gene.diff.expr.probe <- gene.diff.expr.probe[ 
      order( gene.diff.expr.fdr ) ]
    gene.diff.expr.fdr <- sort( gene.diff.expr.fdr )
    
    gene.diff.expr.fprtpr <- sapply( 1 : length( gene.diff.expr.probe ), 
                                     function( i ) {
                                       gene.diff.expr.pos <- gene.diff.expr.probe[ 1 : i ]
                                       c( mean( ma.probe.known.neg %in% gene.diff.expr.pos ), 
                                          mean( ma.probe.known.pos %in% gene.diff.expr.pos ) )
                                     } )
    
    lines( gene.diff.expr.fprtpr[ 1, ], gene.diff.expr.fprtpr[ 2, ], lwd=.5, 
           lty = norm.method.lty[ gdemeth ], col = norm.method.col[ gdemeth ] )
    
    fdr.point <- sapply( plot.fdr.point, function( pfp ) 
      sum( gene.diff.expr.fdr < pfp ) )
    
    points( gene.diff.expr.fprtpr[ 1, fdr.point ], 
            gene.diff.expr.fprtpr[ 2, fdr.point ], 
            pch = norm.method.pch[ gdemeth ], col = norm.method.col[ gdemeth ], 
            bg = norm.method.bg.col[ gdemeth ],cex=.8 )
  }
}
abline(v=c(0.01))

# print probe numbers

cat( sprintf( "number of platform probes: %d\n", 
              sum( ! ma.probe.control.bol ) ) )
cat( sprintf( "number of spike-ins: %d\n", 
              sum( exp.design$fc != -1 & ! ma.probe.control.bol ) ) )
cat( sprintf( "number of positive spike-ins: %d\n", 
              sum( exp.design$fc > 1 & ! ma.probe.control.bol ) ) )
cat( sprintf( "number of negative spike-ins: %d\n", 
              sum( exp.design$fc == 1 & ! ma.probe.control.bol ) ) )

cat( sprintf( "number of probes: %d\n", ma.probe.num ) )
cat( sprintf( "number of unknown probes: %d\n", length( ma.probe.unknown ) ) )
cat( sprintf( "number of known probes: %d\n", length( ma.probe.known ) ) )
cat( sprintf( "number of known positives: %d\n", 
              length( ma.probe.known.pos ) ) )
cat( sprintf( "number of known negatives: %d\n", 
              length( ma.probe.known.neg ) ) )

h0.probe.norm.method <- c( "stdvec.condec","median.condec")

for ( hpnmeth in h0.probe.norm.method )
{
  norm.result <- get( paste0( hpnmeth, ".norm.result" ) )
  norm.h0.probe <- norm.result$h0.probe
  
  cat( sprintf( "%s - h0 probes - number: %d\n", hpnmeth, 
                length( norm.h0.probe ) ) )
  cat( sprintf( "%s - h0 probes - fraction of known probes: %g\n", hpnmeth, 
                mean( norm.h0.probe %in% ma.probe.known ) ) )
  cat( sprintf( 
    "%s - h0 probes - fraction of known negatives among known probes: %g\n", 
    hpnmeth, sum( norm.h0.probe %in% ma.probe.known.neg ) / 
      sum( norm.h0.probe %in% ma.probe.known ) ) )
}



qqplot(log10(stdvec.condec.norm.MCT.gene.diff.expr$spike.vs.constant.fdr),
       log10(stdvec.condec.norm.ttest.gene.diff.expr$spike.vs.constant.fdr), xlab="MCT log10(p-value)", ylab="Welch's t log10(p-value)")
abline(0,1)
############################################################

plot.fdr.point <- c( 0.01)
plot.norm.method <- c( "stdvec.condec","median.condec")
gene.diff.expr.method <- c( "MCT","ttest" )

plot( 1, type="n", cex.lab=0.67, cex.axis=0.67, 
      xlim = c( 0, .2 ), ylim = c( 0, 1 ), 
      xlab="False positive rate", 
      ylab="True positive rate" )


for ( gdemeth in gene.diff.expr.method )
{
  for ( nmeth in plot.norm.method )
  {
    gene.diff.expr <- get( paste0( nmeth, ".norm.", gdemeth, 
                                   ".gene.diff.expr" ) )
    
    gene.diff.expr.probe <- gene.diff.expr$probe
    gene.diff.expr.fdr <- gene.diff.expr$spike.vs.constant.fdr
    
    gene.diff.expr.probe <- gene.diff.expr.probe[ 
      order( gene.diff.expr.fdr ) ]
    gene.diff.expr.fdr <- sort( gene.diff.expr.fdr )
    
    gene.diff.expr.fprtpr <- sapply( 1 : length( gene.diff.expr.probe ), 
                                     function( i ) {
                                       gene.diff.expr.pos <- gene.diff.expr.probe[ 1 : i ]
                                       c( mean( ma.probe.known.neg %in% gene.diff.expr.pos ), 
                                          mean( ma.probe.known.pos %in% gene.diff.expr.pos ) )
                                     } )
    
    lines( gene.diff.expr.fprtpr[ 1, ], gene.diff.expr.fprtpr[ 2, ], lwd=1, 
           lty = norm.method.lty[ gdemeth ], col = norm.method.col[ gdemeth ] )
    
    fdr.point <- sapply( plot.fdr.point, function( pfp ) 
      sum( gene.diff.expr.fdr < pfp ) )
    
    points( gene.diff.expr.fprtpr[ 1, fdr.point ], 
            gene.diff.expr.fprtpr[ 2, fdr.point ], 
            pch = norm.method.pch[ gdemeth ], col = norm.method.col[ gdemeth ], 
            bg = norm.method.bg.col[ gdemeth ],cex=1 )
    
    cat( sprintf( "number of probes: %d\n", fdr.point ) )
    cat( sprintf( "False positive rate: %f\n", gene.diff.expr.fprtpr[ 1, fdr.point ] ) )
    cat( sprintf( "True positive rate: %f\n", gene.diff.expr.fprtpr[ 2, fdr.point ] ) )
  }
}
abline(v=c(0.01))


qqplot(log10(stdvec.condec.norm.MCT.gene.diff.expr$spike.vs.constant.fdr),
       log10(stdvec.condec.norm.ttest.gene.diff.expr$spike.vs.constant.fdr), xlab="MCT log10(p-value)", ylab="Welch's t log10(p-value)")
abline(0,1)
############################################################

#adjusted nominal for Welch's test
plot.fdr.point <- c( 0.018)
plot.norm.method <- c( "stdvec.condec","median.condec")
gene.diff.expr.method <- c( "MCT","ttest" )

plot( 1, type="n", cex.lab=0.67, cex.axis=0.67, 
      xlim = c( 0, .2 ), ylim = c( 0, 1 ), 
      xlab="False positive rate", 
      ylab="True positive rate" )


for ( gdemeth in gene.diff.expr.method )
{
  for ( nmeth in plot.norm.method )
  {
    gene.diff.expr <- get( paste0( nmeth, ".norm.", gdemeth, 
                                   ".gene.diff.expr" ) )
    
    gene.diff.expr.probe <- gene.diff.expr$probe
    gene.diff.expr.fdr <- gene.diff.expr$spike.vs.constant.fdr
    
    gene.diff.expr.probe <- gene.diff.expr.probe[ 
      order( gene.diff.expr.fdr ) ]
    gene.diff.expr.fdr <- sort( gene.diff.expr.fdr )
    
    gene.diff.expr.fprtpr <- sapply( 1 : length( gene.diff.expr.probe ), 
                                     function( i ) {
                                       gene.diff.expr.pos <- gene.diff.expr.probe[ 1 : i ]
                                       c( mean( ma.probe.known.neg %in% gene.diff.expr.pos ), 
                                          mean( ma.probe.known.pos %in% gene.diff.expr.pos ) )
                                     } )
    
    lines( gene.diff.expr.fprtpr[ 1, ], gene.diff.expr.fprtpr[ 2, ], lwd=1, 
           lty = norm.method.lty[ gdemeth ], col = norm.method.col[ gdemeth ] )
    
    fdr.point <- sapply( plot.fdr.point, function( pfp ) 
      sum( gene.diff.expr.fdr < pfp ) )
    
    points( gene.diff.expr.fprtpr[ 1, fdr.point ], 
            gene.diff.expr.fprtpr[ 2, fdr.point ], 
            pch = norm.method.pch[ gdemeth ], col = norm.method.col[ gdemeth ], 
            bg = norm.method.bg.col[ gdemeth ],cex=1 )
    
    cat( sprintf( "number of probes: %d\n", fdr.point ) )
    cat( sprintf( "False positive rate: %f\n", gene.diff.expr.fprtpr[ 1, fdr.point ] ) )
    cat( sprintf( "True positive rate: %f\n", gene.diff.expr.fprtpr[ 2, fdr.point ] ) )
  }
}
abline(v=c(0.01))


########################################################
#adjusted nominal for MCT
plot.fdr.point <- c( 0.01, 0.037)
plot.norm.method <- c( "stdvec.condec","median.condec")
gene.diff.expr.method <- c( "MCT","ttest" )

plot( 1, type="n", cex.lab=0.67, cex.axis=0.67, 
      xlim = c( 0, .2 ), ylim = c( 0, 1 ), 
      xlab="False positive rate", 
      ylab="True positive rate" )


for ( gdemeth in gene.diff.expr.method )
{
  for ( nmeth in plot.norm.method )
  {
    gene.diff.expr <- get( paste0( nmeth, ".norm.", gdemeth, 
                                   ".gene.diff.expr" ) )
    
    gene.diff.expr.probe <- gene.diff.expr$probe
    gene.diff.expr.fdr <- gene.diff.expr$spike.vs.constant.fdr
    
    gene.diff.expr.probe <- gene.diff.expr.probe[ 
      order( gene.diff.expr.fdr ) ]
    gene.diff.expr.fdr <- sort( gene.diff.expr.fdr )
    
    gene.diff.expr.fprtpr <- sapply( 1 : length( gene.diff.expr.probe ), 
                                     function( i ) {
                                       gene.diff.expr.pos <- gene.diff.expr.probe[ 1 : i ]
                                       c( mean( ma.probe.known.neg %in% gene.diff.expr.pos ), 
                                          mean( ma.probe.known.pos %in% gene.diff.expr.pos ) )
                                     } )
    
    lines( gene.diff.expr.fprtpr[ 1, ], gene.diff.expr.fprtpr[ 2, ], lwd=1, 
           lty = norm.method.lty[ gdemeth ], col = norm.method.col[ gdemeth ] )
    
    fdr.point <- sapply( plot.fdr.point, function( pfp ) 
      sum( gene.diff.expr.fdr < pfp ) )
    
    points( gene.diff.expr.fprtpr[ 1, fdr.point ], 
            gene.diff.expr.fprtpr[ 2, fdr.point ], 
            pch = norm.method.pch[ gdemeth ], col = norm.method.col[ gdemeth ], 
            bg = norm.method.bg.col[ gdemeth ],cex=1 )
    
    cat( sprintf( "number of probes: %d\n", fdr.point ) )
    cat( sprintf( "False positive rate: %f\n", gene.diff.expr.fprtpr[ 1, fdr.point ] ) )
    cat( sprintf( "True positive rate: %f\n", gene.diff.expr.fprtpr[ 2, fdr.point ] ) )
  }
}
abline(v=c(0.01,.05))


########################################################
plot.fdr.point <- c( 0.01, 0.05)
gdemeth="ttest"
nmeth="stdvec.condec"

gene.diff.expr <- get( paste0( nmeth, ".norm.", gdemeth, 
                               ".gene.diff.expr" ) )

gene.diff.expr.probe <- gene.diff.expr$probe
gene.diff.expr.fdr <- gene.diff.expr$spike.vs.constant.fdr

gene.diff.expr.probe <- gene.diff.expr.probe[ 
  order( gene.diff.expr.fdr ) ]
gene.diff.expr.fdr <- sort( gene.diff.expr.fdr )

i=555

gene.diff.expr.pos <- gene.diff.expr.probe[ 1 : i ]
c( mean( ma.probe.known.neg %in% gene.diff.expr.pos ), 
   mean( ma.probe.known.pos %in% gene.diff.expr.pos ) )

gene.diff.expr.pos <- gene.diff.expr.probe[ 1 : i ]
c( sum( ma.probe.known.neg %in% gene.diff.expr.pos ), 
   sum( ma.probe.known.pos %in% gene.diff.expr.pos ) )

i=711

gene.diff.expr.pos <- gene.diff.expr.probe[ 1 : i ]
c( mean( ma.probe.known.neg %in% gene.diff.expr.pos ), 
   mean( ma.probe.known.pos %in% gene.diff.expr.pos ) )

gene.diff.expr.pos <- gene.diff.expr.probe[ 1 : i ]
c( sum( ma.probe.known.neg %in% gene.diff.expr.pos ), 
   sum( ma.probe.known.pos %in% gene.diff.expr.pos ) )
##########################################

plot.fdr.point <- c( 0.01, 0.05)
gdemeth="MCT"
nmeth="stdvec.condec"

gene.diff.expr <- get( paste0( nmeth, ".norm.", gdemeth, 
                               ".gene.diff.expr" ) )

gene.diff.expr.probe <- gene.diff.expr$probe
gene.diff.expr.fdr <- gene.diff.expr$spike.vs.constant.fdr

gene.diff.expr.probe <- gene.diff.expr.probe[ 
  order( gene.diff.expr.fdr ) ]
gene.diff.expr.fdr <- sort( gene.diff.expr.fdr )


i=744

gene.diff.expr.pos <- gene.diff.expr.probe[ 1 : i ]
c( mean( ma.probe.known.neg %in% gene.diff.expr.pos ), 
   mean( ma.probe.known.pos %in% gene.diff.expr.pos ) )

c( sum( ma.probe.known.neg %in% gene.diff.expr.pos ), 
   sum( ma.probe.known.pos %in% gene.diff.expr.pos ) )


############################################################

plot.norm.method <- c("stdvec.condec")
gene.diff.expr.method <- c( "MCT")

for ( gdemeth in gene.diff.expr.method )
{
  for ( nmeth in plot.norm.method )
  {
    gene.diff.expr <- get( paste0( nmeth, ".norm.", gdemeth, 
                                   ".gene.diff.expr" ) )
    
    gene.diff.expr.probe <- gene.diff.expr$probe
    gene.diff.expr.fdr <- gene.diff.expr$spike.vs.constant.fdr
    
    gene.diff.expr.probe <- gene.diff.expr.probe[ 
      order( gene.diff.expr.fdr ) ]
    gene.diff.expr.fdr <- sort( gene.diff.expr.fdr )
    
    gene.diff.expr.fprtpr <- sapply( 1 : length( gene.diff.expr.probe ), 
                                     function( i ) {
                                       gene.diff.expr.pos <- gene.diff.expr.probe[ 1 : i ]
                                       c( mean( ma.probe.known.neg %in% gene.diff.expr.pos ), 
                                          mean( ma.probe.known.pos %in% gene.diff.expr.pos ) )
                                     } )
    
    #hist(gene.diff.expr$spike.vs.constant.fdr[gene.diff.expr$probe %in% ma.probe.known.neg])
    aa <- log10(runif(100000))
    mctkn <- log10(gene.diff.expr$spike.vs.constant.fdr[gene.diff.expr$probe %in% ma.probe.known.neg])
    qqplot(aa,mctkn,xlab="log10(unif)",ylab="log10(MCT)")
    abline(0,1)

  }
}

############################################
plot.norm.method <- c("stdvec.condec")
gene.diff.expr.method <- c( "ttest" )

for ( gdemeth in gene.diff.expr.method )
{
  for ( nmeth in plot.norm.method )
  {
    gene.diff.expr <- get( paste0( nmeth, ".norm.", gdemeth, 
                                   ".gene.diff.expr" ) )
    
    gene.diff.expr.probe <- gene.diff.expr$probe
    gene.diff.expr.fdr <- gene.diff.expr$spike.vs.constant.fdr
    
    gene.diff.expr.probe <- gene.diff.expr.probe[ 
      order( gene.diff.expr.fdr ) ]
    gene.diff.expr.fdr <- sort( gene.diff.expr.fdr )
    
    gene.diff.expr.fprtpr <- sapply( 1 : length( gene.diff.expr.probe ), 
                                     function( i ) {
                                       gene.diff.expr.pos <- gene.diff.expr.probe[ 1 : i ]
                                       c( mean( ma.probe.known.neg %in% gene.diff.expr.pos ), 
                                          mean( ma.probe.known.pos %in% gene.diff.expr.pos ) )
                                     } )
    
    #hist(gene.diff.expr$spike.vs.constant.fdr[gene.diff.expr$probe %in% ma.probe.known.neg])
    aa <- log10(runif(100000))
    wkn <- log10(gene.diff.expr$spike.vs.constant.fdr[gene.diff.expr$probe %in% ma.probe.known.neg])
    qqplot(aa,wkn,xlab="log10(unif)",ylab="log10(Welch)")
    abline(0,1)
    
  }
}

############################################################

plot.norm.method <- c("stdvec.condec")
gene.diff.expr.method <- c( "MCT", "ttest")

for ( nmeth in plot.norm.method )
{
  for ( gdemeth in gene.diff.expr.method )
  {
    gene.diff.expr <- get( paste0( nmeth, ".norm.", gdemeth, 
                                   ".gene.diff.expr" ) )
    
    gene.diff.expr.probe <- gene.diff.expr$probe
    gene.diff.expr.fdr <- gene.diff.expr$spike.vs.constant.fdr
    
    gene.diff.expr.probe <- gene.diff.expr.probe[ 
      order( gene.diff.expr.fdr ) ]
    gene.diff.expr.fdr <- sort( gene.diff.expr.fdr )
    
    gene.diff.expr.fprtpr <- sapply( 1 : length( gene.diff.expr.probe ), 
                                     function( i ) {
                                       gene.diff.expr.pos <- gene.diff.expr.probe[ 1 : i ]
                                       c( mean( ma.probe.known.neg %in% gene.diff.expr.pos ), 
                                          mean( ma.probe.known.pos %in% gene.diff.expr.pos ) )
                                     } )
    
    #hist(gene.diff.expr$spike.vs.constant.fdr[gene.diff.expr$probe %in% ma.probe.known.neg])
    assign( paste0( gdemeth, ".gene.diff.expr" ), gene.diff.expr )
  }
  
  mctkp <- log10(MCT.gene.diff.expr$spike.vs.constant.fdr[MCT.gene.diff.expr$probe %in% ma.probe.known.pos])
  wkp <- log10(ttest.gene.diff.expr$spike.vs.constant.fdr[ttest.gene.diff.expr$probe %in% ma.probe.known.pos])
  qqplot(wkp,mctkp,ylab="log10(MCT)",xlab="log10(Welch)")
  abline(0,1)
}


#################################################################

pdf(file="QQGolds.pdf",width=10, height=14)
par(mfrow=c(3,1), mar=c(6, 6, 4, 2))

set.seed(10)
tu <- log10(runif(100000))
qqplot(-tu,-wkn,col=1,lwd=2,xlab="", ylab="", cex.axis=2, las=1)
abline(0,1)

title(main="(a)",cex.main=2)
title(xlab="-log10(uniform)", cex.lab=2, line = 3.5)
title(ylab="-log10(q-value) W", cex.lab=2, line=4)


qqplot(-tu,-mctkn,col=1,lwd=2,xlab="", ylab="", cex.axis=2, las=1)
abline(0,1)

title(main="(b)",cex.main=2)
title(xlab="-log10(uniform)", cex.lab=2, line = 3.5)
title(ylab="-log10(q-value) MCT", cex.lab=2, line=4)


qqplot(-wkp,-mctkp,col=1,lwd=2,xlab="", ylab="", cex.axis=2, las=1)
abline(0,1)

title(main="(c)",cex.main=2)
title(xlab="-log10(q-value) W", cex.lab=2, line = 3.5)
title(ylab="-log10(q-value) MCT", cex.lab=2, line=4)

dev.off()

#################################################################

postscript(file="QQGolda.eps")
par(mar=c(6, 6, 4, 2))
set.seed(10)
tu <- log10(runif(100000))
qqplot(-tu,-wkn,col=1,lwd=2,xlab="", ylab="", cex.axis=2, las=1)
abline(0,1)

title(main="(a)",cex.main=2)
title(xlab="-log10(uniform)", cex.lab=2, line = 3.5)
title(ylab="-log10(q-value) W", cex.lab=2, line=4)
dev.off()

postscript(file="QQGoldb.eps")
par(mar=c(6, 6, 4, 2))
qqplot(-tu,-mctkn,col=1,lwd=2,xlab="", ylab="", cex.axis=2, las=1)
abline(0,1)

title(main="(b)",cex.main=2)
title(xlab="-log10(uniform)", cex.lab=2, line = 3.5)
title(ylab="-log10(q-value) MCT", cex.lab=2, line=4)
dev.off()

postscript(file="QQGoldc.eps")
par(mar=c(6, 6, 4, 2))
qqplot(-wkp,-mctkp,col=1,lwd=2,xlab="", ylab="", cex.axis=2, las=1)
abline(0,1)

title(main="(c)",cex.main=2)
title(xlab="-log10(q-value) W", cex.lab=2, line = 3.5)
title(ylab="-log10(q-value) MCT", cex.lab=2, line=4)
dev.off()