# Copyright (c) 2016, Universitat Rovira i Virgili (Spain), Aarhus University 
# (Denmark) and University of Aveiro (Portugal)
# 
# Written by Carlos P. Roca
# as Research funded by the European Union
# for the research paper by Roca, Gomes, Amorim & Scott-Fordsmand: "Variation-
# preserving normalization unveils blind spots in gene expression profiling".
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


# Implements median condition-decomposition normalization


normalize.median.condec <- function( expression.data, expression.condition, 
    normalize.probe=NULL, convergence.threshold = c( 0.001, 0.1 ), 
    search.h0.probe=TRUE, norm.probability=0.5, p.value.graph=NULL, 
    verbose=FALSE )
{
    # identify samples and conditions to normalize
    normalize.sample <- 
        colnames( expression.data )[ ! is.na( expression.condition ) ]
    normalize.sample.condition <- 
        as.character( na.omit( expression.condition ) )
    normalize.condition <- unique( normalize.sample.condition )
    
    if ( length( normalize.condition ) == 0 )
        stop( "No condition to normalize in normalize.median.condec" )
    else if ( min( table( normalize.sample.condition ) ) < 2 )
        stop( "There must be 2 or more samples for each condition in normalize.median.condec" )
    
    # select probes for normalization
    expression.probe <- rownames( expression.data )
    
    if ( is.null( normalize.probe ) )
        normalize.probe.idx <- 1 : length( expression.probe )
    else
    {
        normalize.probe.idx <- match( normalize.probe, expression.probe )
        if ( any( is.na( normalize.probe.idx ) ) )
            stop( "Bad normalize.probe argument in normalize.median.condec" )
    }
    
    # enforce no missing values in any normalization sample
    all.sample.probe.idx <- which( rowSums( 
        is.na( expression.data[ , normalize.sample ] ) ) == 0 )
    
    normalize.probe.idx <- 
        intersect( normalize.probe.idx, all.sample.probe.idx )
    
    # normalize within conditions

    normalize.expr.data <- matrix( nrow = length( expression.probe ), 
        ncol = length( normalize.sample ) )
    rownames( normalize.expr.data ) <- expression.probe
    colnames( normalize.expr.data ) <- normalize.sample
    
    normalize.within.cond.offset <- rep( 0, length( normalize.sample ) )
    names( normalize.within.cond.offset ) <- normalize.sample
    
    condition.sample.idx <- split( 1 : length( expression.condition ), 
        expression.condition )
    
    for ( sample.idx in condition.sample.idx )
    {
        condition <- as.character( expression.condition[ sample.idx[ 1 ] ] )
        norm.sample.idx <- which( condition == normalize.sample.condition )
        
        within.cond.norm.result <- normalize.median.within.condition( 
            expression.data[ , sample.idx ], condition, normalize.probe.idx, 
            norm.probability, verbose )
        
        normalize.expr.data[ , norm.sample.idx ] <- within.cond.norm.result$data
        
        normalize.within.cond.offset[ norm.sample.idx ] <- 
            within.cond.norm.result$offset
    }
    
    if ( length( normalize.condition ) == 1 )
    {
        # no additional normalization needed
        normalize.offset <- normalize.within.cond.offset[ normalize.sample ]
        
        normalize.result <- list( data = normalize.expr.data, 
            offset = normalize.offset )
        
        return( normalize.result )
    }
    
    # normalize between conditions
    
    between.cond.norm.result <- normalize.median.between.condition( 
        normalize.expr.data, normalize.condition, normalize.sample.condition, 
        normalize.probe.idx, convergence.threshold, search.h0.probe, 
        norm.probability, p.value.graph, verbose )
    
    normalize.expr.data <- between.cond.norm.result$data
    normalize.between.cond.offset <- between.cond.norm.result$offset
    
    if ( search.h0.probe )
    {
        normalize.between.cond.h0.probe <- between.cond.norm.result$h0.probe
        normalize.between.cond.h0.probe.convergence <- 
            between.cond.norm.result$h0.probe.convergence
    }
    
    normalize.offset <- normalize.within.cond.offset[ normalize.sample ] + 
        normalize.between.cond.offset[ normalize.sample.condition ]
    
    normalize.result <- list( data = normalize.expr.data, 
        offset = normalize.offset, 
        within.condition.offset = normalize.within.cond.offset, 
        between.condition.offset = normalize.between.cond.offset )
    
    if ( search.h0.probe )
    {
        normalize.result$h0.probe <- normalize.between.cond.h0.probe
        normalize.result$h0.probe.convergence <- 
            normalize.between.cond.h0.probe.convergence
    }
    
    normalize.result
}


normalize.median.within.condition <- function( edata, condition, 
    norm.probe.idx, norm.prob, verbose )
{
    if ( verbose )
        cat( paste0( condition, "\n" ) )
    
    # calculate normalization
    norm.median.offset <- apply( edata[ norm.probe.idx, ], 2, quantile, 
        probs=norm.prob )
    norm.median.offset <- norm.median.offset - mean( norm.median.offset )
    
    # normalize with obtained offset
    edata <- sweep( edata, 2, norm.median.offset )
    
    norm.within.cond.result <- list( data = edata, offset = norm.median.offset )
    
    norm.within.cond.result
}


normalize.median.between.condition <- function( edata, norm.cond, 
    norm.sample.cond, norm.probe.idx, convergence.threshold, search.h0.probe, 
    norm.prob, p.value.graph, verbose )
{
    # identify samples available per condition
    within.cond.n <- as.vector( table( norm.sample.cond )[ norm.cond ] )
    names( within.cond.n ) <- norm.cond
    
    # calculate balanced within-condition means
    bal.mean.n <- min( within.cond.n )
    
    expr.bal.mean.data <- sapply( norm.cond, function( cond ) 
        if ( within.cond.n[ cond ] == bal.mean.n )
            rowMeans( edata[ norm.probe.idx, cond == norm.sample.cond ] )
        else  # within.cond.n[ cond ] > bal.mean.n
            rowMeans( t( apply( 
                t( edata[ norm.probe.idx, cond == norm.sample.cond ] ), 
                2, function( ed ) sample( ed, bal.mean.n ) ) ) )
    )
    
    if ( search.h0.probe )
    {
        # calculate normalization while looking for h0 probes
        
        # calculate within-condition means for F-statistic
        expr.mean.data <- sapply( norm.cond, function( cond ) 
            rowMeans( edata[ norm.probe.idx, cond == norm.sample.cond ] ) )
        
        # calculate within-condition variances for F-statistic
        within.cond.var <- sweep( sapply( norm.cond, function( cond ) 
            apply( edata[ norm.probe.idx, cond == norm.sample.cond ], 1, var ) ), 
            2, within.cond.n - 1, "*" )
        within.cond.var <- rowSums( within.cond.var ) / 
            ( sum( within.cond.n ) - length( norm.cond ) )
        
        if ( verbose )
            cat( "between.condition.search.h0.probe\n" )
        
        norm.median.search.h0.probe.result <- normalize.median.selection( 
            expr.bal.mean.data, convergence.threshold, expr.mean.data, 
            within.cond.var, within.cond.n, norm.prob, p.value.graph, verbose )
        
        h0.probe <- norm.median.search.h0.probe.result$h0.probe
        h0.probe.convergence <- norm.median.search.h0.probe.result$convergence
        
        norm.probe <- h0.probe
    }
    else
        # uses all probes given by norm.probe.idx for normalization
        norm.probe <- rownames( expr.bal.mean.data )
    
    if ( verbose )
        cat( "between.condition\n" )

    # calculate normalization
    norm.median.offset <- apply( expr.bal.mean.data[ norm.probe, ], 2, 
        quantile, probs=norm.prob )
    norm.median.offset <- norm.median.offset - mean( norm.median.offset )
    
    # normalize each condition with obtained offset
    for( cond in norm.cond )
        edata[ , cond == norm.sample.cond ] <- 
            edata[ , cond == norm.sample.cond ] - norm.median.offset[ cond ]
    
    norm.between.cond.result <- list( data = edata, 
        offset = norm.median.offset )
    
    if( search.h0.probe )
    {
        norm.between.cond.result$h0.probe <- h0.probe
        norm.between.cond.result$h0.probe.convergence <- h0.probe.convergence
    }
    
    norm.between.cond.result
}


normalize.median.selection <- function( edata, convergence.threshold, 
    edata.fstat, within.cond.var, within.cond.n, norm.prob, p.value.graph, 
    verbose )
{
    iter.max <- 200
    single.threshold <- convergence.threshold[ 1 ]
    accum.threshold <- convergence.threshold[ 2 ]
    offset.accum.step.threshold <- 10
    common.h0.probe.accum.step.max <- 10
    
    median.offset <- rep( 0, ncol( edata ) )
    
    norm.median.offset <- matrix( nrow=0, ncol = ncol( edata ) )
    norm.median.offset.sd <- vector( "numeric" )
    norm.median.offset.delta.sd <- vector( "numeric" )
    norm.median.offset.accum.step <- vector( "numeric" )
    norm.median.common.h0.probe.accum.step <- vector( "numeric" )
    norm.median.h0.probe.num <- vector( "numeric" )
    norm.median.common.h0.probe.num <- vector( "numeric" )
    
    last.norm.data <- edata
    last.norm.data.fstat <- edata.fstat
    
    median.offset.accum.step <- 0
    median.common.h0.probe.accum.step <- 0
    
    iter <- 0
    overall.convergence <- FALSE
    
    while ( iter < iter.max && ! overall.convergence )
    {
        iter <- iter + 1
        
        # obtain next step of median offset     
        median.offset.step <- calculate.median.offset( last.norm.data, 
            last.norm.data.fstat, within.cond.var, within.cond.n, norm.prob, 
            p.value.graph, iter )
        
        median.offset.delta <- median.offset.step$value
        
        median.h0.probe <- median.offset.step$h0.probe
        median.h0.probe.num <- length( median.h0.probe )
        
        # check errors
        if ( any( is.nan( median.offset.delta ) ) )
            stop( "NaN error in normalize.median.condec" )
        
        # update total median offset
        median.offset <- median.offset + median.offset.delta
        median.offset <- median.offset - mean( median.offset )
        
        # update data at once with total offset
        last.norm.data <- sweep( edata, 2, median.offset )
        last.norm.data.fstat <- sweep( edata.fstat, 2, median.offset )
        
        # check convergence
        median.offset.sd <- sd( median.offset )
        median.offset.delta.sd <- sd( median.offset.delta )
        
        median.offset.delta.sd.ratio <- ifelse( median.offset.sd > 0, 
            median.offset.delta.sd / median.offset.sd, 1 )
        
        if ( median.offset.delta.sd == 0 )
            median.offset.ratio <- 0
        else if ( median.offset.sd == 0 )
            median.offset.ratio <- 1
        else
            median.offset.ratio <- median.offset.delta.sd / median.offset.sd
        
        median.offset.accum.step <- ifelse( 
            median.offset.ratio < accum.threshold, 
            median.offset.accum.step + 1, 0 )
        
        median.convergence <- median.offset.ratio < single.threshold || 
            median.offset.accum.step > offset.accum.step.threshold
        
        median.common.h0.probe.accum.step <- ifelse( median.convergence, 
            median.common.h0.probe.accum.step + 1, 0 )
        
        overall.convergence <- median.common.h0.probe.accum.step >= 
            common.h0.probe.accum.step.max
        
        # store last results
        last.median.offset <- median.offset
        
        if ( median.common.h0.probe.accum.step == 0 )
            median.common.h0.probe <- NULL
        else if ( median.common.h0.probe.accum.step == 1 )
            median.common.h0.probe <- median.h0.probe
        else
            median.common.h0.probe <- intersect( median.h0.probe, 
                median.common.h0.probe )
        
        median.common.h0.probe.num <- ifelse( 
            is.null( median.common.h0.probe ), NA, 
            length( median.common.h0.probe ) )
        
        # store step results
        norm.median.offset <- rbind( norm.median.offset, median.offset )
        norm.median.offset.sd <- c( norm.median.offset.sd, median.offset.sd )
        norm.median.offset.delta.sd <- c( norm.median.offset.delta.sd, 
            median.offset.delta.sd )
        norm.median.offset.accum.step <- c( norm.median.offset.accum.step, 
            median.offset.accum.step )
        norm.median.common.h0.probe.accum.step <- 
            c( norm.median.common.h0.probe.accum.step, 
                median.common.h0.probe.accum.step )
        norm.median.h0.probe.num <- c( norm.median.h0.probe.num, 
            median.h0.probe.num )
        norm.median.common.h0.probe.num <- c( norm.median.common.h0.probe.num, 
            median.common.h0.probe.num )
        
        if ( verbose )
        {
            cat( sprintf( "  %2d %g %g %02d %02d %d", iter, 
                median.offset.sd, median.offset.delta.sd.ratio, 
                median.offset.accum.step, median.common.h0.probe.accum.step, 
                median.h0.probe.num ) )
            
            if ( ! is.na( median.common.h0.probe.num ) )
                cat( paste0( " ", median.common.h0.probe.num ) )
            
            cat( "\n" )
        }
    }
    
    if ( ! overall.convergence )
        stop( "No convergence in normalize.median.condec" )
    
    # remove sample or condition names from step results
    dimnames( norm.median.offset ) <- NULL
    
    norm.median.convergence <- list( offset = norm.median.offset, 
        offset.sd = norm.median.offset.sd, 
        offset.delta.sd = norm.median.offset.delta.sd, 
        offset.accum.step = norm.median.offset.accum.step, 
        common.h0.probe.accum.step = norm.median.common.h0.probe.accum.step, 
        h0.probe.num = norm.median.h0.probe.num, 
        common.h0.probe.num = norm.median.common.h0.probe.num )
    
    norm.median.result <- list( offset = last.median.offset, 
        convergence = norm.median.convergence, 
        h0.probe = median.common.h0.probe )
    
    norm.median.result
}


calculate.median.offset <- function( edata, edata.fstat, within.cond.var, 
    within.cond.n, norm.prob, p.value.graph, iter )
{
    ks.test.alpha <- 1e-3
    
    # calculate f statistics for each probe
    expr.k <- length( within.cond.n )
    expr.n <- sum( within.cond.n )
    
    expr.grand.mean <- apply( edata.fstat, 1, function( ef ) 
        sum( ef * within.cond.n ) ) / expr.n
    
    between.cond.var <- apply( ( edata.fstat - expr.grand.mean )^2, 1, 
        function( ef2 ) sum( ef2 * within.cond.n ) ) / ( expr.k - 1 )
    
    expr.f <- between.cond.var / within.cond.var
    
    expr.f <- na.omit( expr.f )  # in case of 0/0
    attr( expr.f, "na.action" ) <- NULL
    
    expr.p.value <- pf( expr.f, df1 = expr.k-1, df2 = expr.n-expr.k, 
        lower.tail = FALSE )
    
    # identify h0 probes with one-sided up Kolmogorov-Smirnov test
    ks.test.d <- sqrt( - log( ks.test.alpha ) / 2 )
    
    epv <- sort( expr.p.value )
    epv.n <- length( expr.p.value )
    
    epv.i <- 1
    ks.test.D.up <- 1:epv.n / epv.n - epv
    ks.test.reject <- any( ks.test.D.up > ks.test.d / sqrt( epv.n ) )
    
    while ( ks.test.reject )
    {
        epv.i <- epv.i + 1
        
        ks.test.D.up <- 
            ( epv.i : epv.n - epv.i + 1 ) / ( epv.n - epv.i + 1 ) - 
            ( epv[ epv.i : epv.n ] - epv[ epv.i - 1 ] ) / 
            ( 1 - epv[ epv.i - 1 ] )
        
        ks.test.reject <- 
            any( ks.test.D.up > ks.test.d / sqrt( epv.n - epv.i + 1 ) )
    }
    
    epv.h0.i <- epv.i
    epv.h0.n <- epv.n - epv.i + 1
    epv.h0.p <-  ifelse( epv.i == 1, 0, epv[ epv.i - 1 ] )
    epv.h0.q <- ( epv.i - 1 ) / epv.n
    
    h0.probe.idx <- order( expr.p.value, decreasing=TRUE )[ 1 : epv.h0.n ]
    h0.probe <- names( expr.p.value )[ h0.probe.idx ]
    
    pi0.est <- ( 1 - epv.h0.q ) / ( 1 - epv.h0.p )
    
    # plot graph of p-values
    if ( ! is.null( p.value.graph ) )
    {
        if ( p.value.graph != "" )
        {
            if ( ! file.exists( p.value.graph ) )
                dir.create( p.value.graph, recursive=TRUE )
            
            png.filename <- sprintf( "%s/pvalue_iter%02d.png", 
                p.value.graph, iter )
            
            png( png.filename, width=1280, height=720 )
        }
        
        par.default <- par( no.readonly=TRUE )
        
        par( mfrow = c(1,2), pty="s", mar = c( 3.4, 4.0, 0, 1.6 ), 
            oma = c( 0, 4.2, 3.4, 0 ), mgp = c( 4.0, 1.2, 0 ) )
        
        epv.quant <- 1:epv.n / epv.n
        epv.h0.idx <- epv.h0.i : epv.n
        
        epv.x0 <- epv.h0.p
        epv.y0 <- epv.h0.q
        epv.yd <- ks.test.d * ( sqrt( epv.h0.n ) / epv.n )
        
        xylim <- list( c(0,0), c( epv.h0.p, epv.h0.q ) )
        
        for ( i in 1:2  )
        {
            plot( 0, type="n", 
                xlim = c( xylim[[i]][1], 1 ), ylim = c( xylim[[i]][2], 1 ), 
                xlab="p-value", ylab="", cex.axis=2.5, cex.lab=2.5 )
            
            segments( x0 = c( 0, epv[ - epv.h0.idx ] ), 
                y0 = c( 0, epv.quant[ - epv.h0.idx ] ), 
                x1 = c( epv[ - epv.h0.idx ], epv[ epv.h0.idx ][ 1 ] ), 
                lwd=2.5, col="black" )
            segments( x0 = epv[ epv.h0.idx ], y0 = epv.quant[ epv.h0.idx ], 
                x1 = c( epv[ epv.h0.idx ][ -1 ], 1 ), lwd=2.5, col="red" )
            points( epv.h0.p, epv.h0.q, pch=20, cex=3 )
            
            segments( 0, 1 - pi0.est, 1, 1, lwd=2, lty=2, col="blue" )
            segments( 0, 1 - pi0.est + epv.yd, 1, 1 + epv.yd, lwd=2, lty=3, 
                col="blue" )
        }
        
        mtext( "F( p-value )", side=2, line=1.0, outer=TRUE, cex=2.5 )
        
        graph.title <- substitute( 
            paste( "iter=", iter, "    -    #", H[0], "=", h0.n ), 
            list( iter = sprintf( "%02d", iter ), 
                h0.n = sprintf( "%5d", epv.h0.n ) ) )
        mtext( graph.title, 3, line=-1.2, outer=TRUE, cex=2.7 )
        
        if ( p.value.graph != "" )
            dev.off()
        else
            par( par.default )
    }
    
    # calculate offset
    median.offset <- apply( edata[ h0.probe, ], 2, quantile, probs=norm.prob )
    median.offset <- median.offset - mean( median.offset )
    
    median.result <- list( value = median.offset, h0.probe = h0.probe )
    
    median.result
}

