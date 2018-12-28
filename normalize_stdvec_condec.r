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


# Implements standard-vector condition-decompositon normalization


normalize.stdvec.condec <- function( expression.data, expression.condition, 
    normalize.probe=NULL, convergence.threshold = c( 0.01, 0.1, 0.01, 1 ), 
    search.h0.probe=TRUE, vector.graph=NULL, p.value.graph=NULL, verbose=FALSE )
{
    # identify samples and conditions to normalize
    normalize.sample <- 
        colnames( expression.data )[ ! is.na( expression.condition ) ]
    normalize.sample.condition <- 
        as.character( na.omit( expression.condition ) )
    normalize.condition <- unique( normalize.sample.condition )
    
    if ( length( normalize.condition ) == 0 )
        stop( "No condition to normalize in normalize.stdvec.condec" )
    else if ( min( table( normalize.sample.condition ) ) < 2 )
        stop( "There must be 2 or more samples for each condition in normalize.stdvec.condec" )
    
    # select probes for normalization
    expression.probe <- rownames( expression.data )
    
    if ( is.null( normalize.probe ) )
        normalize.probe.idx <- 1 : length( expression.probe )
    else
    {
        normalize.probe.idx <- match( normalize.probe, expression.probe )
        if ( any( is.na( normalize.probe.idx ) ) )
            stop( "Bad normalize.probe argument in normalize.stdvec.condec" )
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
    
    normalize.within.cond.convergence <- vector( "list", 
        length( normalize.condition ) )
    names( normalize.within.cond.convergence ) <- normalize.condition
    
    condition.sample.idx <- split( 1 : length( expression.condition ), 
        expression.condition )
    
    for ( sample.idx in condition.sample.idx )
    {
        condition <- as.character( expression.condition[ sample.idx[ 1 ] ] )
        norm.sample.idx <- which( condition == normalize.sample.condition )
        
        within.cond.norm.result <- normalize.stdvec.within.condition(
            expression.data[ , sample.idx ], condition, normalize.probe.idx, 
            convergence.threshold, vector.graph, verbose )
        
        normalize.expr.data[ , norm.sample.idx ] <- within.cond.norm.result$data
        
        normalize.within.cond.offset[ norm.sample.idx ] <-
            within.cond.norm.result$offset
        
        normalize.within.cond.convergence[[ condition ]] <-
            within.cond.norm.result$convergence
    }
    
    # remove condition names from convergence data
    names( normalize.within.cond.convergence ) <- NULL
    
    if ( length( normalize.condition ) == 1 )
    {
        # no additional normalization needed
        normalize.offset <- normalize.within.cond.offset[ normalize.sample ]
        
        normalize.result <- list( data = normalize.expr.data, 
            offset = normalize.offset, 
            convergence = normalize.within.cond.convergence )
        
        return( normalize.result )
    }
    
    # normalize between conditions
    
    between.cond.norm.result <- normalize.stdvec.between.condition( 
        normalize.expr.data, normalize.condition, normalize.sample.condition, 
        normalize.probe.idx, convergence.threshold, search.h0.probe, 
        vector.graph, p.value.graph, verbose )
    
    normalize.expr.data <- between.cond.norm.result$data
    normalize.between.cond.offset <- between.cond.norm.result$offset
    normalize.between.cond.convergence <- between.cond.norm.result$convergence
    
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
        between.condition.offset = normalize.between.cond.offset, 
        within.condition.convergence = normalize.within.cond.convergence, 
        between.condition.convergence = normalize.between.cond.convergence )
    
    if ( search.h0.probe )
    {
        normalize.result$h0.probe <- normalize.between.cond.h0.probe
        normalize.result$h0.probe.convergence <- 
            normalize.between.cond.h0.probe.convergence
    }
    
    normalize.result
}


normalize.stdvec.within.condition <- function( edata, condition, 
    norm.probe.idx, convergence.threshold, vector.graph, verbose )
{
    # calculate normalization
    norm.stdvec.result <- normalize.standard.vector( edata[ norm.probe.idx, ], 
        condition, convergence.threshold, FALSE, NULL, NULL, NULL, 
        vector.graph, NULL, verbose )
    
    # normalize with obtained offset
    edata <- sweep( edata, 2, norm.stdvec.result$offset )
    
    norm.within.cond.result <- list( data = edata, 
        offset = norm.stdvec.result$offset, 
        convergence = norm.stdvec.result$convergence )
    
    norm.within.cond.result
}


normalize.stdvec.between.condition <- function( edata, norm.cond, 
    norm.sample.cond, norm.probe.idx, convergence.threshold, search.h0.probe, 
    vector.graph, p.value.graph, verbose )
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
        
        if ( ! is.null( vector.graph ) && vector.graph != "" )
            vector.graph.h0.probe <- paste0( vector.graph, ".h0.probe" )
        else
            vector.graph.h0.probe <- vector.graph
        
        norm.stdvec.search.h0.probe.result <- normalize.standard.vector( 
            expr.bal.mean.data, "between.condition.search.h0.probe", 
            convergence.threshold, TRUE, expr.mean.data, within.cond.var, 
            within.cond.n, vector.graph.h0.probe, p.value.graph, verbose )
        
        h0.probe <- norm.stdvec.search.h0.probe.result$h0.probe
        h0.probe.convergence <- norm.stdvec.search.h0.probe.result$convergence
        
        norm.probe <- h0.probe
    }
    else
        # uses all probes given by norm.probe.idx for normalization
        norm.probe <- rownames( expr.bal.mean.data )

    # calculate normalization
    norm.stdvec.result <- normalize.standard.vector( 
        expr.bal.mean.data[ norm.probe, ], "between.condition", 
        convergence.threshold, FALSE, NULL, NULL, NULL, vector.graph, NULL, 
        verbose )
    
    # normalize each condition with obtained offset
    for( cond in norm.cond )
        edata[ , cond == norm.sample.cond ] <- 
            edata[ , cond == norm.sample.cond ] - 
            norm.stdvec.result$offset[ cond ]
    
    norm.between.cond.result <- list( data = edata, 
        offset = norm.stdvec.result$offset, 
        convergence = norm.stdvec.result$convergence )
    
    if( search.h0.probe )
    {
        norm.between.cond.result$h0.probe <- h0.probe
        norm.between.cond.result$h0.probe.convergence <- h0.probe.convergence
    }
    
    norm.between.cond.result
}


normalize.standard.vector <- function( edata, condition, convergence.threshold, 
    search.h0.probe, edata.fstat, within.cond.var, within.cond.n, vector.graph, 
    p.value.graph, verbose )
{
    if ( ! search.h0.probe )
    {
        iter.max <- 100
        single.threshold <- convergence.threshold[ 1 ]
        accum.threshold <- convergence.threshold[ 2 ]
        offset.accum.step.threshold <- 10
    }
    else
    {
        iter.max <- 200
        single.threshold <- convergence.threshold[ 3 ]
        accum.threshold <- convergence.threshold[ 4 ]
        offset.accum.step.threshold <- 10
        common.h0.probe.accum.step.max <- 10
    }
    
    vector.graph.probe.num <- 10000
    
    if ( verbose )
        cat( paste0( condition, "\n" ) )
    
    stdvec.offset <- rep( 0, ncol( edata ) )
    
    norm.stdvec.offset <- matrix( nrow=0, ncol = ncol( edata ) )
    norm.stdvec.offset.sd <- vector( "numeric" )
    norm.stdvec.offset.stderr <- vector( "numeric" )
    norm.stdvec.offset.delta.sd <- vector( "numeric" )
    norm.stdvec.offset.accum.step <- vector( "numeric" )
    norm.stdvec.numerical.demand <- vector( "numeric" )
    
    if ( ncol( edata ) >= 3 )
        norm.stdvec.watson.u2 <- matrix( nrow=0, 
            ncol = ( ncol( edata ) - 1 ) %/% 3 + 1 )
    else
        norm.stdvec.watson.u2 <- NULL
    
    if ( search.h0.probe )
    {
        norm.stdvec.common.h0.probe.accum.step <- vector( "numeric" )
        norm.stdvec.h0.probe.num <- vector( "numeric" )
        norm.stdvec.common.h0.probe.num <- vector( "numeric" )
    }
    
    if ( ! is.null( vector.graph ) )
    {
        # select a sample of probes for plotting standardized sample vectors
        edata.probe <- rownames( edata )
        edata.probe.num <- length( edata.probe )
        
        if ( edata.probe.num <= vector.graph.probe.num )
            vector.graph.probe <- edata.probe
        else
            vector.graph.probe <- sample( edata.probe, vector.graph.probe.num )
    }
    else
        vector.graph.probe <- NULL
    
    last.norm.data <- edata
    stdvec.offset.accum.step <- 0
    
    if ( search.h0.probe )
    {
        last.norm.data.fstat <- edata.fstat
        stdvec.common.h0.probe.accum.step <- 0
    }
    
    iter <- 0
    overall.convergence <- FALSE
    
    while ( iter < iter.max && ! overall.convergence )
    {
        iter <- iter + 1
        
        # obtain next step of standard vector offset
        if ( ! search.h0.probe )
            stdvec.offset.step <- calculate.stdvec.offset( last.norm.data, 
                condition, FALSE, NULL, NULL, NULL, vector.graph, 
                vector.graph.probe, NULL, iter )
        else
            stdvec.offset.step <- calculate.stdvec.offset( last.norm.data, 
                condition, TRUE, last.norm.data.fstat, within.cond.var, 
                within.cond.n, vector.graph, vector.graph.probe, p.value.graph, 
                iter )
        
        stdvec.offset.delta <- stdvec.offset.step$value
        stdvec.offset.stderr <- stdvec.offset.step$stderr
        stdvec.numerical.demand <- stdvec.offset.step$numerical.demand
        stdvec.watson.u2 <- stdvec.offset.step$watson.u2
        
        if ( search.h0.probe )
        {
            stdvec.h0.probe <- stdvec.offset.step$h0.probe
            stdvec.h0.probe.num <- length( stdvec.h0.probe )
        }
        
        # check errors
        if ( any( is.nan( stdvec.offset.delta ) ) || 
                is.nan( stdvec.offset.stderr ) )
            stop( "NaN error in normalize.stdvec.condec" )
        
        if ( stdvec.numerical.demand < .Machine$double.eps * 10^3 )
            stop( "Numerical error in normalize.stdvec.condec" )
        
        # update total standard vector offset
        stdvec.offset <- stdvec.offset + stdvec.offset.delta
        stdvec.offset <- stdvec.offset - mean( stdvec.offset )
        
        # update data at once with total offset
        last.norm.data <- sweep( edata, 2, stdvec.offset )
        if ( search.h0.probe )
            last.norm.data.fstat <- sweep( edata.fstat, 2, stdvec.offset )
        
        # check convergence
        stdvec.offset.sd <- sd( stdvec.offset )
        stdvec.offset.delta.sd <- sd( stdvec.offset.delta )
        
        stdvec.offset.stderr.ratio <- ifelse( stdvec.offset.sd > 0, 
            stdvec.offset.stderr / stdvec.offset.sd, 1 )
        stdvec.offset.delta.sd.ratio <- ifelse( stdvec.offset.sd > 0, 
            stdvec.offset.delta.sd / stdvec.offset.sd, 1 )
        
        if ( stdvec.offset.delta.sd == 0 )
            stdvec.offset.ratio <- 0
        else if ( stdvec.offset.stderr == 0 )
            stdvec.offset.ratio <- 1
        else
            stdvec.offset.ratio <- stdvec.offset.delta.sd / stdvec.offset.stderr
        
        stdvec.offset.accum.step <- ifelse( 
            stdvec.offset.ratio < accum.threshold, 
            stdvec.offset.accum.step + 1, 0 )
        
        stdvec.convergence <- stdvec.offset.ratio < single.threshold || 
            stdvec.offset.accum.step > offset.accum.step.threshold
        
        if ( ! search.h0.probe )
            overall.convergence <- stdvec.convergence
        else
        {
            stdvec.common.h0.probe.accum.step <- ifelse( stdvec.convergence, 
                stdvec.common.h0.probe.accum.step + 1, 0 )
            
            overall.convergence <- stdvec.common.h0.probe.accum.step >= 
                common.h0.probe.accum.step.max
        }
        
        # store last results
        last.stdvec.offset <- stdvec.offset
        
        if ( search.h0.probe )
        {
            if ( stdvec.common.h0.probe.accum.step == 0 )
                stdvec.common.h0.probe <- NULL
            else if ( stdvec.common.h0.probe.accum.step == 1 )
                stdvec.common.h0.probe <- stdvec.h0.probe
            else
                stdvec.common.h0.probe <- intersect( stdvec.h0.probe, 
                    stdvec.common.h0.probe )
            
            stdvec.common.h0.probe.num <- ifelse( 
                is.null( stdvec.common.h0.probe ), NA, 
                length( stdvec.common.h0.probe ) )
        }
        
        # store step results
        norm.stdvec.offset <- rbind( norm.stdvec.offset, stdvec.offset )
        norm.stdvec.offset.sd <- c( norm.stdvec.offset.sd, stdvec.offset.sd )
        norm.stdvec.offset.stderr <- c( norm.stdvec.offset.stderr, 
            stdvec.offset.stderr )
        norm.stdvec.offset.delta.sd <- c( norm.stdvec.offset.delta.sd, 
            stdvec.offset.delta.sd )
        norm.stdvec.offset.accum.step <- c( norm.stdvec.offset.accum.step, 
            stdvec.offset.accum.step )
        norm.stdvec.numerical.demand <- c( norm.stdvec.numerical.demand, 
            stdvec.numerical.demand )
        norm.stdvec.watson.u2 <- rbind( norm.stdvec.watson.u2, 
            stdvec.watson.u2 )
        
        if ( search.h0.probe )
        {
            norm.stdvec.common.h0.probe.accum.step <- 
                c( norm.stdvec.common.h0.probe.accum.step, 
                    stdvec.common.h0.probe.accum.step )
            norm.stdvec.h0.probe.num <- c( norm.stdvec.h0.probe.num, 
                stdvec.h0.probe.num )
            norm.stdvec.common.h0.probe.num <- 
                c( norm.stdvec.common.h0.probe.num, stdvec.common.h0.probe.num )
        }
        
        if ( verbose )
        {
            if ( ! is.null( stdvec.watson.u2 ) )
                stdvec.watson.u2.char <- paste0( signif( stdvec.watson.u2, 6 ), 
                    collapse=" " )
            else
                stdvec.watson.u2.char <- ""
            
            cat( sprintf( "  %2d %g %g %g %02d %g [%s]", iter, 
                stdvec.offset.sd, stdvec.offset.stderr.ratio, 
                stdvec.offset.delta.sd.ratio, stdvec.offset.accum.step, 
                stdvec.numerical.demand, stdvec.watson.u2.char ) )
            
            if ( search.h0.probe )
            {
                cat( sprintf( " %02d %d", stdvec.common.h0.probe.accum.step, 
                    stdvec.h0.probe.num ) )
                
                if ( ! is.na( stdvec.common.h0.probe.num ) )
                    cat( paste0( " ", stdvec.common.h0.probe.num ) )
            }
            
            cat( "\n" )
        }
    }
    
    if ( verbose )
        cat( "\n" )
    
    if ( ! overall.convergence )
        stop( "No convergence in normalize.stdvec.condec" )
    
    # remove sample or condition names from step results
    dimnames( norm.stdvec.offset ) <- NULL
    if ( ! is.null( norm.stdvec.watson.u2 ) )
        dimnames( norm.stdvec.watson.u2 ) <- NULL
    
    norm.stdvec.convergence <- list( offset = norm.stdvec.offset, 
        offset.sd = norm.stdvec.offset.sd, 
        offset.stderr = norm.stdvec.offset.stderr, 
        offset.delta.sd = norm.stdvec.offset.delta.sd, 
        offset.accum.step = norm.stdvec.offset.accum.step, 
        numerical.demand = norm.stdvec.numerical.demand, 
        watson.u2 = norm.stdvec.watson.u2 )
    
    norm.stdvec.result <- list( offset = last.stdvec.offset, 
        convergence = norm.stdvec.convergence )
    
    if ( search.h0.probe )
    {
        norm.stdvec.result$h0.probe <- stdvec.common.h0.probe
        norm.stdvec.result$convergence$common.h0.probe.accum.step <- 
            norm.stdvec.common.h0.probe.accum.step
        norm.stdvec.result$convergence$h0.probe.num <- norm.stdvec.h0.probe.num
        norm.stdvec.result$convergence$common.h0.probe.num <- 
            norm.stdvec.common.h0.probe.num
    }
    
    norm.stdvec.result
}


calculate.stdvec.offset <- function( edata, condition, search.h0.probe, 
    edata.fstat, within.cond.var, within.cond.n, vector.graph, 
    vector.graph.probe, p.value.graph, iter )
{
    stdvec.trim <- 0.01
    ks.test.alpha <- 1e-3
    
    if ( ! search.h0.probe )
    {
        # identify probes for normalization
        expr.var <- apply( edata, 1, var )
        expr.var[ expr.var < max( expr.var ) * .Machine$double.eps  ] <- NA
        
        stdvec.probe <- names( which( 
            expr.var > quantile( expr.var, stdvec.trim/2, na.rm=TRUE ) & 
            expr.var < quantile( expr.var, 1 - stdvec.trim/2, na.rm=TRUE ) ) )
    }
    else
    {
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
        
        # identify probes for normalization
        expr.var <- apply( edata.fstat[ h0.probe, ], 1, var )
        expr.var[ expr.var < max( expr.var ) * .Machine$double.eps  ] <- NA
        
        stdvec.probe <- h0.probe[ 
            expr.var > quantile( expr.var, stdvec.trim/2, na.rm=TRUE ) & 
            expr.var < quantile( expr.var, 1 - stdvec.trim/2, na.rm=TRUE ) ]
    }
    
    # center and scale expression data
    expr.mean <- rowMeans( edata[ stdvec.probe, ] )
    expr.centered <- sweep( edata[ stdvec.probe, ], 1, expr.mean, "-" )
    
    expr.sd.inv <- 1 / apply( expr.centered, 1, sd )
    expr.scaled <- sweep( expr.centered, 1, expr.sd.inv, "*" )
    
    # calculate offset
    expr.sd.inv.sum <- sum( expr.sd.inv )
    stdvec.offset <- apply( expr.scaled, 2, sum ) / expr.sd.inv.sum
    stdvec.offset <- stdvec.offset - mean( stdvec.offset )
    
    # estimate error and numerical demand
    stdvec.offset.stderr <- sqrt( length( stdvec.probe ) ) / expr.sd.inv.sum
    expr.sd.inv.min <- min( expr.sd.inv )
    stdvec.numerical.demand <- expr.sd.inv.min / expr.sd.inv.sum
    
    # calculate and plot density distribution of standard vector angles
    
    dimension.num <- ncol( expr.scaled )
    
    if ( dimension.num < 3 )
        theta.watson.u2 <- NULL
    else
    {
        # identify condition groups
        dimension.group.num <- ( dimension.num - 1 ) %/% 3 + 1
        
        dimension.group <- lapply( 1 : dimension.group.num, function ( g )
            if ( g < dimension.group.num )
                ( 3*g - 2 ) : ( 3*g )
            else
                ( dimension.num - 2 ) : dimension.num )
        
        theta.watson.u2 <- vector( "numeric", dimension.group.num )
        
        uv <- matrix( c( 0, -1/sqrt(2), 1/sqrt(2), 
            2/sqrt(6), -1/sqrt(6), -1/sqrt(6) ), nrow=3 )
        
        for ( dim.group.idx in 1:dimension.group.num )
        {
            # select expression values for each condition group
            expr.dim.group <- expr.scaled[ , 
                dimension.group[[ dim.group.idx ]] ]
            
            if ( dimension.group.num > 1 )
            {
                # re-standardize again for this group
                expr.dim.group.mean <- rowMeans( expr.dim.group )
                expr.dim.group.centered <- sweep( expr.dim.group, 1, 
                    expr.dim.group.mean, "-" )
                
                expr.dim.group.sd <- apply( expr.dim.group.centered, 1, sd )
                expr.dim.group.sel <- expr.dim.group.sd != 0
                expr.dim.group <- sweep( 
                    expr.dim.group.centered[ expr.dim.group.sel, ], 1, 
                    expr.dim.group.sd[ expr.dim.group.sel ], "/" )
            }
            
            expr.uv <- expr.dim.group %*% uv
            expr.u <- expr.uv[ , 1 ]
            expr.v <- expr.uv[ , 2 ]
            
            # calculate density distribution of angles
            theta.density.n <- 2^11
            theta.density.adjust <- 0.5
            
            expr.theta <- atan2( expr.v, expr.u )
            expr.theta.ex <- c( expr.theta, 
                expr.theta + ifelse( expr.theta > 0, -2*pi, 2*pi ) )
            
            expr.theta.density <- density( expr.theta.ex, 
                adjust=theta.density.adjust, n=theta.density.n )
            expr.theta.density.sel <- expr.theta.density$x > -pi & 
                expr.theta.density$x <= pi
            
            # calculate density distribution of angles after permutations
            expr.theta.permu <- cbind( expr.theta.ex, 
                expr.theta.ex + (2*pi)/3, 
                expr.theta.ex - (2*pi)/3, 
                - expr.theta.ex + pi, 
                - expr.theta.ex + pi/3, 
                - expr.theta.ex - pi/3 )
            
            invar.theta <- expr.theta.permu[ expr.theta.permu > -pi &
                expr.theta.permu <= pi ]
            invar.theta.ex <- c( invar.theta, 
                invar.theta + ifelse( invar.theta > 0, -2*pi, 2*pi ) )
            
            invar.theta.density <- density( invar.theta.ex, 
                adjust=theta.density.adjust, n=theta.density.n )
            invar.theta.density.sel <- invar.theta.density$x > -pi & 
                invar.theta.density$x <= pi
            
            # calculate Watson U2 statistic
            theta.watson.u2[ dim.group.idx ] <- 
                watson.u2( expr.theta, sample( invar.theta, 
                    length( expr.theta ), replace=TRUE ) )
            
            if ( ! is.null( vector.graph ) && require( plotrix, quietly=TRUE ) )
            {
                if ( vector.graph != "" )
                {
                    if ( ! file.exists( vector.graph ) )
                        dir.create( vector.graph, recursive=TRUE )
                    
                    png.filename <- sprintf( 
                        "%s/stdvec_%s%s_iter%02d.png", vector.graph, condition, 
                        ifelse( dimension.group.num > 1, 
                            sprintf( "_dg%02d", dim.group.idx ), "" ), 
                        iter )
                    
                    png( png.filename, width=1280, height=720 )
                }
                
                par.default <- par( no.readonly=TRUE )
                
                par( mfrow = c(1,2), pty="s", xpd=FALSE, cex.lab=2, 
                    mar = c( 1, 0.5, 2.2, 0.5 ), oma = c( 0, 0, 0, 0 ) )
                
                # select offset values for condition group
                stdvec.offset.uv <- stdvec.offset[ 
                    dimension.group[[ dim.group.idx ]] ] %*% uv
                stdvec.offset.u <- stdvec.offset.uv[ , 1 ]
                stdvec.offset.v <- stdvec.offset.uv[ , 2 ]
                
                # plot a sample of standardized sample vectors
                uv.lim <- c( -1.5, 1.5 )
                plot( 0, type="n", xlim=uv.lim, ylim=uv.lim, axes=FALSE, 
                    ann=FALSE, frame.plot=FALSE, asp=1 )
                
                expr.uv.probe <- names( expr.u )
                if ( length( expr.uv.probe ) > length( vector.graph.probe ) ) {
                    expr.uv.probe <- intersect( expr.uv.probe, 
                        vector.graph.probe )
                }
                expr.uv.probe.num <- length( expr.uv.probe )
                
                expr.uv.factor <- 1.06
                expr.uv.color <- gray( 0.3 )
                expr.uv.width <- ifelse( expr.uv.probe.num > 1000, 0.1, 0.2 )
                segments( 0, 0, expr.uv.factor * expr.u[ expr.uv.probe ], 
                    expr.uv.factor * expr.v[ expr.uv.probe ], 
                    lwd=expr.uv.width, col=expr.uv.color )
                
                grid.pos.x <- c( 0, -sqrt(3/4), sqrt(3/4) )
                grid.pos.y <- c( 1, -1/2, -1/2 )
                
                grid.length <- expr.uv.factor * sqrt( 2 )
                segments( 0, 0, grid.length * grid.pos.x, 
                    grid.length * grid.pos.y, lwd=2.5, col="black" )
                
                stdvec.offset.uv.factor <- 10
                stdvec.offset.uv.color <- "red"
                segments( 0, 0, stdvec.offset.uv.factor * stdvec.offset.u, 
                    stdvec.offset.uv.factor * stdvec.offset.v, lwd=3, 
                    col="red" )
                
                par( xpd=TRUE )
                
                grid.label.length <- 1.63
                grid.labels <- c( "s1", "s2", "s3" )
                text( grid.label.length * grid.pos.x, 
                    grid.label.length * grid.pos.y, grid.labels, cex=2.5 )
                
                par( xpd=FALSE )
                
                stdvec.offset.mag <- sqrt( sum( stdvec.offset[ 
                    dimension.group[[ dim.group.idx ]] ]^2 ) )
                mtext( paste0( "||offset|| = ", sprintf( "%.3e", 
                    stdvec.offset.mag ) ), side=1, line=0.7, cex=2.5 )
                
                # plot polar distributions of standard vector angles
                polar.expr.theta <- expr.theta.density$x[ 
                    expr.theta.density.sel ]
                polar.expr.rho <- 2 * expr.theta.density$y[ 
                    expr.theta.density.sel ]

                polar.invar.theta <- invar.theta.density$x[ 
                    invar.theta.density.sel ]
                polar.invar.rho <- 2 * invar.theta.density$y[ 
                    invar.theta.density.sel ]
                
                theta.labels <- c( "", "", "" )
                rho.grid <- seq( 0, 3/(2*pi), length.out=4 )
                rho.labels <- c( "", "", expression(1/pi), "" )
                
                radial.plot( polar.expr.rho, polar.expr.theta - pi/2, 
                    rp.type="p", start=pi/2, radial.lim=rho.grid, 
                    show.grid.labels = length( theta.labels ), 
                    labels=theta.labels, radial.labels=rho.labels, lwd=2.5, 
                    mar = par( "mar" ) )
                
                radial.plot( polar.invar.rho, polar.invar.theta - pi/2, 
                    rp.type="p", start=pi/2, radial.lim=rho.grid, lwd=2, lty=2, 
                    line.col="blue", add=TRUE )
                
                grid.label.length <- 0.52
                text( grid.label.length * grid.pos.x, 
                    grid.label.length * grid.pos.y, grid.labels, cex=2.5 )
                
                mtext( substitute( paste( "Watson U"^"2", " = ", wu2 ), 
                    list( wu2 = sprintf( "%.3e", 
                        theta.watson.u2[ dim.group.idx ] ) ) ), 
                    side=1, line=0.8, cex=2.5 )
                
                graph.title <- sprintf( "%s%s    -    iter=%02d", condition, 
                    ifelse( dimension.group.num > 1, 
                        sprintf( ":%02d", dim.group.idx ), "" ), 
                    iter )
                title( main=graph.title, outer=TRUE, line=-3, font.main=1, 
                    cex.main=2.7 )
                
                if ( vector.graph != "" )
                    dev.off()
                else
                    par( par.default )
            }
        }
    }
    
    stdvec.result <- list( value = stdvec.offset, 
        stderr = stdvec.offset.stderr, 
        numerical.demand = stdvec.numerical.demand, 
        watson.u2 = theta.watson.u2 )
    
    if ( search.h0.probe )
        stdvec.result$h0.probe = h0.probe
    
    stdvec.result
}

