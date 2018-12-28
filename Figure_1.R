
BF.sim <- function(n1, n2, mu1, mu2, sd1, sd2, size, n.replicate){
  v.sim <- replicate(n.replicate,{
    y1 <- rnorm(n=n1, mean=mu1, sd=sd1)
    y2 <- rnorm(n=n2, mean=mu2, sd=sd2)
    
    y1.bar <- mean(y1)
    y2.bar <- mean(y2)
    
    var1 <- var(y1)
    var2 <- var(y2)
    
    Tn <- ( y1.bar-y2.bar )/sqrt(var1/n1 + var2/n2)
    
    f <- (var1/n1 + var2/n2)^2 / (var1^2/(n1^3-n1^2) + var2^2/(n2^3-n2^2))
    pT1 <- 2*pt(abs(Tn), f, lower.tail = FALSE)
    
    alpha <- (var1/n1)/(var1/n1+var2/n2)
    TNull <- rnorm(10000)/sqrt(alpha*rchisq(10000,n1-1)/(n1-1)+(1-alpha)*rchisq(10000,n2-1)/(n2-1))
    pT <- mean(abs(TNull) >= abs(Tn))

    c(pT1,pT)   
  }
  )
  
  T1p <- mean(v.sim[1,] < size)
  
  TP <- mean(v.sim[2,] < size)
  
  c(T1p,TP)
}

random.seed <- 23452
set.seed( random.seed )

n1 <- 10
n2 <- 4

mu1 <- 1

sd1 <- 1
sd2 <- 2

size <- 0.05
delta <- c(0,1,2,3)
tau <- delta*sqrt(sd1^2/n1+sd2^2/n2)
n.replicate <- 10000

M1 <- array(NA, c(2, length(tau)))

for (i in 1:length(tau)){
  mu2 <- mu1+tau[i]
  pow.v1 <- BF.sim(n1, n2, mu1, mu2, sd1, sd2,size, n.replicate)
  M1[,i] <- pow.v1
}

#######################################################
n1 <- 4
n2 <- 10

mu1 <- 1

sd1 <- 1
sd2 <- 2

size <- 0.05
delta <- c(0,1,2,3)
tau <- delta*sqrt(sd1^2/n1+sd2^2/n2)
n.replicate <- 10000

M2 <- array(NA, c(2, length(tau)))

for (i in 1:length(tau)){
  mu2 <- mu1+tau[i]
  pow.v1 <- BF.sim(n1, n2, mu1, mu2, sd1, sd2,size, n.replicate)
  M2[,i] <- pow.v1
}
############################################################
n1 <- 10
n2 <- 4

mu1 <- 1

sd1 <- 1
sd2 <- 4

size <- 0.05
delta <- c(0,1,2,3)
tau <- delta*sqrt(sd1^2/n1+sd2^2/n2)
n.replicate <- 10000

M3 <- array(NA, c(2, length(tau)))

for (i in 1:length(tau)){
  mu2 <- mu1+tau[i]
  pow.v1 <- BF.sim(n1, n2, mu1, mu2, sd1, sd2,size, n.replicate)
  M3[,i] <- pow.v1
}
###############################################################
n1 <- 4
n2 <- 10

mu1 <- 1

sd1 <- 1
sd2 <- 4

size <- 0.05
delta <- c(0,1,2,3)
tau <- delta*sqrt(sd1^2/n1+sd2^2/n2)
n.replicate <- 10000

M4 <- array(NA, c(2, length(tau)))

for (i in 1:length(tau)){
  mu2 <- mu1+tau[i]
  pow.v1 <- BF.sim(n1, n2, mu1, mu2, sd1, sd2,size, n.replicate)
  M4[,i] <- pow.v1
}

################################################################
################################################################
pdf(file="Fig1.pdf",width=14, height=14)
par(mfrow=c(2,2), mar=c(6, 6, 4, 2))
n1 <- 10
n2 <- 4
plot(delta, M1[1 ,]
     ,ylim=c(0,.8),  ylab="", xaxt='n'
     , xlab=""
     ,xlim=c(0, 3),type ="l",col=2,lwd=2, cex.axis=2
)
axis(1,at=c(0, 1, 2, 3),labels=c( 0, 1, 2, 3),line=0,tck=0, cex.axis=2)

lines(delta, M1[2,] ,col=4,lwd=2)
abline(h=0.05, col="gray",lwd=2)

mtext("(a)",at=0, line = 1, cex =2.5)
title(main=bquote(paste("(",n[1],", ",n[2], ")", " = ","(",.(n1),", ",.(n2),")",", ",rho," = ","1/4" )),cex.main=2)
title(xlab=expression(eta), cex.lab=2.5, line = 3.5)
title(ylab="power", cex.lab=2.5, line=3.5)
legend("topleft", legend=c("W", "MCT"), 
       lty=c(1,1), col=c(2,4), 
       cex = 2,lwd=c(2,2))

#############################################
n1 <- 4
n2 <- 10
plot(delta, M2[1 ,]
     ,ylim=c(0,.8),  ylab="", xaxt='n'
     , xlab=""
     ,xlim=c(0, 3),type ="l",col=2,lwd=2, cex.axis=2
)
axis(1,at=c(0, 1, 2, 3),labels=c( 0, 1, 2, 3),line=0,tck=0, cex.axis=2)

lines(delta, M2[2,] ,col=4,lwd=2)
abline(h=0.05, col="gray",lwd=2)

mtext("(b)",at=0, line = 1, cex =2.5)
title(main=bquote(paste("(",n[1],", ",n[2], ")", " = ","(",.(n1),", ",.(n2),")",", ",rho," = ","1/4" )),cex.main=2)
title(xlab=expression(eta), cex.lab=2.5, line = 3.5)
title(ylab="power", cex.lab=2.5, line=3.5)
legend("topleft", legend=c("W", "MCT"), 
       lty=c(1,1), col=c(2,4), 
       cex = 2,lwd=c(2,2))
#############################################
n1 <- 10
n2 <- 4
plot(delta, M3[1 ,]
     ,ylim=c(0,.8),  ylab="", xaxt='n'
     , xlab=""
     ,xlim=c(0, 3),type ="l",col=2,lwd=2, cex.axis=2
)
axis(1,at=c(0, 1, 2, 3),labels=c( 0, 1, 2, 3),line=0,tck=0, cex.axis=2)

lines(delta, M3[2,] ,col=4,lwd=2)
abline(h=0.05, col="gray",lwd=2)

mtext("(c)",at=0, line = 1, cex =2.5)
title(main=bquote(paste("(",n[1],", ",n[2], ")", " = ","(",.(n1),", ",.(n2),")",", ",rho," = ","1/16" )),cex.main=2)
title(xlab=expression(eta), cex.lab=2.5, line = 3.5)
title(ylab="power", cex.lab=2.5, line=3.5)
legend("topleft", legend=c("W", "MCT"), 
       lty=c(1,1), col=c(2,4), 
       cex = 2,lwd=c(2,2))
#############################################
n1 <- 4
n2 <- 10
plot(delta, M4[1 ,]
     ,ylim=c(0,.8),  ylab="", xaxt='n'
     , xlab=""
     ,xlim=c(0, 3),type ="l",col=2,lwd=2, cex.axis=2
)
axis(1,at=c(0, 1, 2, 3),labels=c( 0, 1, 2, 3),line=0,tck=0, cex.axis=2)

lines(delta, M4[2,] ,col=4,lwd=2)
abline(h=0.05, col="gray",lwd=2)

mtext("(d)",at=0, line = 1, cex =2.5)
title(main=bquote(paste("(",n[1],", ",n[2], ")", " = ","(",.(n1),", ",.(n2),")",", ",rho," = ","1/16" )),cex.main=2)
title(xlab=expression(eta), cex.lab=2.5, line = 3.5)
title(ylab="power", cex.lab=2.5, line=3.5)
legend("topleft", legend=c("W", "MCT"), 
       lty=c(1,1), col=c(2,4), 
       cex = 2,lwd=c(2,2))
dev.off()
##############################################



################################################################
################################################################
postscript(file="Fig1a.eps")
par(mar=c(8, 8, 4, 2))
n1 <- 10
n2 <- 4
plot(delta, M1[1 ,]
     ,ylim=c(0,.8),  ylab="", xaxt='n'
     , xlab=""
     ,xlim=c(0, 3),type ="l",col=2,lwd=2, cex.axis=2
)
axis(1,at=c(0, 1, 2, 3),labels=c( 0, 1, 2, 3),line=0,tck=0, cex.axis=2)

lines(delta, M1[2,] ,col=4,lwd=2)
abline(h=0.05, col="gray",lwd=2)

mtext("(a)",at=0, line = 1, cex =2.5)
title(main=bquote(paste("(",n[1],", ",n[2], ")", " = ","(",.(n1),", ",.(n2),")",", ",rho," = ","1/4" )),cex.main=2)
title(xlab=expression(eta), cex.lab=2.5, line = 3.5)
title(ylab="power", cex.lab=2.5, line=3.5)
legend("topleft", legend=c("W", "MCT"), 
       lty=c(1,1), col=c(2,4), 
       cex = 2,lwd=c(2,2))
dev.off()
#############################################
postscript(file="Fig1b.eps")
par(mar=c(8, 8, 4, 2))
n1 <- 4
n2 <- 10
plot(delta, M2[1 ,]
     ,ylim=c(0,.8),  ylab="", xaxt='n'
     , xlab=""
     ,xlim=c(0, 3),type ="l",col=2,lwd=2, cex.axis=2
)
axis(1,at=c(0, 1, 2, 3),labels=c( 0, 1, 2, 3),line=0,tck=0, cex.axis=2)

lines(delta, M2[2,] ,col=4,lwd=2)
abline(h=0.05, col="gray",lwd=2)

mtext("(b)",at=0, line = 1, cex =2.5)
title(main=bquote(paste("(",n[1],", ",n[2], ")", " = ","(",.(n1),", ",.(n2),")",", ",rho," = ","1/4" )),cex.main=2)
title(xlab=expression(eta), cex.lab=2.5, line = 3.5)
title(ylab="power", cex.lab=2.5, line=3.5)
legend("topleft", legend=c("W", "MCT"), 
       lty=c(1,1), col=c(2,4), 
       cex = 2,lwd=c(2,2))
dev.off()
#############################################

postscript(file="Fig1c.eps")
par(mar=c(8, 8, 4, 2))
n1 <- 10
n2 <- 4
plot(delta, M3[1 ,]
     ,ylim=c(0,.8),  ylab="", xaxt='n'
     , xlab=""
     ,xlim=c(0, 3),type ="l",col=2,lwd=2, cex.axis=2
)
axis(1,at=c(0, 1, 2, 3),labels=c( 0, 1, 2, 3),line=0,tck=0, cex.axis=2)

lines(delta, M3[2,] ,col=4,lwd=2)
abline(h=0.05, col="gray",lwd=2)

mtext("(c)",at=0, line = 1, cex =2.5)
title(main=bquote(paste("(",n[1],", ",n[2], ")", " = ","(",.(n1),", ",.(n2),")",", ",rho," = ","1/16" )),cex.main=2)
title(xlab=expression(eta), cex.lab=2.5, line = 3.5)
title(ylab="power", cex.lab=2.5, line=3.5)
legend("topleft", legend=c("W", "MCT"), 
       lty=c(1,1), col=c(2,4), 
       cex = 2,lwd=c(2,2))
dev.off()
#############################################

postscript(file="Fig1d.eps")
par(mar=c(8, 8, 4, 2))
n1 <- 4
n2 <- 10
plot(delta, M4[1 ,]
     ,ylim=c(0,.8),  ylab="", xaxt='n'
     , xlab=""
     ,xlim=c(0, 3),type ="l",col=2,lwd=2, cex.axis=2
)
axis(1,at=c(0, 1, 2, 3),labels=c( 0, 1, 2, 3),line=0,tck=0, cex.axis=2)

lines(delta, M4[2,] ,col=4,lwd=2)
abline(h=0.05, col="gray",lwd=2)

mtext("(d)",at=0, line = 1, cex =2.5)
title(main=bquote(paste("(",n[1],", ",n[2], ")", " = ","(",.(n1),", ",.(n2),")",", ",rho," = ","1/16" )),cex.main=2)
title(xlab=expression(eta), cex.lab=2.5, line = 3.5)
title(ylab="power", cex.lab=2.5, line=3.5)
legend("topleft", legend=c("W", "MCT"), 
       lty=c(1,1), col=c(2,4), 
       cex = 2,lwd=c(2,2))
dev.off()
##############################################

