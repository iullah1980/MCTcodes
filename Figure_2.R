require(tidyr)
require(ggplot2)


n1 <- 6
n2 <- 3
mu1 <- 1
mu2 <- 1
sd1 <- 1
sd2 <- 1.2 #sqrt(sd1^2*(n1+n2*lambda)/(n1*lambda))

random.seed <- 23452
set.seed( random.seed )

res <- replicate(1000,{
  
  lambda <- (sd1^2/n1)/(sd1^2/n1+sd2^2/n2)
  tt <-  rnorm(100000)/sqrt(lambda*rchisq(100000,n1-1)/(n1-1)+(1-lambda)*rchisq(100000,n2-1)/(n2-1))
  d1 <- density(tt,from=-5,to=5,n=1000)
  
  y1 <- rnorm(n=n1, mean=mu1, sd=sd1)
  y2 <- rnorm(n=n2, mean=mu2, sd=sd2)
  y1.bar <- mean(y1)
  y2.bar <- mean(y2)
  var1 <- var(y1)
  var2 <- var(y2)
  
  Tn <- ( y1.bar-y2.bar )/sqrt(var1/n1 + var2/n2)
  f <- (var1/n1 + var2/n2)^2 / (var1^2/(n1^3-n1^2) + var2^2/(n2^3-n2^2))
  
  tw <- rt(100000,f)
  d2 <- density(tw,from=-5,to=5,n=1000)
  
  alpha <- (var1/n1)/(var1/n1+var2/n2)
  ty <-  rnorm(100000)/sqrt(alpha*rchisq(100000,n1-1)/(n1-1)+(1-alpha)*rchisq(100000,n2-1)/(n2-1))
  d3 <- density(ty,from=-5,to=5,n=1000)
  
  rbind(d1$x,d1$y,d2$y,d3$y)
  
})

aa <- apply(res,c(1,2),mean)
dat <- data.frame(t(aa))
dd <-  gather(dat, lab, pp,-X1)


fulltheme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   panel.border = element_rect(colour = "black", fill=NA, size=1),
                   plot.margin = unit(c(0,1,5,5),"mm"),
                   axis.text.x =  element_text(size = 15 , lineheight = 2, colour = "black", vjust = 1),
                   axis.text.y =  element_text(size = 15, lineheight = 2, colour = "black", hjust = 1),
                   axis.title = element_text(size = 20)
)

zoomtheme <- theme(legend.position="none", 
                   axis.title.x=element_blank(),axis.title.y=element_blank(),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                   panel.background = element_rect(color='black', fill="white"),
                   plot.margin = unit(c(0,0,0,0),"mm"))


xlim <- c(-5,-3.5) ; ylim <- c(0.002061365,0.010957709)
p.zoom <- ggplot(data=dat, aes_string(x = "X1")) + 
  geom_line(aes_string(y = "X2"),col=1) +
  geom_line(aes_string(y = "X3"),col=2) +
  geom_line(aes_string(y = "X4"),col=3) +
  coord_cartesian(xlim=xlim, ylim=ylim) + zoomtheme
# put them together
g1 <- ggplotGrob(p.zoom)

xlim <- c(3.5,5) ; ylim <- c(0.002061365,0.010957709)
p.zoom <- ggplot(data=dat, aes_string(x = "X1")) + 
  geom_line(aes_string(y = "X2"),col=1) +
  geom_line(aes_string(y = "X3"),col=2) +
  geom_line(aes_string(y = "X4"),col=3) +
  coord_cartesian(xlim=xlim, ylim=ylim) + zoomtheme
# put them together
g2 <- ggplotGrob(p.zoom)


p1 <- ggplot(data=dd, aes(x=X1,y=pp,color=lab))+
  geom_line()+fulltheme+
  labs(title = "") + xlab("T") + ylab("density\n")+ 
  theme(legend.position=c(.10,.9),legend.background = element_rect(color= "black"),
        legend.key = element_rect(colour = NA, fill = 'white'))+
  scale_colour_manual(name="",values=c("black","red","green"),labels=c("TRUE", "Welch's Approximation", "MC Approximation"), 
                      guide = guide_legend(override.aes=aes(fill=NA)))+xlim(-5,5)

postscript(file="Fig2.eps")
p1 + annotation_custom(grob = g1, xmin = -5.3, xmax = -2, ymin = 0.08, ymax = 0.19)+
  annotation_custom(grob = g2, xmin = 2, xmax = 5.3, ymin = 0.08, ymax = 0.19)
dev.off()

