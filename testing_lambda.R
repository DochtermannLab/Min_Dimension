library(viridisLite)
#load data from Dochtermann et al. 2023 for simulated G matrices
load("C:/Users/neddo/Dropbox/Working/Projects/HoleyLandscapes/ManuscriptMaterials (Code, etc)/IVCsimulation_6combs.RData")

#convert correlation matrix to covariance matrix
cor2cov<-function(vars,cormat){   
  sdMat<-diag(sqrt(vars))
  corMat<-cormat
  mat<-sdMat %*% corMat %*% t(sdMat)
  return(mat)
}

#### Distribution of min, max, and mean correlations via onion method ####
#LKJ onion method for generating random correlation matrices,
#third-hand from McElreath's "rethinking" package
rlkjcorr <- function (n, K, eta = 1) {
  
  stopifnot(is.numeric(K), K >= 2, K == as.integer(K))
  stopifnot(eta > 0)
  #if (K == 1) return(matrix(1, 1, 1))
  
  f <- function() {
    alpha <- eta + (K - 2)/2
    r12 <- 2 * rbeta(1, alpha, alpha) - 1
    R <- matrix(0, K, K) # upper triangular Cholesky factor until return()
    R[1,1] <- 1
    R[1,2] <- r12
    R[2,2] <- sqrt(1 - r12^2)
    if(K > 2) for (m in 2:(K - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m / 2, alpha)
      
      # Draw uniformally on a hypersphere
      z <- rnorm(m, 0, 1)
      z <- z / sqrt(crossprod(z)[1])
      
      R[1:m,m+1] <- sqrt(y) * z
      R[m+1,m+1] <- sqrt(1 - y)
    }
    return(crossprod(R))
  }
  R <- replicate( n , f() )
  if ( dim(R)[3]==1 ) {
    R <- R[,,1]
  } else {
    # need to move 3rd dimension to front, so conforms to array structure that Stan uses
    R <- aperm(R,c(3,1,2))
  }
  return(R)
}

#### Generate random matrices and calculate relevant metrics
sims <- 5000
ks <- c(3,5,10,15,20)
stor <- matrix(NA,nrow = sims*length(ks),ncol=3)
counter <- 0
for(j in ks){
for(i in 1:sims)
{
  counter <- counter+1
  x <- rlkjcorr(1,j,eta=1)
  vars <- rgamma(j,10)
  x_cov <- cor2cov(vars,x)
  stor[counter,1] <- j
  stor[counter,2] <- min(eigen(x_cov)$values)
  stor[counter,3] <- max(eigen(x_cov)$values)
  
}}
stor <- as.data.frame(stor)
names(stor) <- c("k","min","max")
stor$ratio <- stor$min/stor$max

#means and sds
ratio_summ <- cbind(ks,aggregate(stor$ratio,by=list(stor$k),mean)[,2],
                    (aggregate(stor$ratio,by=list(stor$k),mean)-aggregate(stor$ratio,by=list(stor$k),sd))[,2],
      (aggregate(stor$ratio,by=list(stor$k),mean)+aggregate(stor$ratio,by=list(stor$k),sd))[,2])


k_3 <- stor[which(stor$k==3),]
k_5 <- stor[which(stor$k==5),]
k_10 <- stor[which(stor$k==10),]
k_15 <- stor[which(stor$k==15),]
k_20 <- stor[which(stor$k==20),]

median(stor$ratio)
max(stor$ratio)

which(stor$ratio==max(stor$ratio)) #if less than 5000, k = 3

ratio_summ

cut_offs <- cbind(c(3,5,10,15,20),rbind(quantile(k_3$ratio,.05),
quantile(k_5$ratio,.05),
quantile(k_10$ratio,.05),
quantile(k_15$ratio,.05),
quantile(k_20$ratio,.05)))

cut_offs

#### Generate disibtrution of corrs and vars
stor1 <- matrix(NA,100000,3)
for(i in 1:100000)
{
  x <- rlkjcorr(1,10,eta=1)
  low_tri <- x[which(lower.tri(x)==TRUE)]
  stor1[i,1] <- min(abs(low_tri))
  stor1[i,2] <- max(abs(low_tri))
  stor1[i,3] <- mean(abs(low_tri))
}

r_vars <- rgamma(100000,10)

Wright_ratio <- eigen.stor.Wright.1[,10]/eigen.stor.Wright.1[,1]
sum(Wright_ratio>0.0047213272)

#### FIGURES: A: random correlations, B: random variances,  ####
#### C: observed values, D: simulated populations
#pal <- viridis(6)
pal <- c("#440154FF", "#414487FF", "#2A788EFF",
         "#22A884FF" ,"#7AD151FF" ,"#FDE725FF")

par(mfrow=c(2,2),mar=c(5, 6.5, 4, 2) + 0.1,
    oma = c(0, 0, 0, 0),pty='m')

#A
plot(x=density(stor1[,1])$x, 
     y=density(stor1[,1])$y/(max(density(stor1[,1])$y)),
     xlim=c(0,1),type='l',xlab="|r|",ylab="Density (standardized)",
     cex.axis=1.25,cex.lab=1.5)
polygon(x=density(stor1[,1])$x,y=density(stor1[,1])$y/(max(density(stor1[,1])$y)),
        col=adjustcolor(pal[1],alpha.f = .5))
abline(v=median(stor1[,1]),lty=2)

lines(x=density(stor1[,2])$x,y=density(stor1[,2])$y/(max(density(stor1[,2])$y)))
polygon(x=density(stor1[,2])$x,y=density(stor1[,2])$y/(max(density(stor1[,2])$y)),
        col=adjustcolor(pal[4],alpha.f = .5))
abline(v=median(stor1[,2]),lty=2)

lines(x=density(stor1[,3])$x,y=density(stor1[,3])$y/(max(density(stor1[,3])$y)))
polygon(x=density(stor1[,3])$x,y=density(stor1[,3])$y/(max(density(stor1[,3])$y)),
        col=adjustcolor(pal[6],alpha.f = .5))
abline(v=median(stor1[,3]),lty=2)
mtext("A",side=3,line=-.25,cex=1.25,adj=-.3)

#B
plot(x=density(r_vars)$x,y=density(r_vars)$y/(max(density(r_vars)$y)),
     type='l',xlab="Variance",ylab="Density (standardized)",
     cex.axis=1.25,cex.lab=1.5)
polygon(x=density(r_vars)$x,y=density(r_vars)$y/(max(density(r_vars)$y)),
        col=adjustcolor(pal[2],alpha.f = .5))
abline(v=median(r_vars),lty=2) 
mtext("B",side=3,line=-.25,cex=1.25,adj=-.3)

#C
plot(NA,xlim=c(0,20),ylim=c(0,max(ratio_summ[,4])),
     xlab="Dimensions (k)",
     ylab=expression(paste(frac(lambda[min],lambda[max]))),
     cex.axis=1.25,cex.lab=1.5)
lines(x=ratio_summ[,1],y=ratio_summ[,2])
arrows(x0=ratio_summ[,1],y0=ratio_summ[,3],
       x1=ratio_summ[,1],y1=ratio_summ[,4],
       code=3,angle=90,length=.1)
points(x=ratio_summ[,1],y=ratio_summ[,2],pch=21,bg=pal[3],cex=3)
mtext("C",side=3,line=-.25,cex=1.25,adj=-.3)

#D
plot(x=density(Wright_ratio)$x,
     y=density(Wright_ratio)$y/(max(density(Wright_ratio)$y)),
     type='l',ylab="Density (standardized)",
     xlab="",
     cex.axis=1.25,cex.lab=1.5)
polygon(x=density(Wright_ratio)$x,
        y=density(Wright_ratio)$y/(max(density(Wright_ratio)$y)),
        col=adjustcolor(pal[5],alpha.f = .5))
abline(v=cut_offs[3,2],lty=2) 
mtext(expression(paste(frac(lambda[min],lambda[max]))),side=1,line=4)
mtext("D",side=3,line=-.25,cex=1.25,adj=-.3)
