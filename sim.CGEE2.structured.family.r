# Updated: 30 June, 2017
# ------------------------------------------------------------------------------
# Conditional second-order estimating equation for ascertained family data
# ------------------------------------------------------------------------------
# Methods proposed in Zhong and Cook (2017). Second-order estimating equations for the analysis of current status data from response-dependent biased family studies. Statistics in Biosciences.
# ------------------------------------------------------------------------------

# Event time is available for probands, non-probands only provide current status data
# ------------------------------------------------------------------------------
# Marginal distribution for event time: Weibull distribution
# Within-family association: Gaussian Copula
# ------------------------------------------------------------------------------
# Consider 2 generation family: 2 parents and 2 children.
# The pairwise relationships are: between parents, parent-child, and siblings
# ------------------------------------------------------------------------------

library(MASS)
library(splines)
library(survival)
library(copula)
library(mnormt)

# From Ktau to rho
Ktau.to.rho.f <- function(Ktau){
  rho <- sin((pi*Ktau)/2)
  return(rho)
}

# From Ktau.vec to alp.vec
Ktau.vec.to.alp.vec.f <- function(Ktau.vec){
  alp.vec <- log((1+Ktau.vec)/(1-Ktau.vec))
  return(alp.vec)
}

# From Ktau.vec to gam.vec (in Second order regression)
Ktau.vec.to.gam.vec.f <- function(Ktau.vec){
  alp.vec <- log((1+Ktau.vec)/(1-Ktau.vec))
  gam.vec <- alp.vec - c(0, alp.vec[1], alp.vec[1])
  return(gam.vec)
}

# From gam.vec to Ktau.vec
gam.vec.to.Ktau.vec.f <- function(gam.vec){
  alp.vec <- gam.vec + c(0, gam.vec[1], gam.vec[1])
  Ktau.vec <- (exp(alp.vec)-1)/(exp(alp.vec)+1)
  return(Ktau.vec)
}

# ----------------------------------------------------------
# Setup Parameters for structure gaussian
# ----------------------------------------------------------
getstats.str.f <- function(kap, beta, Ktau.pp, Ktau.ss, Ktau.ps, ad.cens, med.age) {
  lam <- ((-log(0.5))^(1/kap))/med.age
  tau.pp <- Ktau.pp
  tau.ps <- Ktau.ps
  tau.ss <- Ktau.ss
  rho.pp <- Ktau.to.rho.f(tau.pp)
  rho.ps <- Ktau.to.rho.f(tau.ps)
  rho.ss <- Ktau.to.rho.f(tau.ss)
  alp.vec <- Ktau.vec.to.alp.vec.f(c(tau.pp, tau.ss, tau.ps))
  gam.vec <- Ktau.vec.to.gam.vec.f(c(tau.pp, tau.ss, tau.ps))
  out <- NULL
  out$lam <- lam
  out$kap <- kap
  out$beta <- beta
  out$ad.cens <- ad.cens
  out$tau.pp <- tau.pp
  out$tau.ps <- tau.ps
  out$tau.ss <- tau.ss
  out$rho.pp <- rho.pp
  out$rho.ps <- rho.ps
  out$rho.ss <- rho.ss
  out$alp.vec <- alp.vec
  out$gam.vec <- gam.vec
  return(out)
}

# ----------------------------------------------------------
# Generate biased family data, structure Gaussian
# ----------------------------------------------------------
generate.str.f <- function(nsize, stats){
  n <- nsize
  lam <- stats$lam
  kap <- stats$kap
  beta <- matrix(stats$beta, ncol=1)
  ad.cens <- stats$ad.cens
  tau.pp <- stats$tau.pp
  tau.ps <- stats$tau.ps
  tau.ss <- stats$tau.ss
  rho.pp <- sin((pi*tau.pp)/2)
  rho.ss <- sin((pi*tau.ss)/2)
  rho.ps <- sin((pi*tau.ps)/2)
  mycop <- normalCopula(para=c(rho.pp, rep(rho.ps, 4), rho.ss), dim=4, dispstr = "un")

  clinic.entry <- rnorm(n, mean=50, sd=sqrt(20))
  ith <- 1
  outdata <- NULL
  while(ith <= n){
    proband <- sample(seq(1, 4), size=1, prob=rep(1/4, 4), replace=TRUE)
    nrej <- 0
    t1 <- 999999
    while(t1 > clinic.entry[ith]){
      uu <- as.numeric(rCopula(1, mycop))
      xx <- matrix(rbinom(4, size=1, prob=0.5), ncol=1)
      tt <- (((-log(1-uu))*exp((-1)*as.numeric(xx %*% beta)))^(1/kap))/lam
      exam.age <- c(rnorm(2, mean=60, sd=sqrt(10)), rnorm(2, mean=30, sd=sqrt(10)))
      exam.age <- ifelse(exam.age <= ad.cens, exam.age, ad.cens)
      exam.age[proband] <- clinic.entry[ith]
      data.tmp <- data.frame(fam.id=rep(ith, 4), rel.id=c(1, 1, 2, 2), event.time=tt, cens=exam.age, x1 = xx, proband=rep(0, 4))
      data.tmp$proband[proband] <- 1
      t1 <- data.tmp$event.time[proband]
      if(t1 < clinic.entry[ith]){
        data.tmp$delta <- ifelse(data.tmp$event.time <= data.tmp$cens, 1, 0)
        data.tmp$time <- rep(999, 4)
        data.tmp$time[data.tmp$proband == 1] <- data.tmp$event.time[data.tmp$proband == 1]
        data.tmp <- data.tmp[order(data.tmp$proband, decreasing=TRUE), ]
        data.tmp$ind.id <- seq(0, 3, by=1)
        outdata <- rbind(outdata, data.tmp)
        ith <- ith + 1
        break
      }else{
        nrej <- nrej + 1
      }
    }
  }
  return(outdata)
}

# ----------------------------------------------------------
# Marginal function of Weibull distribution
# ----------------------------------------------------------

marg.f <- function(lam, kap, beta, xx, tt){
  p <- length(beta)
  xx <- matrix(xx, ncol=p)
  beta <- matrix(beta, ncol=1)
  tmp <- as.numeric(xx %*% beta)
  St <- exp(-((lam*tt)^kap)*exp(tmp))
  ft <- lam*kap*((lam*tt)^(kap-1))*exp(tmp)*exp(-((lam*tt)^kap)*exp(tmp))
  Ft <- 1-exp(-((lam*tt)^kap)*exp(tmp))
  out <- data.frame(St = St, Ft = Ft, ft = ft)
  return(out)
}

# ----------------------------------------------------------
# F(t1, t2) via Gaussian copula
# ----------------------------------------------------------

biv.f <- function(q1, q2, rho){
  value <- pmnorm(x=c(q1, q2), mean=rep(0, 2), varcov = p2P(param=rho, d=2))[1]
  return(value)
}

# ----------------------------------------------------------
# F(t1, t2, t3) via Gaussian copula
# ----------------------------------------------------------

triv.f <- function(q1, q2, q3, rho12, rho13, rho23){
  covmat <- p2P(param=c(rho12, rho13, rho23), d=3)
  value <- pmnorm(x=c(q1, q2, q3), mean=rep(0, 3), varcov=covmat)[1]
  return(value)
}

# ----------------------------------------------------------
# F(t1, t2, t3, t4) via Gaussian copula
# ----------------------------------------------------------

fourv.f <- function(q1, q2, q3, q4, rho12, rho13, rho14,
                    rho23, rho24, rho34){
  covmat <- p2P(param=c(rho12, rho13, rho14, rho23, rho24, rho34), d=4)
  value <- pmnorm(x=c(q1, q2, q3, q4), mean=rep(0, 4), varcov=covmat)[1]
  return(value)}

# ----------------------------------------------------------
# Construct relationship code for two datasets.
# ----------------------------------------------------------

bi.rel.func <- function(rho.vec, data1, data2){
  rho.pp <- rho.vec[1]
  rho.ss <- rho.vec[2]
  rho.ps <- rho.vec[3]
  rel.id1 <- data1$rel.id
  rel.id2 <- data2$rel.id
  rho <- rep(rho.pp, dim(data1)[1])
  rho <- ifelse(rel.id1 == rel.id2 & rel.id1 == 2, rho.ss, rho)
  rho <- ifelse(rel.id1 != rel.id2, rho.ps, rho)
  return(rho)
}

# ----------------------------------------------------------
# Conditional Expectation of E[Yi|t0] and  E[Zi|t0].
# ----------------------------------------------------------
cexp.givt0.f <- function(rho.vec, Ft.vec, data){
  q.vec <- qnorm(Ft.vec, mean=0, sd=1)
  msize <- max(data$ind.id)
  rho0j.list <- sapply(1:msize, function(r){bi.rel.func(rho.vec, data1=data[data$ind.id == 0, ], data2=data[data$ind.id == r, ])})
  mu.vec <- sapply(1:msize, function(k){
    pnorm(q=(q.vec[k+1] - q.vec[1]*rho0j.list[k])/sqrt(1-rho0j.list[k]^2), mean=0, sd=1)})
  mat2 <- data.frame(j = rep(seq(1, msize), each=msize), k=rep(seq(1, msize), msize))
  mat2 <- mat2[mat2$j < mat2$k, ]
  eta.vec <- sapply(1:dim(mat2)[1], function(r){
    j <- mat2[r, ]$j
    k <- mat2[r, ]$k
    rhojk <- bi.rel.func(rho.vec, data[data$ind.id==j, ], data[data$ind.id == k, ])
    rho0j <- rho0j.list[j]
    rho0k <- rho0j.list[k]
    tmp1 <- c(q.vec[j+1], q.vec[k+1]) - c(rho0j, rho0k)*q.vec[1]
    tmp2 <- matrix(c(1-rho0j^2, rhojk - rho0j*rho0k, rhojk - rho0j*rho0k, 1-rho0k^2), ncol=2)
    value <- pmnorm(x=tmp1, mean=rep(0, 2), varcov=tmp2)[1]
  })
  mean.vec <- c(mu.vec, eta.vec)
  return(mean.vec)
}

# ----------------------------------------------------------
# Full covariance matrix for (Yi, Zi | Ti0)
# ----------------------------------------------------------
cov.givt0.f <- function(rho.vec, Ft.vec, mean.vec, mat2, data){
  q.vec <- qnorm(Ft.vec, mean=0, sd=1)
  msize <- max(data$ind.id)
  rho0j.list <- sapply(1:msize, function(r){bi.rel.func(rho.vec, data1=data[data$ind.id == 0, ], data2=data[data$ind.id == r, ])})
  mu.vec <- mean.vec[1:msize]
  eta.vec <- mean.vec[-seq(1, msize)]
  cov11 <- p2P(param=eta.vec, d=msize) - diag(1, nrow=msize, ncol=msize) + diag(mu.vec, nrow=msize, ncol=msize) - matrix(mu.vec, ncol=1)%*%matrix(mu.vec, nrow=1)
  mat2.value <- mat2
  mat2.value$value <- eta.vec

  mat3 <- data.frame(id=seq(1, msize*dim(mat2)[1]), s=rep(seq(1, msize), each=dim(mat2)[1]), j=rep(mat2[,1], msize), k=rep(mat2[,2], msize))

  value3.tmp <- sapply(1:dim(mat3)[1], function(r){
    iid <-  sort(unique(as.numeric(mat3[r, c("s", "j", "k")])))
    if(length(iid) == 2){
      value <- mat2.value[mat2.value$j == iid[1] & mat2.value$k == iid[2], ]$value
    }
    if(length(iid) == 3){
      s <- iid[1]; j <- iid[2]; k <- iid[3]
      rho.sj <- bi.rel.func(rho.vec, data[data$ind.id==s, ], data[data$ind.id==j, ])
      rho.sk <- bi.rel.func(rho.vec, data[data$ind.id==s, ], data[data$ind.id==k, ])
      rho.jk <- bi.rel.func(rho.vec, data[data$ind.id==j, ], data[data$ind.id==k, ])
      rho.0s <- rho0j.list[s]
      rho.0j <- rho0j.list[j]
      rho.0k <- rho0j.list[k]
      tmp1 <- c(q.vec[(s+1)], q.vec[(j+1)], q.vec[(k+1)]) - q.vec[1]*c(rho.0s, rho.0j, rho.0k)
      tmp2 <- matrix(c(1, rho.sj, rho.sk, rho.sj, 1, rho.jk, rho.sk, rho.jk, 1), ncol=3, byrow=TRUE) - matrix(c(rho.0s^2, rho.0s*rho.0j, rho.0s*rho.0k, rho.0j*rho.0s, rho.0j^2, rho.0j*rho.0k, rho.0k*rho.0s, rho.0j*rho.0k, rho.0k^2), ncol=3, byrow=TRUE)
      value <- pmnorm(x=tmp1, mean=rep(0, 3), varcov=tmp2)[1]
    }
    return(value)})
  mat3$value <- as.numeric(value3.tmp)
  cov12 <- matrix(mat3$value, nrow=msize, ncol=dim(mat2)[1], byrow=TRUE) - matrix(mu.vec, ncol=1)%*%matrix(eta.vec, nrow=1)

  mat4 <- data.frame(id=seq(1, (dim(mat2)[1])^2), s=rep(mat2$j, each=dim(mat2)[1]), l=rep(mat2$k, each=dim(mat2)[1]), j=rep(mat2$j, dim(mat2)[1]), k=rep(mat2$k, dim(mat2)[1]))

  value4.tmp <- sapply(1:dim(mat4)[1], function(r){
    iid <- sort(unique(as.numeric(mat4[r, c("s", "l", "j", "k")])))
    if(length(iid) == 2){
      value <- unique(mat2.value$value[mat2.value$j == iid[1] & mat2.value$k == iid[2]])
    }
    if(length(iid) == 3){
      mat3.tmp <- mat3[, c("s", "j", "k")]
      mat3.tmp <- t(apply(mat3.tmp, 1, sort))
      mat3.tmp <- cbind(mat3.tmp, mat3$value)
      value <- unique(mat3.tmp[mat3.tmp[,1] == iid[1] & mat3.tmp[, 2] == iid[2] & mat3.tmp[,3] == iid[3], 4])
    }
    if(length(iid) == 4){
      s <- iid[1]; l <- iid[2]; j <- iid[3]; k <- iid[4]
      rho.sl <- bi.rel.func(rho.vec, data[data$ind.id==s, ], data[data$ind.id==l, ])
      rho.sj <- bi.rel.func(rho.vec, data[data$ind.id==s, ], data[data$ind.id==j, ])
      rho.sk <- bi.rel.func(rho.vec, data[data$ind.id==s, ], data[data$ind.id==k, ])
      rho.lj <- bi.rel.func(rho.vec, data[data$ind.id==l, ], data[data$ind.id==j, ])
      rho.lk <- bi.rel.func(rho.vec, data[data$ind.id==l, ], data[data$ind.id==k, ])
      rho.jk <- bi.rel.func(rho.vec, data[data$ind.id==j, ], data[data$ind.id==k, ])
      rho.0s <- rho0j.list[s]; rho.0l <- rho0j.list[l]
      rho.0j <- rho0j.list[j]; rho.0k <- rho0j.list[k]
      tmp1 <- c(q.vec[s+1], q.vec[l+1], q.vec[j+1], q.vec[k+1]) - c(rho.0s, rho.0l, rho.0j, rho.0k)*q.vec[1]
      tmp2 <- p2P(c(rho.sl - rho.0s*rho.0l, rho.sj - rho.0s*rho.0j, rho.sk-rho.0s*rho.0k, rho.lj-rho.0l*rho.0j, rho.lk-rho.0l*rho.0k, rho.jk-rho.0j*rho.0k), d=4)
      tmp2 <- tmp2 - diag(c(rho.0s^2, rho.0l^2, rho.0j^2, rho.0k^2), ncol=4, nrow=4)
      value <- pmnorm(x=tmp1, mean=rep(0, 4), varcov=tmp2)[1]
    }
    return(value)})
  mat4$value <- as.numeric(value4.tmp)
  cov22 <- matrix(mat4$value, nrow=dim(mat2)[1], ncol=dim(mat2)[1], byrow=TRUE) - matrix(eta.vec, ncol=1)%*%matrix(eta.vec, nrow=1)
  cov.mat <- rbind(cbind(cov11, cov12), cbind(t(cov12), cov22))
  return(cov.mat)
}

# ----------------------------------------------------------
# Working partial independence (WPI) matrix : W_12 = W_21 = 0
# ----------------------------------------------------------
cov.givt0.simple2.f <- function(rho.vec, Ft.vec, mean.vec, mat2, data){
  q.vec <- qnorm(Ft.vec, mean=0, sd=1)
  msize <- max(data$ind.id)
  rho0j.list <- sapply(1:msize, function(r){bi.rel.func(rho.vec, data1=data[data$ind.id == 0, ], data2=data[data$ind.id == r, ])})
  mu.vec <- mean.vec[1:msize]
  eta.vec <- mean.vec[-seq(1, msize)]
  cov11 <- p2P(param=eta.vec, d=msize) - diag(1, nrow=msize, ncol=msize) + diag(mu.vec, nrow=msize, ncol=msize) - matrix(mu.vec, ncol=1)%*%matrix(mu.vec, nrow=1)

  cov12 <- matrix(0, nrow=msize, ncol=length(eta.vec))
  cov22 <- diag(eta.vec*(1-eta.vec), nrow=length(eta.vec), ncol=length(eta.vec))
  cov.mat <- rbind(cbind(cov11, cov12), cbind(t(cov12), cov22))
  return(cov.mat)
}

# ----------------------------------------------------------
# Numerical gradient of mean vector
# based on Structure Gaussian
# Here para=(loglam, logkap, beta, gam0, gam1, gam2)
# ----------------------------------------------------------
num.grad.str.f <- function(para, data){
  mean.f <- function(par, data){
    lam <- exp(par[1])
    kap <- exp(par[2])
    beta <- par[3]
    gam.vec <- par[4:6]
    alp.vec <- gam.vec + c(0, rep(gam.vec[1], 2))
    tau.vec <- (exp(alp.vec) - 1)/(exp(alp.vec) + 1)
    rho.vec <- sin(pi*tau.vec/2)
    Ft.vec <- marg.f(lam=lam, kap=kap, beta=beta, xx=data$x1, tt=c(data[data$ind.id==0,]$time, data[data$ind.id!=0, ]$cens))$Ft
    mean.vec <- cexp.givt0.f(rho.vec, Ft.vec, data)
    return(mean.vec)}
  eps <- 	1e-05
  out <- sapply(1:length(para), function(r){
    eps.vec <- rep(0, length(para))
    eps.vec[r] <- eps

    para.l <- para - eps.vec
    para.r <- para + eps.vec
    value <- (mean.f(par=para.r, data=data) - mean.f(par=para.l, data=data))/(2*eps)
    return(value)})
  return(matrix(out, ncol=length(para)))
}

# ----------------------------------------------------------
# Conditional GEE2 for nonproband|T0
# Based on Structure Gaussian
# ----------------------------------------------------------
cgee2.ind.str.f <- function(para, data, der.method=c("full", "simple1", "simple2"), cov.method=c("full",  "simple2")){
  der.method <- match.arg(der.method)
  cov.method <- match.arg(cov.method)
  lam <- exp(para[1])
  kap <- exp(para[2])
  beta <- para[3]
  gam.vec <- para[4:6]
  alp.vec <- gam.vec + c(0, rep(gam.vec[1], 2))
  tau.vec <- (exp(alp.vec) - 1)/(exp(alp.vec) + 1)
  rho.vec <- sin(pi*tau.vec/2)
  y.vec <- data[data$ind.id!=0, ]$delta
  msize <- max(data$ind.id)
  mat2 <- data.frame(j = rep(seq(1, msize), each=msize), k=rep(seq(1, msize), msize), yj=rep(y.vec, each=msize), yk=rep(y.vec, msize))
  mat2 <- mat2[mat2$j < mat2$k, ]

  Ft.vec <- marg.f(lam=lam, kap=kap, beta=beta, xx=data$x1, tt=c(data[data$ind.id==0,]$time, data[data$ind.id!=0, ]$cens))$Ft
  mean.vec <- cexp.givt0.f(rho.vec, Ft.vec, data)
  # Full derivative matrix
  der.mat1 <- der.mat2 <-  der.mat3 <- num.grad.str.f(para=para, data=data)
  # der.mat2: G12=G21=0 and der.mat3: G12 = 0
  der.mat2[1:msize, (4:6)] <- matrix(0, nrow=msize, ncol=3)
  der.mat2[(msize + 1):length(mean.vec), 1:3] <- matrix(0, nrow=(length(mean.vec)-msize), ncol=3)
  der.mat3[(msize + 1):length(mean.vec), 1:3] <- matrix(0, nrow=(length(mean.vec)-msize), ncol=3)

  if(der.method == "full"){
    der.mat <- der.mat1
  }
  if(der.method == "simple1"){
    der.mat <- der.mat3
  }
  if(der.method == "simple2"){
    der.mat <- der.mat2
  }
  if(cov.method == "full"){
    cov.mat <- cov.givt0.f(rho.vec, Ft.vec, mean.vec, mat2, data)
  }
  if(cov.method == "simple2"){
    cov.mat <- cov.givt0.simple2.f(rho.vec, Ft.vec, mean.vec, mat2, data)
  }
  z.vec <- mat2$yj*mat2$yk
  icov.mat <- ginv(cov.mat)
  out <- t(der.mat)%*%icov.mat%*%matrix(c(y.vec, z.vec) - mean.vec, ncol=1)
  output <- list(NULL)
  output[[1]] <- out
  output[[2]] <- der.mat
  output[[3]] <- icov.mat
  output[[4]] <- mean.vec
  output[[5]] <- der.mat1
  return(output)
}

# ----------------------------------------------------------
# E[T0|Y0=1]
# ----------------------------------------------------------
mu0.f <- function(para, xx, cens0){
  lam <- exp(para[1])
  kap <- exp(para[2])
  beta <- para[3]
  integrand1 <- function(tt){
    tt*(lam*kap)*((lam*tt)^(kap-1))*exp(beta*xx) * exp( (-1)*((lam*tt)^kap)*exp(beta*xx) )
  }
  Fcens <- marg.f(lam=lam, kap=kap, beta=beta, xx=xx, tt=cens0)$Ft
  value <- integrate(function(tt){sapply(tt, function(tt){integrand1(tt)})}, lower=0, upper=cens0)$value/Fcens
  return(value)
}

# ----------------------------------------------------------
# E[T0^2|Y0 = 1]
# ----------------------------------------------------------
eta0.f <- function(para, xx, cens0){
  lam <- exp(para[1])
  kap <- exp(para[2])
  beta <- para[3]
  Fcens <- marg.f(lam=lam, kap=kap, beta=beta, xx=xx, tt=cens0)$Ft
  integrand2 <- function(tt){
    (tt^2)*(lam*kap)*((lam*tt)^(kap-1))*exp(beta*xx) * exp( (-1)*((lam*tt)^kap)*exp(beta*xx) )
  }
  value <- integrate(function(tt){sapply(tt, function(tt){integrand2(tt)})}, lower=0, upper=cens0)$value/Fcens
  return(value)
}

# ----------------------------------------------------------
# CGEE for T0|Y0
# ----------------------------------------------------------
cgee.t0.ind.f <- function(para, data){
  p <- length(para)
  lam <- exp(para[1])
  kap <- exp(para[2])
  beta <- para[3]
  t0 <- data[data$ind.id==0,]$time
  cens0 <- data[data$ind.id == 0, ]$cens
  xx <- data[data$ind.id == 0, ]$x1
  mu0 <- mu0.f(para=para, xx=xx, cens0=cens0)
  eta0 <- eta0.f(para=para, xx=xx, cens0=cens0)
  var0 <- eta0 - mu0^2
  eps <- 1e-05
  der0 <- sapply(1:3, function(r){
    eps.vec <- rep(0, p)
    eps.vec[r] <- eps
    term <- (mu0.f(para=para+eps.vec, xx=xx, cens0=cens0) - mu0.f(para=para-eps.vec, xx=xx, cens0=cens0))/(2*eps)
    return(term)})
  der0 <- c(der0, rep(0,(p-3)))
  out <- as.numeric(der0*(1/var0)*(t0 - mu0))
  output <- list(NULL)
  output[[1]] <- out
  output[[2]] <- der0
  output[[3]] <- 1/var0
  output[[4]] <- mu0
  return(output)
}

# ----------------------------------------------------------
# Newton Raphson for CGEE2 + CGEE for T0|Y0=1.
# With Structure Gaussian
# ----------------------------------------------------------
NR.adt0.str.f <- function(para, alldata, der.method=c("full", "simple1", "simple2"), cov.method=c("full", "simple2")){
  der.method <- match.arg(der.method)
  cov.method <- match.arg(cov.method)
  p <- length(para)
  fam.id <- sort(unique(alldata$fam.id))
  tmp <- sapply(1:length(fam.id), function(i){
    inner1 <- cgee2.ind.str.f(para=para, data=alldata[alldata$fam.id == fam.id[i], ], der.method, cov.method)
    inner2 <- cgee.t0.ind.f(para=para, data=alldata[alldata$fam.id == fam.id[i], ])
    der.mat1 <- inner1[[2]]
    icov.mat1 <- inner1[[3]]
    eq.value1 <- inner1[[1]]
    der.matf <- inner1[[5]]

    der.mat0 <- matrix(inner2[[2]], ncol=p)
    icov.mat0 <- inner2[[3]]
    eq.value0 <- inner2[[1]]

    eq.value <- as.numeric(eq.value1) + as.numeric(eq.value0)

    tmp2 <- t(der.mat1)%*%icov.mat1%*%der.matf + icov.mat0*t(der.mat0)%*%der.mat0

    result <- c(as.numeric(eq.value), as.numeric(tmp2))
    return(result)})
  out <- apply(tmp, 1, sum)
  term1 <- matrix(out[1:p], ncol=1)
  term2 <- matrix(out[(p+1): (p^2 + p)], ncol=p)
  update <- para + as.numeric(ginv(term2)%*%term1)
  output <- list(NULL)
  output[[1]] <- update
  output[[2]] <- tmp[1:p, ]
  output[[3]] <- term2
  return(output)
}


# ----------------------------------------------------------
# Newton raphson for CGEE2 + GEE for T0|Y0=1.
# With Structure Gaussian
# ----------------------------------------------------------
script.addprob.str.f <- function(instats, indata, der.method=c("full", "simple1", "simple2"), cov.method=c("full", "simple2")){
  der.method <- match.arg(der.method)
  cov.method <- match.arg(cov.method)
  diff.new <- diff.old <- m.score.new <- 1000
  ith <- 0
  break.id <- 0
  para.new <- c(log(instats$lam), log(instats$kap), instats$beta, as.numeric(instats$gam.vec))
  while(diff.new > 1e-06|m.score.new > 1e-05){
    diff.old <- diff.new
    para.old <- para.new
    m.score.old <- m.score.new
    out <- NR.adt0.str.f(para.old, indata, der.method, cov.method)
    para.new <- out[[1]]
    diff.new <- max(abs(para.new - para.old))
    ith <- ith + 1
    m.score.new <- max(abs(apply(out[[2]],1,sum)))
    if((m.score.new-m.score.old) > 7 |diff.new - diff.old > 5|ith >= 50){
      break.id <- 1; break
    }
    cat('ith is: ', ith, '\n')
    cat('diff is: ', diff.new, '\n')
    cat('score is: ', m.score.new, '\n')
    cat('new para: ', para.new, '\n')
  }
  if(break.id == 1){
    return(c(rep(999, 27)))
  }else{
    iAmat <- ginv(out[[3]])
    score.mat <- out[[2]]
    Bmat <- sapply(1:dim(score.mat)[2], function(j){as.numeric(matrix(score.mat[, j], ncol=1)%*%matrix(score.mat[, j], nrow=1))})
    Bmat <- matrix(apply(Bmat, 1, sum), ncol=length(para.new))
    R.var <- iAmat%*%Bmat%*%t(iAmat)
    alp.est <- para.new[4:6] + c(0, rep(para.new[4], 2))
    tau.est <- (exp(alp.est)-1)/(exp(alp.est)+1)
    Hmat <- matrix(c(1, 0, 0, 1, 1, 0, 1, 0, 1), ncol=3, byrow=TRUE)
    alp.var <- Hmat%*%R.var[4:6, 4:6]%*%t(Hmat)
    tau.var <- (((1+tau.est)*(1-tau.est)/2)^2)*diag(alp.var)
    return(c(as.numeric(para.new), c(alp.est, tau.est), diag(R.var), diag(alp.var), tau.var, max(abs(apply(score.mat, 1, sum))), diff.new, ith))
  }
}
