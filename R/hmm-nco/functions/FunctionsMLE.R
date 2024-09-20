# Various functions for MLE fitting
library(VGAM)

################################################################
# Gamma distribution
################################################################
GammaMLE.moments <- function(mean.X, mean.logX, tol = 1e-4){
   # fits MLE estimates for a Gamma distribution based on the sufficient statistics
   target = log(mean.X) - mean.logX
   ainf = 1; while(log(ainf)-digamma(ainf)<target){ainf = ainf/2}
   asup = 1; while(log(asup)-digamma(asup)>target){asup = asup*2}
   a = (asup+ainf)/2
   while(asup - ainf > tol){
#       cat(ainf, a, asup, asup - ainf, '\n')
      if (log(a)-digamma(a)>target){ainf=a}else{asup=a}   
      a = (asup+ainf)/2
   }
   shape = a
   rate = shape/mean.X
   invisible(list(shape=shape, rate=rate))
}

GammaMLE <- function(X){
   # fits MLE estimates for a Gamma distribution based on the original sample
   MLE = GammaMLE.moments(mean(X), mean(log(X)))
   shape = MLE$shape
   rate = MLE$rate
   invisible(list(shape=shape, rate=rate))
}

################################################################
# Dirichlet distribution
################################################################
DirichletMLE.moments <- function(mean.logX, n, a.init, tol=1e-6, max.iter = 100){
    # fits the MLE estimates for a Dirichlet distribution using Minka (00) fix-point algo
   
   diff = 2*tol; a = a.init; K = length(a.init); iter = 0
   while((diff > tol) & (iter < max.iter)){
      iter = iter+1
#       cat(a, diff, '\n')
      psi = digamma(sum(a)) + mean.logX
      a.new = rep(0, K)
      invisible(lapply(1:K, function(k){a.new[k] <<- inv.digamma(psi[k])}))
      diff = max(abs(a.new-a))
      a = a.new
   }
   InvVar = matrix(trigamma(sum(a)), K, K)
   diag(InvVar) = diag(InvVar) - trigamma(a)
   Var = -solve(InvVar)/n
   return(list(a = a, Var = Var))
}
   
inv.digamma <- function(y, tol=1e-6){
   # invert the digamma functiuon using dichotimy
   
   xsup = 1; while(digamma(xsup) < y){xsup = 2*xsup}
   xinf = 1; while(digamma(xinf) > y){xinf = xinf/2}
   x = sqrt(xinf*xsup)
   while((xsup - xinf)/x > tol){
#       cat(xinf, x, xsup, digamma(x), y, '\n')
      if (digamma(x) > y){xsup = x}else{xinf=x}
      x = sqrt(xinf*xsup)
   }
   invisible(x)
}

################################################################
# Beta-binomial distribution
################################################################
# Beta-binomial distribution BB(R, alpha, beta):
# P_i ~ Beta(alpha, beta); Y_i|P_i ~ B(N_i, P_i)
# R = vector of number of trials, Y = vector of number of successes, W = vector of weights
dbetabinom <- function(Y, R, alpha, beta, log=F){
   prob = (lgamma(R+1) - lgamma(Y+1) - lgamma(R-Y+1)) +
      (lgamma(Y+alpha) + lgamma(R-Y+beta) - lgamma(R+alpha+beta)) +
      (lgamma(alpha+beta) - lgamma(alpha) - lgamma(beta)) 
   if (log==F){prob = exp(prob)}
   invisible(prob)
}

negloglik.betabinom <- function(alphabeta, Y, R, W=rep(1, length(Y))){
   return(-sum(W * dbetabinom(Y, R, alphabeta[1], alphabeta[2], log=T)))
}

# BetaBinomial.MLE <- function(Ymat, W=rep(1, length(Y)), tolCvg=1e-3, tolAlpha = 1e-6){
#    R = rowSums(Ymat); n = nrow(Ymat); 
#    Sel = which(W > 0); Ymat = Ymat[Sel, ]; W = W[Sel]
#    Wcum = sum(W)
#    Pmat = (Ymat +1) / (R + 1)
#    alpha = DirichletMLE.moments(mean.logX=t(W)%*%log(Pmat), n=Wcum, a.init=colMeans(Pmat))$a
#    diff = 2*tolCvg; iter = 0
#    # cat('a', n, dim(Ymat), length(alpha), length(W), '\n')
#    while (diff > tolCvg){
#       iter = iter +1
#       Numer.vec = t(W) %*% digamma(Ymat + outer(rep(1, n), alpha)) - Wcum*digamma(alpha)
#       Denom = t(W) %*% digamma(R + sum(alpha)) - Wcum*digamma(sum(alpha))
#       alpha1 = as.vector((alpha * Numer.vec) / Denom[1, 1])
#       alpha1[which(alpha1<tolAlpha)] = tolAlpha
#       
#       # Numer.vec = t(W) %*% digamma(Ymat + outer(rep(1, n), alpha1)) - digamma(alpha1)
#       # Denom = t(W) %*% digamma(R + sum(alpha1)) - digamma(sum(alpha1))
#       # alpha2 = as.vector((alpha1 * Numer.vec) / Denom[1, 1])
#       # alpha2[which(alpha2<tolAlpha)] = tolAlpha
#       # 
#       # Coef =   sqrt(sum(alpha1 - alpha)^2) / sqrt(sum(alpha2 - 2*alpha1 + alpha)^2)
#       # alpha.new = alpha + Coef*(alpha1-alpha)
#       # while(min(alpha.new)<tolAlpha){
#       #    Coef = Coef/2
#       #    alpha.new = alpha + Coef*(alpha1-alpha)
#       # }
#       alpha.new = alpha1
# 
#       diff = max(abs((alpha.new-alpha)/alpha))
#       alpha = alpha.new
#       # if(iter%%10==0){cat(iter, alpha, diff, '\n')}
#    }
#    # cat(iter, alpha, diff, '\n')
#    return(alpha)
# }

# BetaBinomial.MLOOE <- function(Ymat, W=rep(1, length(R)), tolCvg=1e-6){
#    # Maximum leave-one-out likelihood estimates (Fix-point relation, Minka, 00)
#    R = rowSums(Ymat); n = nrow(Ymat); W = W / sum(W)
#    alpha = colMeans(Ymat)
#    diff = 2*tolCvg
#    # cat('a', n, dim(Ymat), length(alpha), length(W), '\n')
#    while (diff > tolCvg){
#       Numer.vec = Ymat + outer(rep(1, n), alpha) - 1
#       # cat('b', n, dim(Ymat), length(alpha), length(W), dim(Numer.vec), '\n')
#       Numer.vec = as.vector(t(W) %*% (Ymat / Numer.vec))
#       # cat('c', n, dim(Ymat), length(alpha), length(W), dim(Numer.vec), '\n')
#       Denom = t(W) %*% (R / (R - 1 + sum(alpha)))
#       # cat('d', n, dim(Ymat), length(alpha), length(W), dim(Numer.vec))
#       alpha.new = as.vector((alpha * Numer.vec) / Denom)
#       diff = max(abs((alpha.new-alpha)))
#       alpha = alpha.new
#       cat(alpha, diff, '\n')
#    }
#    return(alpha)
# }

BetaBinomial.VGAM.MLE <- function(Ymat, W=rep(1, nrow(Ymat)), 
                                  alphabeta.init=c(1, 1), traceVGAM=F){
   murho.init = c(alphabeta.init[1]/sum(alphabeta.init), 1/(1+sum(alphabeta.init)))
   VGLM = vglm(Ymat ~ 1, betabinomial, weights=W, 
               coefstart=qlogis(murho.init), trace=traceVGAM)
   # print(VGLM); summary(VGLM)
   mu = Coef(VGLM)[1]; rho = Coef(VGLM)[2]; 
   return((1-rho)/rho*c(mu, 1-mu))
}
   