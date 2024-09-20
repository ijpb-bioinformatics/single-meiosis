source('/home/dcharif/save/bioinfo-dev/r-dev/hmm-nco/functions/FunctionsMLE.R')
library(SDMTools)

EstepHMM <- function(Mstep){
  FB = ForBack(Mstep$Phi, Mstep$nu, Mstep$Pi, nrow(Mstep$Phi), ncol(Mstep$Phi))
  return(list(Tau = FB$Tau, logL = FB$logL, Eta = FB$Eta))
}
 
CalPhiBetaBinom <- function(Y, R, gamma ,log=F){
  Phi = matrix(0, length(Y), nrow(gamma))
  for (k in 1:nrow(gamma)){Phi[, k] = dbetabinom(Y, R, gamma[k, 1], gamma[k, 2], log=log)}
  # Modif 22/9/16 : plantage si Phi = 0
  if (log==T){
     Phi[which(Phi < -300)] = -300
  }else{
     Phi[which(Phi==0)] = exp(-300)
  }
  return(Phi)
}

MstepBetaBinom <- function(Y, R, Tau, gamma.init=matrix(1, ncol(Tau), 2), traceVGAM=F, tolSd = 1e-4, parmBetaMax = 1e3){
   K = nrow(gamma.init); gamma = gamma.init; X = Y/R
   for (k in 1:K){
      cat(k, ':')
      mXk = wt.mean(X, Tau[, k])
      cat(mXk, wt.sd(X, Tau[, k]), '')
      if (mXk*(1-mXk)  < tolSd){
         if (mXk < .5){
            gamma[k, ] = c(1/parmBetaMax, parmBetaMax)
         }else{
            gamma[k, ] = c(parmBetaMax, 1/parmBetaMax)
         }
      }else{
         # gamma[k, ] = BetaBinomial.VGAM.MLE(cbind(Y, R-Y), Tau[, k], alphabeta.init=gamma[k, ], traceVGAM=traceVGAM)
         # Modif 23/9/16 : pb with null Tau's
         Selk = which(Tau[, k]>0); Yk = Y[Selk]; Rk = R[Selk]; Tauk = Tau[Selk, k]
         gamma[k, ] = BetaBinomial.VGAM.MLE(cbind(Yk, Rk-Yk), Tauk, alphabeta.init=gamma[k, ], traceVGAM=traceVGAM)
         gamma[k, which(gamma[k, ] > parmBetaMax)] = parmBetaMax
         gamma[k, which(gamma[k, ] < 1/parmBetaMax)] = parmBetaMax
      }
      # gamma[k, ] = BetaBinomial.MLE(cbind(Y, R-Y), Tau[, k], alpha.init=gamma.init[k, ])
      # gamma[k, ] = BetaBinomial.MLOOE(cbind(Y, R-Y), R, Tau[, k])
      cat(gamma[k, ], '/ ')
   }
   cat('\n')
   Phi = CalPhiBetaBinom(Y, R, gamma)
   logPhi = CalPhiBetaBinom(Y, R, gamma, log=T)
   logPhi[which(logPhi < -300)] = -300
   return(list(gamma=gamma, Phi=Phi, logPhi=logPhi))
}   

MstepHMMBetaBinom <- function(Y, R, Estep, gamma.init=matrix(1, ncol(Estep$Tau), 2), traceVGAM=F){
   cat('E-step: ')
   Pi = Estep$Eta; Pi = Pi / rowSums(Pi)
   nu = StatDist(Pi); K = length(nu)
   MstepBBParms = MstepBetaBinom(Y, R, Estep$Tau, gamma.init, traceVGAM=F)
   gamma = MstepBBParms$gamma; Phi = MstepBBParms$Phi; logPhi = MstepBBParms$logPhi
   alpha.tau = gamma[, 1] / rowSums(gamma)
   order = order(alpha.tau)
   gamma = gamma[order, ]; Pi = Pi[order, order];  nu = nu[order]
   return(list(Pi=Pi, nu=nu, gamma=gamma, Phi=Phi))
}

FitHMMBetaBinom <- function(Y, R, K, tolCvg=1e-5, epsTau=1e-4, iterMax=50, traceVGAM=F){
  # R = nb trials, Y = nb success
  n = length(Y)
  # Init
  library(mclust); Estep = list()
  GMM = Mclust(asin(sqrt(Y/R)), G=K, modelNames='E')
  Estep$Tau = GMM$z
  Estep$Eta = t(Estep$Tau[1:(n-1), ]) %*% Estep$Tau[2:n, ]
  Mstep = MstepHMMBetaBinom(Y, R, Estep, traceVGAM=traceVGAM)
  Theta = c(as.vector(Mstep$Pi), as.vector(t(Mstep$gamma)))
  
  # E-M
  diff = 2*tolCvg; iter = 0
  while((diff > tolCvg) & (iter < iterMax)){
     iter = iter + 1
     # Plot beta distributions
     p = seq(0, 1, by=1e-3); par(mfrow=c(K, 1))
     plot(p, dbeta(p, Mstep$gamma[1, 1], Mstep$gamma[1, 2]), type='l', lwd=2, col=2, ylab='', xlab='')
     lines(p, dbeta(p, Mstep$gamma[2, 1], Mstep$gamma[2, 2]), lwd=2, col=4)
   
     # Reguler EM  
     Estep = EstepHMM(Mstep)
     Mstep = MstepHMMBetaBinom(Y, R, Estep, gamma.init=Mstep$gamma, traceVGAM=traceVGAM)
     Theta.new = c(as.vector(Mstep$Pi), as.vector(t(Mstep$gamma)))
     
     # # Accelated EM
     # Estep = EstepHMM(Mstep)
     # Mstep = MstepHMMBetaBinom(Y, R, Estep, gamma.init=Mstep$gamma)
     # # Theta1 = c(as.vector(Mstep$Pi), as.vector(t(Mstep$gamma))); # cat(Theta1)
     # Estep = EstepHMM(Mstep)
     # Mstep = MstepHMMBetaBinom(Y, R, Estep)
     # Theta2 = c(as.vector(Mstep$Pi), as.vector(t(Mstep$gamma))); # cat(Theta2)
     # alpha = + sqrt(sum(Theta1 - Theta)^2) / sqrt(sum(Theta2 - 2*Theta1 + Theta)^2)
     # Theta.new = Theta + alpha*(Theta1 - Theta)
     # while(min(Theta.new)<0){alpha = alpha/2; Theta.new = Theta + alpha*(Theta1 - Theta)}
     # Mstep$Pi = matrix(Theta.new[1:K^2], K, K); Mstep$nu = StatDist(Mstep$Pi)
     # Mstep$gamma = matrix(Theta.new[(K^2+1):((K^2)+(2*K))], K, 2, byrow=T); 
     # Mstep$Phi = CalPhiBetaBinom(Y, R, Mstep$gamma)
     
     # Test & mise a jour
     diff = max(abs(Theta.new - Theta))
     Theta = Theta.new; # cat(Theta)
     
     cat(iter, as.vector(Mstep$Pi), as.vector(t(Mstep$gamma)), Estep$logL, diff, '\n')
     }

  HMM = list(tau=Estep$Tau, Pi=Mstep$Pi, nu=Mstep$nu, gamma=Mstep$gamma, Phi=Mstep$Phi, logL=Estep$logL)
  return(HMM)
}
  
