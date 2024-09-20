########################################################
F_EMBinContHMMBinom1D <- function(Y, R, diffTime, RateInit, tolCvg=1e-5, iterMax=100){
  # Fits a binary continous time HMM to beta-binomial data via EM
  # Y[i] = nb success at time t_i (n*1), 
  # R[i] = nb trials at time t_i (n*1)
  # diffTime[i] = t_(i+1) - t_i
  # RateInit = initial guess for the transition rates (1->2 and 2->1)
  # Rate = transition rates, gamma = parms of the (latent) beta distribution in each state
  # Init
  Tmax = sum(diffTime)
  diffTime = diffTime/Tmax
  Tau = Mclust(asin(sqrt(Y/R)), G=2)$z; Tau1 = Tau[1, ]
  Eta = t(sapply(2:n, function(t){as.vector(outer(Tau[t-1,], Tau[t,]))}))
  
  diff = 2*tolCvg; iter = 0; Rate = RateInit*Tmax
  while((diff > tolCvg) & (iter < iterMax)){
    iter = iter+1
    
    # M-step
    cat('Probas ')
    Mstep = MstepBinom1D(Y, R, Tau)
    cat('Rates \n')
    Rate = MstepRateBinContHMM(Rate, diffTime, Tau1, Eta)$Rate
    Order = order(Mstep$gamma)
    Mstep$gamma = Mstep$gamma[Order]
    Rate = Rate[Order]
    Mstep$logPhi = Mstep$logPhi[, Order]
    
    # E step
    PiNu = F_ComputeBinaryPiNu(Rate, diffTime)
    cat('Forward ')
    FW = F_HetForward(PiNu$Pi, PiNu$Nu, Mstep$logPhi)
    cat('Backward ')
    BW = F_HetBackward(PiNu$Pi, FW$Fw)
    
    # Test & update
    cat('Test \n')
    diff = max(abs(Tau - BW$Tau))
    Tau = BW$Tau; Eta = BW$Eta
    cat(iter, Mstep$gamma, Rate/Tmax, FW$logL, diff, '\n')
  }
  Rate = Rate/Tmax
  
  return(list(gamma=Mstep$gamma, Rate=Rate, Tau=Tau, Eta=Eta, logL=FW$logL, 
              logPhi=Mstep$logPhi, Pi = PiNu$Pi, Nu = PiNu$Nu))
}

########################################################
F_EMJointBinContHMMBinom1D_DiffRates <- function(Y, R, diffTime, RateInit, TauInit, tolCvg=1e-5, iterMax=100){
   # Fits a common binary continous time HMM to binomial data via EM on I profiles
   # Model woth same proba (gamma) for all profiles
   # but profil-specific rates
   # Y[u, i] = nb success at time t_u in profile u (n*I), 
   # R[u, i] = nb trials at time t_u in profile i (n*I)
   # diffTime[u] = t_(u+1) - t_u
   # RateInit = initial guess for the transition rates (1->2 and 2->1)
   # Rate = transition rates 
   # gamma = succes probability in each state  
   # Init
   n = nrow(Y); I = ncol(Y); K = 2
   Tmax = sum(diffTime)
   diffTime = diffTime/Tmax
   Tau = TauInit; Eta = list(); 
   invisible(sapply(1:I, function(i){
      Eta[[i]] <<- t(sapply(2:n, function(t){as.vector(outer(Tau[[i]][t-1,], Tau[[i]][t,]))}))
   }))
   
   diff = 2*tolCvg; iter = 0; 
   Rate = lapply(1:I, function(i){RateInit[[i]]*Tmax})
   while((diff > tolCvg) & (iter < iterMax)){
      iter = iter+1
      
      # M-step
      cat('Probas & Phi ')
      Ytau = Rtau = rep(0, K); 
      invisible(sapply(1:I, function(i){
         Ytau <<- Ytau + t(Tau[[i]])%*%Y[, i]
         Rtau <<- Rtau + t(Tau[[i]])%*%R[, i]
      }))
      gamma = as.vector(Ytau / Rtau)
      logPhi = list()
      invisible(sapply(1:I, function(i){
         logPhi[[i]] <<- matrix(0, n, K)
         sapply(1:K, function(k){logPhi[[i]][, k] <<- dbinom(Y[, i], R[, i], gamma[k], log=T)})
         logPhi[[i]][which(logPhi[[i]] < -300)] <<- -300
      }))
      cat('Rates ')
      Tau1 = lapply(1:I, function(i){Tau[[i]][1, ]})
      Rate = lapply(1:I, function(i){cat(i, ''); MstepRateBinContHMM(Rate[[i]], diffTime, Tau1[[i]], Eta[[i]])$Rate})
      
      # E step
      PiNu = lapply(1:I, function(i){F_ComputeBinaryPiNu(Rate[[i]], diffTime)})
      cat('Forward-Backward ')
      FW = BW = list(); logL = 0
      invisible(sapply(1:I, function(i){
         FW[[i]] <<- F_HetForward(PiNu[[i]]$Pi, PiNu[[i]]$Nu, logPhi[[i]])
         BW[[i]] <<- F_HetBackward(PiNu[[i]]$Pi, FW[[i]]$Fw)
         logL <<- logL + FW[[i]]$logL
      }))
      
      # Test & update
      cat('Test \n')
      diff = max(unlist(lapply(1:I, function(i){max(abs(Tau[[i]] - BW[[i]]$Tau))})))
      invisible(sapply(1:I, function(i){Tau[[i]] <<- BW[[i]]$Tau; Eta[[i]] <<- BW[[i]]$Eta}))
      cat(iter, gamma, unlist(lapply(1:I, function(i){Rate[[i]]/Tmax})), logL, diff, '\n')
   }
   Rate = lapply(1:I, function(i){Rate[[i]]/Tmax})
   Pi = Nu = list()
   lapply(1:I, function(i){Pi[[i]] <<- PiNu[[i]]$Pi; Nu[[i]] <<- PiNu[[i]]$Nu})
   
   return(list(gamma=gamma, Rate=Rate, Tau=Tau, Eta=Eta, logL=logL, 
               logPhi=logPhi, Pi = Pi, Nu = Nu))
}

########################################################
F_EMJointBinContHMMBinom1D_SameRates <- function(Y, R, diffTime, RateInit, TauInit, tolCvg=1e-5, iterMax=100){
   # Fits a common binary continous time HMM to binomial data via EM on I profiles
   # Model woth same proba (gamma) and rates
   # Y[u, i] = nb success at time t_u in profile u (n*I),
   # R[u, i] = nb trials at time t_u in profile i (n*I)
   # diffTime[u] = t_(u+1) - t_u
   # RateInit = initial guess for the transition rates (1->2 and 2->1)
   # Rate = transition rates (same for all profiles)
   # gamma = succes probability in each state  (same for all profiles)
   # Init
   n = nrow(Y); I = ncol(Y); K = 2
   Tmax = sum(diffTime)
   diffTime = diffTime/Tmax
   Tau = TauInit; Eta = list()
   invisible(sapply(1:I, function(i){
      Eta[[i]] <<- t(sapply(2:n, function(t){as.vector(outer(Tau[[i]][t-1,], Tau[[i]][t,]))}))
   }))
   
   diff = 2*tolCvg; iter = 0; Rate = RateInit*Tmax
   while((diff > tolCvg) & (iter < iterMax)){
      iter = iter+1
      
      # M-step
      cat('Probas & Phi ')
      Ytau = Rtau = rep(0, K);
      invisible(sapply(1:I, function(i){
         Ytau <<- Ytau + t(Tau[[i]])%*%Y[, i]
         Rtau <<- Rtau + t(Tau[[i]])%*%R[, i]
      }))
      gamma = (as.vector(Ytau / Rtau))
      logPhi = list()
      invisible(sapply(1:I, function(i){
         logPhi[[i]] <<- matrix(0, n, K)
         sapply(1:K, function(k){logPhi[[i]][, k] <<- dbinom(Y[, i], R[, i], gamma[k], log=T)})
         logPhi[[i]][which(logPhi[[i]] < -300)] <<- -300
      }))
      cat('Rates \n')
      Tau1 = lapply(1:I, function(i){Tau[[i]][1, ]})
      Rate = MstepRateJointBinContHMM(Rate, diffTime, Tau1, Eta)$Rate
      
      # E step
      PiNu = F_ComputeBinaryPiNu(Rate, diffTime)
      cat('Forward-Backward ')
      FW = BW = list(); logL = 0
      invisible(sapply(1:I, function(i){
         FW[[i]] <<- F_HetForward(PiNu$Pi, PiNu$Nu, logPhi[[i]])
         BW[[i]] <<- F_HetBackward(PiNu$Pi, FW[[i]]$Fw)
         logL <<- logL + FW[[i]]$logL
      }))
      
      # Test & update
      cat('Test \n')
      diff = max(unlist(lapply(1:I, function(i){max(abs(Tau[[i]] - BW[[i]]$Tau))})))
      invisible(sapply(1:I, function(i){Tau[[i]] <<- BW[[i]]$Tau; Eta[[i]] <<- BW[[i]]$Eta}))
      cat(iter, gamma, Rate/Tmax, logL, diff, '\n')
   }
   Rate = Rate/Tmax
   
   return(list(gamma=gamma, Rate=Rate, Tau=Tau, Eta=Eta, logL=logL,
               logPhi=logPhi, Pi = PiNu$Pi, Nu = PiNu$Nu))
}


