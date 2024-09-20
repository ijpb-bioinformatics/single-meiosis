########################################################
LikelihoodBinContHMM <- function(Rate, diffTime, State){
  # Log-likelihood of a binary cont. time MC observed a sequence of times
  # EigQ = eigen(Q); U = solve(t(EigQ$vectors)); diagD = EigQ$values; D = diag(diagD); V = t(EigQ$vectors)
  # tau = lambda+mu; diagD = c(0, -tau); D = diag(diagD); V = rbind(c(1, 1), c(-lambda, mu)); U = rbind(c(mu, -1), c(lambda, 1))/tau
  # EigtQ = eigen(t(Q)); U = solve(t(EigtQ$vectors)); diagD = EigtQ$values; D = diag(diagD); V = t(EigtQ$vectors)
  # EigtQ = eigen(t(Q)); U = EigtQ$vectors; diagD = EigtQ$values; D = diag(diagD); V = solve(U)
  lambda = Rate[1]; mu = Rate[2]; 
  tau = lambda+mu; a = lambda/tau; b = 1-a; Nu = c(b, a)
  Pcst = matrix(c(a, a, b, b), 2, 2)
  Pcoef = matrix(c(b, -a, -b, a), 2, 2)
  logL = log(Nu[State[1]]); 
  invisible(sapply(2:(length(diffTime)+1), function(i){
    logL <<- logL + log(Pcst[State[i-1], State[i]] + Pcoef[State[i-1], State[i]]*exp(-tau*diffTime[i-1]))
  }))
  return(logL)
}

########################################################
CondExpLikelihoodBinContHMM <- function(Rate, diffTime, Tau1, Eta){
  # Conditional expectation of the complete log-likelihood
  # of a binary cont. time MC observed a sequence of times
  # diffTime = sequence of interval observations times : n x 1
  # Eta = cond. probability for each possible intial 
  # state (0, 1) : 1 x 2
  # Eta = cond. probability for each possible pair 
  # of successive states (00, 10, 01, 11) : n x 4
  lambda = Rate[1]; mu = Rate[2]; 
  tau = lambda+mu; a = lambda/tau; b = 1-a; Nu = c(b, a)
  Pcst = c(a, a, b, b)
  Pcoef = c(b, -a, -b, a)
  logL = Tau1 %*% log(Nu); 
  invisible(sapply(1:(length(diffTime)), function(t){
    logL <<- logL + Eta[t, ] %*% log(Pcst + Pcoef*exp(-tau*diffTime[t]))
    }))
  return(logL)
}

########################################################
InitRateFromTransProba <- function(TransProb, diffTime){
  # Initial guess for the transition rates based on the transition proba from 
  # the discrete time HMM
  
  sumRate = -log(1 - sum(TransProb))/mean(diffTime)
  stop = F
  while(stop==F){
    RateInit = sumRate*TransProb/(1 - mean(exp(-sumRate*diffTime)))
    sumRate.new = sum(RateInit)
    if(abs(sumRate.new-sumRate)/sumRate < 1e-4){stop=T}
    sumRate = sumRate.new
  }
  return(RateInit)
}

########################################################
MstepRateBinContHMM <- function(RateInit, diffTime, Tau1, Eta){
  # M step of the estimation of the rates of an hidden continuous time Markov model
  
  negLogL <- function(Rate, diffTime, Tau1, Eta){
    -CondExpLikelihoodBinContHMM(Rate, diffTime, Tau1, Eta)
  }
  Opt = optim(RateInit, negLogL, gr=NULL, diffTime, Tau1, Eta, method='BFGS', control=list(trace=2))
  
  return(list(Rates=Opt$par))
}

########################################################
MstepRateJointBinContHMM <- function(RateInit, diffTime, Tau1, Eta){
  # M step of the estimation of the rates of an hidden continuous time Markov model
  
  negLogL <- function(Rate, diffTime, Tau1, Eta){
    L = 0; I = length(Tau1)
    invisible(sapply(1:I, function(i){
      L <<- L - CondExpLikelihoodBinContHMM(Rate, diffTime, Tau1[[i]], Eta[[i]])
    }))
    return(L)
  }
  Opt = optim(RateInit, negLogL, gr=NULL, diffTime, Tau1, Eta, method='BFGS', control=list(trace=2))
  
  return(list(Rates=Opt$par))
}

########################################################
F_ComputeBinaryPiNu <- function(Rate, diffTime){
  # Calculation of the initial dist. and heterogenous transition natrices
  lambda = Rate[1]; mu = Rate[2]; 
  tau = lambda+mu; a = lambda/tau; b = 1-a; 
  Nu = c(b, a)
  Pcst = c(a, a, b, b)
  Pcoef = c(b, -a, -b, a)
  Pi = outer(rep(1, length(diffTime)), Pcst) + outer(exp(-tau*diffTime), Pcoef)
  return(list(Pi=Pi, Nu=Nu))
}


source('/home/dcharif/save/bioinfo-dev/r-dev/hmm-nco/functions/FunctionsBinaryContTimeHMM-Binom1D.R')
source('/home/dcharif/save/bioinfo-dev/r-dev/hmm-nco/functions/FunctionsBinaryContTimeHMM-BetaBinom.R')

