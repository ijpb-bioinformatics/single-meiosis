########################################################
StatDistGenerator <- function(Q){
  # Stationary distribution of a cont. time MC with generator Q
  EigtQ = eigen(t(Q))
  Nu = EigtQ$vectors[, which.min(abs(EigtQ$values))]
  Nu = Nu / sum(Nu)
  return(Nu)
}

########################################################
SimContTimeMC <- function(Q, N){
  # Simulates a stationary continuous time MC until N events are observed
  K = nrow(Q)
  nu = StatDistGenerator(Q)
  
  # State simulation
  P = Q - diag(diag(Q)); P = P / rowSums(P)
  State = SimMC(N, nu, P) %*% (1:K)
  
  # Time simulation
  Time = rep(0, N)
  for (i in 2:N){
    Time[i] = Time[i-1] + rexp(1, -Q[State[i-1], State[i-1]])
    }

  return(list(State=State, Time=Time))
}

########################################################
MstepGauss <- function(Y, Tau, modelNames='E'){
  # M step of the estimation of the emission parameters
  
  K = ncol(Tau)
  N = colSums(Tau)
  gamma = Y%*%Tau / N
  sigma = Y^2%*%Tau / N - gamma^2
  if(modelNames=='E'){
    sigma = rep(sqrt(sum(N*sigma) / n), K)
  }else{
    sigma = sqrt(sigma)
  }
  logPhi = sapply(1:K, function(k){dnorm(Y, mean=gamma[k], sd=sigma[k], log=T)})
  logPhi[which(logPhi < -300)] = -300
  
  return(list(gamma=gamma, sigma=sigma, logPhi=logPhi))
}

########################################################
MstepBinom1D <- function(Y, R, Tau){
  # M step of the estimation of the emission parameters
  
  K = ncol(Tau)
  gamma = (t(Tau)%*%Y) / (t(Tau)%*%R)
  logPhi = sapply(1:K, function(k){dbinom(Y, R, prob=gamma[k], log=T)})
  logPhi[which(logPhi < -300)] = -300
  
  # return(list(gamma=gamma, sigma=sigma, logPhi=logPhi))
  return(list(gamma=gamma, logPhi=logPhi))
}

########################################################
F_HetForward <- function(Pi, Nu, logPhi)
{
  # Forward recursion for heterogeneous HMM
  # Nu = initial distribution : 1 x k
  # Pi = matrix of transition matrices : (n-1) x K^2
  # such that Pi(t) = matrix(Pi[t, ], K, K)
  # logPhi = log-emission densities : n x k
  n = nrow(logPhi); K = ncol(logPhi)
  Phi = exp(logPhi)
  Fw = matrix(0, n, K)
  # if (n > 1e4){pr.step = round(sqrt(n))}else{pr.step = n+1}
  
  logL = 0
  Fw[1,] = Phi[1,] * Nu
  logL = logL + log(sum(Fw[1,]))
  Fw[1,] = Fw[1,] / sum(Fw[1,])
  for (t in (2:n)){
    # if(t %% pr.step==0){cat(t, '')}
    Fw[t,] = Fw[t-1,] %*% matrix(Pi[t-1, ], K, K)
    Fw[t,] = Fw[t,] * Phi[t,]
    logL = logL + log(sum(Fw[t,]))
    Fw[t,] = Fw[t,] / sum(Fw[t,])
  }
  # if (n > 1e6){cat('\n')}
  return(list(Fw=Fw, logL=logL))
}

########################################################
F_HetBackward <- function(Pi, Fw)
{
  # Eta = matrix with same structure as Pi
  n = dim(Fw)[1]
  K = dim(Fw)[2]
  # Tau(t,k) = E(Z(t,k) | X)
  Tau = matrix(0, n, K)
  G = Tau
  Tau[n,] = Fw[n, ]
  # Eta[t, ] = E(Z(t, k) Z(t+1, l) | X)
  Eta = matrix(0, (n-1), K^2)
  for (t in ((n-1):1)){
    Pi.t = matrix(Pi[t, ], K, K)
    G[t+1, ] = Fw[t,] %*% Pi.t
    B = as.vector(Tau[t+1, ] / G[t+1, ])
    Tau[t,] = Fw[t,] * (Pi.t %*% B)
    Tau[t,] = Tau[t,] / sum(Tau[t,])
    Eta[t, ] = Pi[t, ] * as.vector(outer(Fw[t,], B))
  }
  return(list(Tau=Tau, Eta=Eta, G=G))
}

########################################################
F_HetViterbi <- function(Pi, Nu, logPhi)
{
  # Viterbi recursion for heterogeneous HMM
  # Nu = initial distribution : 1 x k
  # Pi = matrix of transition matrices : (n-1) x K^2
  # such that Pi(t) = matrix(Pi[t, ], K, K)
  # logPhi = log-emission densities : n x k
  n = nrow(logPhi); K = ncol(logPhi); Phi = exp(logPhi)
  
  # Forward recursion
  V = S = matrix(0, n, K)
  V[1,] = Phi[1,] * Nu; V[1,] = V[1,] / sum(V[1,])
  for (t in (2:n)){
    # if(t %% pr.step==0){cat(t, '')}
    Mtmp = outer(rep(1, K), V[t-1, ]) * matrix(Pi[t-1, ], K, K)
    V[t,] = apply(Mtmp, 1, max) * Phi[t, ]
    V[t,] = V[t,] / sum(V[t,])
    S[t-1,] = apply(Mtmp, 1, which.max)
  }

  # Backward recursion
  Zhat = rep(0, n)
  Zhat[n] = which.max(V[n, ])
  invisible(sapply((n-1):1, function(t){Zhat[t] <<- S[t, Zhat[t+1]]}))
  return(Zhat)
}

########################################################
F_MergePhiPiNu <- function(EMres.list){
  # Merges the results of I independent K-state cont. time HMM 
  # into a single cont. time HMM with K^I states
  I = length(EMres.list)
  n = nrow(EMres.list[[1]]$logPhi)
  K = ncol(EMres.list[[1]]$logPhi)
  State = expand.grid((1:K), (1:K), (1:K), (1:K))
  S = nrow(State)
  
  # Merge and shape logPhi and Nu
  logPhi = matrix(0, n, S); Nu = rep(1, S); 
  invisible(sapply(1:S, function(s){
    sapply(1:I, function(i){
      logPhi[, s] <<- logPhi[, s]+ EMres.list[[i]]$logPhi[, State[s, i]]
      Nu[s] <<- Nu[s] * EMres.list[[i]]$Nu[State[s, i]]
    })
  }))
  logPhi[which(logPhi < -300)] = -300
  
  # Shape each EMres.list[[i]]$Pi as a (n-1)*K*K array
  Pi.mat.list = list()
  invisible(lapply(1:I, function(i){
    Pi.mat.list[[i]] <<- array(dim=c((n-1), K, K))
    sapply(1:K, function(k){
      Pi.mat.list[[i]][, k, ] <<- EMres.list[[i]]$Pi[, k+K*(0:(K-1))]
    })}))
  # Compute the enlarged transition matrices as a (n-1)*S*S array
  Pi.mat = array(dim=c((n-1), S, S))
  invisible(sapply(1:S, function(s){sapply(1:S, function(r){
    Pi.mat[, s, r] <<- 1
    sapply(1:I, function(i){
      Pi.mat[, s, r] <<- Pi.mat[, s, r] * Pi.mat.list[[i]][, State[s, i], State[r, i]]
    })
  })}))
  Pi = matrix(0, (n-1), S^2)
  invisible(sapply(1:S, function(s){sapply(1:S, function(r){
    Pi[, s + S*(r-1)] <<- Pi.mat[, s, r]
  })}))
  invisible(list(Nu=Nu, Pi=Pi, logPhi=logPhi))
}


########################################################
F_ClassifS <- function(Parms.list){
  # Classification of the joint hidden state S = (Z1, Z2,  ... ZI)
  # ParmList = list of I element, each being a list(logPhi, Pi, Nu)
  
  cat('Check results: Merge parms ')
  PhiPiNu = F_MergePhiPiNu(Parms.list)
  PhiPiNu$logPhi[which(PhiPiNu$logPhi < -300)] = -300
  cat('MAP ')
  JointFW = F_HetForward(PhiPiNu$Pi, PhiPiNu$Nu, PhiPiNu$logPhi)
  JointBW = F_HetBackward(PhiPiNu$Pi, JointFW$Fw)
  JointMAP = apply(JointBW$Tau, 1, which.max);
  cat('Viterbi ')
  JointVit = F_HetViterbi(PhiPiNu$Pi, PhiPiNu$Nu, PhiPiNu$logPhi)
  
  return(list(MAPS = JointMAP, TauS = JointBW$Tau, VitS = JointVit, logL = JointFW$logL))
}


########################################################
F_VitMaskS <- function(Parms.list, Mask.list){
  # Classification of the joint hidden state S = (Z1, Z2,  ... ZI)
  # ParmList = list of I element, each being a list(logPhi, Pi, Nu)
  # Mask.list = list of matrix of allowed transtions
  
  cat('Classif with masks: Merge parms ')
  PhiPiNu = F_MergePhiPiNu(Parms.list)
  S = sqrt(ncol(PhiPiNu$Pi)); n1 = nrow(PhiPiNu$Pi)
  MaskClassif = list()
  for (d in 1:length(Mask.list)){
    cat(d, '')
    Mask.vec = as.vector(Mask.list[[d]])
    Pi = PhiPiNu$Pi * outer(rep(1, n1), Mask.vec)
    invisible(sapply(1:K, function(k){
      k.vec = k + K*(0:(K-1))
      Pi[, k.vec] <<- Pi[, k.vec] / rowSums(Pi[, k.vec])
    }))
    cat('MAP ')
    MaskFW = F_HetForward(Pi, PhiPiNu$Nu, PhiPiNu$logPhi)
    logL = MaskFW$logL
    MaskBW = F_HetBackward(Pi, MaskFW$Fw)
    TauS = MaskBW$Tau
    MAPS = apply(MaskBW$Tau, 1, which.max);
    cat('Viterbi ')
    VitS = F_HetViterbi(Pi, PhiPiNu$Nu, PhiPiNu$logPhi)
    MaskClassif[[d]] = list(MAPS=MAPS, VitS=VitS, TauS=TauS, logL=logL)
    cat(logL ,'')
  }
  return(MaskClassif)
}
