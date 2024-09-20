#library(HMM)

#######################################################
# Forward recursion
Forward <- function(Phi, nu, Pi)
{
  n = dim(Phi)[1]
  K = dim(Phi)[2]
  Fw = matrix(0, n, K)
  
  logL = 0
  Fw[1,] = Phi[1,] * nu
  logL = logL + log(sum(Fw[1,]))
  Fw[1,] = Fw[1,] / sum(Fw[1,])
  for (t in (2:n)){
    Fw[t,] = Fw[t-1,] %*% Pi
    Fw[t,] = Fw[t,] * Phi[t,]
    logL = logL + log(sum(Fw[t,]))
    Fw[t,] = Fw[t,] / sum(Fw[t,])
  }
  return(list(Fw=Fw, logL=logL))
}

#######################################################
# Backward recursion
Backward <- function(Fw, Pi, epsTau = 1e-6)
{
  n = dim(Fw)[1]
  K = dim(Fw)[2]
  Fw = Fw + epsTau; Fw = Fw / rowSums(Fw)
  # Tau(t,k) = E(Z(t,k) | X)
  Tau = matrix(0, n, K)
  G = Tau
  Tau[n,] = Fw[n, ]
  # Eta(k, l) = sum_t E(Z(t, k) Z(t+1, l) | X)
  Eta = matrix(0, K, K)
  for (t in ((n-1):1)){
    # t = t-1
    G[t+1, ] = Fw[t,] %*% Pi
    B = as.vector(Tau[t+1, ] / G[t+1, ])
    Tau[t,] = Fw[t,] * (Pi %*% B); Tau[t,] = Tau[t,] / sum(Tau[t,])
    Tau[t,] = Tau[t,] + epsTau; Tau[t,] = Tau[t,] / sum(Tau[t,])
    Eta = Eta + Pi * outer(Fw[t,], B)
    # cat(t, ':', Tau[t, ], '\n')
  }
  return(list(Tau=Tau, Eta=Eta, G=G))
}

#######################################################
# Viterbi path
Viterbi <- function(Phi, nu, Pi)
{
  n = dim(Phi)[1]
  K = dim(Phi)[2]
  # Prob of the most probable path
  Lik = matrix(0, n, K)
  TransMax = Lik
  Lik[1,] = Phi[1,] * nu
  logLikTot = log(sum(Lik[1, ]))
  Lik[1, ] = Lik[1, ] / sum(Lik[1, ])
  for (t in (2:n)){
    LikTrans = matrix(rep(Lik[t-1,], each=K), K, K) * Pi
    Lik[t,] = Phi[t,] * apply(LikTrans, 1, max)
    logLikTot = logLikTot + log(sum(Lik[t, ]))
    Lik[t, ] = Lik[t, ] / sum(Lik[t, ])
    TransMax[t, ] = apply(LikTrans, 1, which.max)   
  }
  # Most probable path
  Z = rep(0, n)
  Z[n] = which.max(Lik[n,])
  for (t in ((n-1):1)){
    Z[t] = TransMax[(t+1), Z[t+1]]
  }
  LikTot = exp(logLikTot)
  
  return(list(Z=Z, logLik=logLikTot))
}

#######################################################
# Stationary distribution
StatDist <- function(Pi)
{
  Eig = eigen(t(Pi))
  k = which(abs(Eig$values-1) < 1e-8)
  nu = Re(as.vector(Eig$vectors[,k]))
  nu = nu / sum(nu)
  return(nu)
}

#######################################################
# Simulation of Markov chain
SimMC <- function(n, nu, Pi)
{
  K = length(nu)
  Z = matrix(0, n, K)
  Z[1, ] = t(rmultinom(1, 1, nu))
  for (t in (2:n)){Z[t, ] = t(rmultinom(1, 1, Pi[Z[(t-1), ]%*%(1:K), ]))}
  return(Z)
}

#######################################################
# Include C function 'Forward'
ForwardR = function(Phi.vect,muHMM,Mat.trans.norm.vect,n,K)
{
  .C("Forward",as.double(Phi.vect),as.double(muHMM),as.double(Mat.trans.norm.vect),as.integer(n),as.integer(K),F=double(n*K), PACKAGE="HMM")

}

#######################################################
# Include C function 'Backward'
BackwardR = function(F,Mat.trans.norm.vect,n,K)
{
  .C("Backward", as.double(F),as.double(Mat.trans.norm.vect),as.integer(n),as.integer(K),tau=double(n*K),G=double(n*K), PACKAGE="HMM")
}

#######################################################
# ForwardR, BackwardR outputs into 'HMM' format
ForBack = function(Phi, nu, Pi, n, K, epsTau=1e-4)
{
  # C version
#   Fw = ForwardR(as.vector(Phi), as.vector(nu), as.vector(Pi), n, K)
#   Bw = BackwardR(Fw$F, as.vector(Pi), n, K)
#   Fw = matrix(Fw$F, n, K, byrow=F)
#   Tau = matrix(Bw$tau, n, K, byrow=F)
#   G = matrix(Bw$G, n, K, byrow=F)
#   Eta = Pi * (t(Fw[-n,])%*% (Tau[-1,] / G[-1,]))
#   logL = NaN
  
  # R version
  Fw = Forward(Phi, nu, Pi)
  Bw = Backward(Fw$Fw, Pi)
  logL = Fw$logL
  Fw = Fw$Fw
  Tau = Bw$Tau
  G = Bw$G
  Eta = Bw$Eta

  # Smoothing
  Tau = Tau + epsTau;   Tau = Tau / rowSums(Tau)
  
  return(list(Tau=Tau, Fw=Fw, G=G, Eta=Eta, logL=logL))
}
