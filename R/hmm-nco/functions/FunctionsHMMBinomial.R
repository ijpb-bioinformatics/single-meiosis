EstepHMM <- function(Mstep){
  FB = ForBack(Mstep$Phi, Mstep$nu, Mstep$Pi, nrow(Mstep$Phi), ncol(Mstep$Phi))
  return(list(Tau = FB$Tau, logL = FB$logL, Eta = FB$Eta))
}
 
CalPhiBinom1D <- function(Y, R, gamma){
  Phi = matrix(0, length(Y), length(gamma))
  for (k in 1:length(gamma)){Phi[, k] = dbinom(Y, R, prob=gamma[k])}
  # Modif 22/9/16 : plantage si Phi = 0
  Phi[which(Phi==0)] = exp(-300)
  return(Phi)
}

MstepHMMBinom1D <- function(Y, R, Estep){
  Pi = Estep$Eta; Pi = Pi / rowSums(Pi)
  nu = StatDist(Pi); K = length(nu)
  gamma = (t(Estep$Tau)%*%Y) / (t(Estep$Tau)%*%R)
  order = order(gamma)
  gamma = gamma[order]; Pi = Pi[order, order];  nu = nu[order]
  return(list(Pi = Pi, nu = nu, gamma = gamma, Phi = CalPhiBinom1D(Y, R, gamma)))
}

FitHMMBinom1D <- function(Y, R, K, tolCvg=1e-5, epsTau=1e-4, iterMax=100)
{
  # R = nb trials, Y = nb success
  n = length(Y)
  # Init
  library(mclust); Estep = list()
  GMM = Mclust(asin(sqrt((Y+.5)/(R+1))), G=K, modelNames='E')
  Estep$Tau = GMM$z
  Estep$Eta = t(Estep$Tau[1:(n-1), ]) %*% Estep$Tau[2:n, ]
  Mstep = MstepHMMBinom1D(Y, R, Estep)
  Theta = c(as.vector(Mstep$Pi), Mstep$gamma)
  
  # E-M
  diff = 2*tolCvg
  while(diff > tolCvg){
    
    # Accelated EM
    Estep = EstepHMM(Mstep)
    Mstep = MstepHMMBinom1D(Y, R, Estep)
    # Theta1 = c(as.vector(Mstep$Pi), Mstep$gamma); # cat(Theta1)
    # Estep = EstepHMM(Mstep)
    # Mstep = MstepHMMBinom1D(Y, R, Estep)
    # Theta2 = c(as.vector(Mstep$Pi), Mstep$gamma); # cat(Theta2)
    # alpha = + sqrt(sum(Theta1 - Theta)^2) / sqrt(sum(Theta2 - 2*Theta1 + Theta)^2)
    # Theta.new = Theta + alpha*(Theta1 - Theta)
    # while((min(Theta.new)<0) | (max(Theta.new)>1)){alpha = alpha/2; Theta.new = Theta + alpha*(Theta1 - Theta)}
    # Mstep$Pi = matrix(Theta.new[1:K^2], K, K); Mstep$nu = StatDist(Mstep$Pi)
    # Mstep$gamma = Theta.new[(K^2+1):(K^2+K)]; Mstep$Phi = CalPhiBinom1D(Y, R, Mstep$gamma)
    Theta.new = c(as.vector(Mstep$Pi), Mstep$gamma); 
    
    # Test & mise a jour
    diff = max(abs(Theta.new - Theta))
    Theta = Theta.new; # cat(Theta)
    
    cat(as.vector(Mstep$Pi), Mstep$gamma, Estep$logL, diff, '\n')
  }

  HMM = list(tau=Estep$Tau, Pi=Mstep$Pi, nu=Mstep$nu, gamma=Mstep$gamma, Phi=Mstep$Phi, logL=Estep$logL)
  return(HMM)
}
  