# EM algo for the MLE of a (BINARY) continuous time Markov process
# lambda : 0 -> 1; mu : 1 -> 0; Q = (-\lambda \lambda; \mu -\mu)
# Analysis of the 4 tetrades independetly

rm(list=ls())
source('/home/dcharif/save/bioinfo-dev/r-dev/hmm-nco/functions/FunctionsContTimeHMM.R'); 
source('/home/dcharif/save/bioinfo-dev/r-dev/hmm-nco/functions/FunctionsBinaryContTimeHMM.R'); 
source('/home/dcharif/save/bioinfo-dev/r-dev/hmm-nco/functions/FunctionsHMMGeneral.R')
source('/home/dcharif/save/bioinfo-dev/r-dev/hmm-nco/functions/FunctionsHMMBinomial.R')

library(mclust)
# Parameters
args <- commandArgs(trailingOnly=TRUE)


DataDir = args[1]	
ResDir = args[2]
DataName = args[3]

ColorCode = c(4, 1)

   # Fit
load(paste(DataDir,"/", DataName, '.Rdata', sep=''))
tetrade = get(paste(DataName))
data.M = tetrade[[1]]
Chr.list = levels(data.M$Chrom); Chr.nb = length(Chr.list)
M.list = names(tetrade); I = length(M.list)
for (c in  1:Chr.nb){
    Chr = Chr.list[c]
   #Chr = Chr.lisChr.nb-c+1]
  cat('******************************************************\n')
  cat(DataName, Chr, ': initialization \n')
  Y = c(); R = c(); RateInit = rep(0, 2)
  for (i in 1:I){
    # Data & Init
    load(paste(DataDir,"/", DataName, '.Rdata', sep=''))
    data.M = tetrade[[i]]
    Yi = data.M$Y; Y = cbind(Y, Yi[which(data.M$Chrom==Chr)])
    Ri = data.M$R; R = cbind(R, Ri[which(data.M$Chrom==Chr)])
    Pos = data.M$Pos; Pos = Pos[which(data.M$Chrom==Chr)]
  }
  # Data reduction
  # n = 20000; Y = Y[1:n, ]; R = R[1:n, ]; Pos = Pos[1:n]
  diffPos = diff(Pos); 

  # Init
  RateInit = TauInit = list()
  for (i in 1:I){
    cat('EMdisc', i, '\n')
    EMdisc = FitHMMBinom1D(Y[, i], R[, i], K=2)
    TauInit[[i]] = EMdisc$tau
    RateInit[[i]] = InitRateFromTransProba(c(EMdisc$Pi[1, 2], EMdisc$Pi[2, 1]), diffPos)
  }
  print(unlist(lapply(1:I, function(i){RateInit[[i]]})))
  CommonRateInit = matrix(0, I, 2)
  invisible(sapply(1:I, function(i){CommonRateInit[i, ] <<- RateInit[[i]]}))
  CommonRateInit = exp(apply(log(CommonRateInit), 2, median))
  
  # EM diff rates  cat('******************************************************\n')
  #cat(DataName, Chr, ': different rates:\n')
  #EMdiff = F_EMJointBinContHMMBinom1D_DiffRates(Y, R, diffPos, RateInit, TauInit)
  #EMdiff$MAP = lapply(1:I, function(i){apply(EMdiff$Tau[[i]], 1, which.max)})
  #EMdiff$Vit = lapply(1:I, function(i){F_HetViterbi(EMdiff$Pi[[i]], EMdiff$Nu[[i]], EMdiff$logPhi[[i]])})
  #save(EMdiff, file=paste(ResDir, DataName, '-', Chr, '-JointCHMM-DiffRates.Rdata', sep=''))
  
  # EM same rate
  cat('******************************************************\n')
  cat(DataName, Chr, ': same rates:\n')
  EMsame = F_EMJointBinContHMMBinom1D_SameRates(Y, R, diffPos, CommonRateInit, TauInit)
  EMsame$MAP = lapply(1:I, function(i){apply(EMsame$Tau[[i]], 1, which.max)})
  EMsame$Vit = lapply(1:I, function(i){F_HetViterbi(EMsame$Pi, EMsame$Nu, EMsame$logPhi[[i]])})
  save(EMsame, file=paste(ResDir, "/",DataName, '-', Chr, '-JointCHMM-SameRates.Rdata', sep=''))
}

