# EM algo for the MLE of a (BINARY) continuous time Markov process
# lambda : 0 -> 1; mu : 1 -> 0; Q = (-\lambda \lambda; \mu -\mu)
# Analysis of the 4 tetrades independetly

rm(list=ls())
source('FunctionsContTimeHMM.R'); 
source('FunctionsBinaryContTimeHMM.R'); 
source('FunctionsHMMGeneral.R')
source('FunctionsHMMBinomial.R')
# library(expm); 
library(mclust)

# Parms
DataDir = '../Data/'; ResDir = '../Res/'; 
ColorCode = c(4, 1)
DataList = c('tetrade1', 'tetrade2LUHQ1bis', 'tetrade3bis', 'tetrade3', 'tetrade4')

##############################
# boucle sur toutes les tetrades
for (DataName in DataList){

# Analyses
load(paste(DataDir, DataName, '.Rdata', sep=''))
tetrade = get(paste(DataName))
M.list = names(tetrade); I = length(M.list)
for (i in 1:I){
   M = M.list[i]
  data.M = tetrade[[i]]
  Chr.list = levels(data.M$Chrom)
  for (Chr in  Chr.list){
        ######################################################
    # Data extraction
    cat('******************************************************\n')
    cat(DataName, M, Chr, '\n')
    cat('******************************************************\n')
    DataSubName = paste(DataName, '-', M, '-', Chr, sep='')
    Y = data.M$Y; Y = Y[which(data.M$Chrom==Chr)]
    R = data.M$R; R = R[which(data.M$Chrom==Chr)]
    Pos = data.M$Pos; Pos = Pos[which(data.M$Chrom==Chr)]
    # Data reduction
    pdf(paste(ResDir, DataSubName, '-IndCHMM-Data.pdf', sep=''))
    par(pch=20, mfrow=c(2, 1), mex=.6)
    # Interval = (1:20000); Y = Y[Interval]; R = R[Interval]; Pos = Pos[Interval]
    diffPos = diff(Pos); n = length(Y)
    # Plot data
    plot(Pos, Y/R, xlab='', ylab='', main = paste(M, Chr)); 
    hist(log10(diffPos), breaks=sqrt(n), main='', ylab='')
    dev.off()
    
    ######################################################
    # EM discrete time
    cat('\nEM for discrete time HMM \n')
    pdf(paste(ResDir, DataSubName, '-IndCHMM-Classif.pdf', sep=''))
    par(mfrow=c(3, 1), pch=20, mex=.6)
    EMdisc = FitHMMBinom1D(Y, R, K=2)
    cat('Transitions:', EMdisc$Pi[1, 2], EMdisc$Pi[2, 1], '\n')
    cat('Proba:', EMdisc$gamma, '\n')
    EMdisc$MAP = apply(EMdisc$tau, 1, which.max);
    EMdisc$Vit = Viterbi(EMdisc$Phi, EMdisc$nu, EMdisc$Pi)$Z
    plot(Pos, Y/R, col=ColorCode[EMdisc$MAP], xlab='', ylab='')
    
    ######################################################
    # Init EM continous time
    cat('\nFast approximation of the continuous time HMM \n')
    TransProba = c(EMdisc$Pi[1, 2], EMdisc$Pi[2, 1])
    RateInit = InitRateFromTransProba(TransProba, diffPos)
    cat(RateInit, '\n')
    EMfast = list(); EMfast$Rate = RateInit
    
    ######################################################
    # Continous time classif without parm update
    PiNu = F_ComputeBinaryPiNu(EMfast$Rate, diffPos)
    EMfast$Pi = PiNu$Pi; EMfast$Nu = PiNu$Nu
    FW = F_HetForward(PiNu$Pi, PiNu$Nu, log(EMdisc$Phi))
    BW = F_HetBackward(PiNu$Pi, FW$Fw)
    EMfast$Tau = BW$Tau; EMfast$logL = FW$logL
    EMfast$MAP = apply(BW$Tau, 1, which.max); 
    EMfast$Vit = F_HetViterbi(EMfast$Pi, EMfast$Nu, log(EMdisc$Phi))
    plot(Pos, Y/R, col=ColorCode[EMfast$MAP], xlab='', ylab='')
    table(EMdisc$MAP, EMfast$MAP); MAPdiff = which(EMfast$MAP != EMdisc$MAP)
    points(Pos[MAPdiff], Y[MAPdiff]/R[MAPdiff], pch=21, col='red')
    
    ######################################################
    # EM continous time
    cat('\nEM for continuous time HMM \n')
    EMcont = F_EMBinContHMMBinom1D(Y, R, diff(Pos), RateInit)
    cat('Rates:', EMcont$Rate, '\n') 
    cat('Proba:', EMcont$gamma, '\n')
    EMcont$MAP = apply(EMcont$Tau, 1, which.max); 
    EMcont$Vit = F_HetViterbi(EMcont$Pi, EMcont$Nu, EMcont$logPhi)
    
    plot(Pos, Y/R, col=ColorCode[EMcont$MAP], xlab='', ylab='')
    table(EMcont$MAP, EMdisc$MAP); MAPdiff = which(EMcont$MAP != EMdisc$MAP)
    points(Pos[MAPdiff], Y[MAPdiff]/R[MAPdiff], pch=21, col='red')
    table(EMcont$MAP, EMfast$MAP)
    dev.off()
    
    Res = list(EMdisc=  EMdisc, EMfast=EMfast, EMcont=EMcont)
    save(Res, file=paste(ResDir, DataSubName, '-IndCHMM-Results.Rdata', sep=''))
  }
}

##############################
# fin boucle sur toutes les tetrades
}
   