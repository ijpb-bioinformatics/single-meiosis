# Merge results from various models in view of marke classification

rm(list=ls())

args <- commandArgs(trailingOnly=TRUE)


DataDir = args[1]
ResDir = args[2]
DataName = args[3]

ColorCode = c(4, 1)


# States
K = 2; I = 4; S = K^I
State = expand.grid((1:K),(1:K),(1:K),(1:K))
StateH = 4 - rowSums(State-1)
NCO = which((StateH==1) | (StateH==3))
DistState = as.matrix(dist(State, diag=T)^2)
Mask.list = lapply(1:(max(DistState)-1), function(d){1*(DistState <= d)})

# Functions
F_Z2S <- function(Zvec){1+sum((Zvec-1)*(2^(0:(length(Zvec)-1))))}

# Plot #  Modif tmp pour tetrade1bis car Chr5 pb
Chr.nb = 5; Chr.list = paste('Chr', (1:Chr.nb), sep='')
I = 4

Model.list = c('SameProbSameRate')
for (c in  1:Chr.nb){
  Chr = Chr.list[c]
  cat(Chr)
  

  # Results
  DataChrName = paste(DataName, '-', Chr, sep='')
  load(paste(ResDir,"/", DataChrName, '-CHMM-SameProbSameRate-Results.Rdata', sep=''))
  #Model.name = 'SameProbSpecRate'; Fit = Res$SameProbSpecRate
  Model.name = 'SameProbSameRate'; Fit = Res$SameProbSameRate
  # Fit = get(paste('Res$', Model.name, sep=''))

  # Z classif
  TauZ.thres = .9999999999999
  TauZ.MAP = sapply(1:I, function(i){apply(Fit$TauZ[[i]], 1, max)})
  MAPZ.strong = which(apply(TauZ.MAP, 1, min) > TauZ.thres)
  
  # NCO classif
  VitNCO = which(is.element(Fit$VitS, NCO))
  MAPNCO = which(is.element(Fit$MAPS, NCO))
  print(table(is.element(Fit$VitS, NCO), is.element(Fit$MAPS, NCO)))
  ProbNCO = rowSums(Fit$TauS[, NCO])
  MAPNCO.strong = intersect(MAPNCO, MAPZ.strong)
  MAPNCO.weak = setdiff(MAPNCO, MAPZ.strong)
  
  # Plots
   pdf(paste(ResDir,"/", DataChrName, '-CHMM-', Model.name, '.pdf', sep=''))
  par(mfrow=c(6, 1), pch=20, mex=.3, cex=.5)
  Xlim = c(min(Res$Pos), max(Res$Pos))
  for (i in 1:I){
     plot(Res$Pos, Res$Y[, i]/Res$R[, i], xlab='', ylab='', xlim=Xlim, ylim=c(0, 1));
     # points(Res$Pos[MAPNCO], Res$Y[MAPNCO, i]/Res$R[MAPNCO, i], pch=21, col=2)
     points(Res$Pos[MAPNCO.strong], Res$Y[MAPNCO.strong, i]/Res$R[MAPNCO.strong, i], pch=21, col=2)
     points(Res$Pos[MAPNCO.weak], Res$Y[MAPNCO.weak, i]/Res$R[MAPNCO.weak, i], pch=21, col=4)
     abline(h=Fit$Prob, col=3)
  }
  CumRatio = rowSums(Res$Y/Res$R)
  plot(Res$Pos, CumRatio, xlab='', ylab='', xlim=Xlim); 
  points(Res$Pos[MAPNCO.strong], CumRatio[MAPNCO.strong], col=2, pch=21)
  points(Res$Pos[MAPNCO.weak], CumRatio[MAPNCO.weak], col=4, pch=21)
  abline(h=3, col=3) 
  # Attention à Inf
  LProbNCO = -log10(1-ProbNCO)
  LProbNCO[which(LProbNCO==Inf)]=20
  plot(Res$Pos, LProbNCO, type='l')
  # plot(ProbNCO, CumRatio)
  # Write file of Position

   dfRes= data.frame("Chromosome"=rep(Chr,length(Res$Pos)), "NCOPosition"=Res$Pos, "NCOlogProba"=LProbNCO, "NCOProba"=ProbNCO,"MAPNCO"=Fit$MAPS) 
  write.table(dfRes,file=paste(ResDir,"/", DataChrName, '-CHMM-', Model.name, '-MAP-All.txt', sep=''),row.names=FALSE)
  write.table(dfRes[MAPNCO,],file=paste(ResDir,"/", DataChrName, '-CHMM-', Model.name, '-MAP-NCO.txt', sep=''),row.names=FALSE)
  

  
  dev.off()
}
