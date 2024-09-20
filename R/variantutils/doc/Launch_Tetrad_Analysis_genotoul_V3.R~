# Get Args

args <- commandArgs(trailingOnly=TRUE)


#
# Positionnement des variables
#


# Path to the tetrad
TetradPath <- args[1]

# Nom de l'analyse
AnalysisName <- "SCHNEE"

# Fichiers de polymorphismes
GoldVariantPath <- args[2]
sdiFilePath <- args[3]
ColPositionFilePathPref <- args[4]
LerPositionFilePathPref <- args[5]

# List of individu 
# ListIndiv <- c("F1","ParentColLer")
ListIndiv <- c("M1","M2","M3","M4")

VcfSuffix <- "sorted_rmdup_reformated"
Chromosom <- paste("Chr",1:5,sep="")

#
# Chargement des librairies
#

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

# Load libraries

library(ggplot2)
library(stringr)
library(seqinr)
library(data.table)
library(doParallel)
library(VariantAnnotation)

# Source the R code of the package variantutils

sourceDir("/home/dcharif/save/bioinfo-dev/r-dev/variantutils/R")

# Genomes de reference
TAIR10 <- read.fasta("/save/project/ijpb/bank/GENOME/Ath-Col0-tair10-WG/Ath-Col0-tair10-WG.fa")
RefGenomes <- list("TAIR10"=TAIR10)


######################################
# Creation des fichiers au format SM
# -> Fusion des vcf Col et Ler grâce à synthénie
# -> Calcul des couvertures pour les indels
######################################



sdiDT <- fread(sdiFilePath)
sdiDT$ID <- paste(sdiDT$Chrom,sdiDT$colPos,sep="_")
setkey(sdiDT,ID)

# Création des directories output

# Les dossiers ANALYSES, GoldVariant
if(! file.exists(paste(TetradPath,"Analysis",sep="/"))){
    stopifnot(dir.create(paste(TetradPath,"Analysis",sep="/"), showWarnings = TRUE, recursive = FALSE, mode = "0770"))
}

if(! file.exists(paste(TetradPath,"Analysis",basename(GoldVariantPath),sep="/"))){
    stopifnot(dir.create(paste(TetradPath,"Analysis",basename(GoldVariantPath),sep="/"), showWarnings = TRUE, recursive = FALSE, mode = "0770"))
}

if(! file.exists(paste(TetradPath,"Analysis",basename(GoldVariantPath),"SM",sep="/"))){
    stopifnot(dir.create(paste(TetradPath,"Analysis",basename(GoldVariantPath),"SM",sep="/"), showWarnings = TRUE, recursive = FALSE, mode = "0770"))
}

for(i in 1:length(ListIndiv)){
    print(paste("Construct SM format for individu",ListIndiv[i],"for", basename(TetradPath),"and for", basename(GoldVariantPath),sep=" "))

    for(j in 1:length(Chromosom)){

        print(Chromosom[j])
        
        colPosDT <- fread(paste(ColPositionFilePathPref,"_",Chromosom[j],".txt",sep=""))
        setnames(colPosDT, names(colPosDT), c("Chrom","ColPos","ID"))
        stopifnot(file.exists(paste(ColPositionFilePathPref,"_",Chromosom[j],".txt",sep="")))

        colM1DT <- fread(paste(TetradPath,"/",ListIndiv[i],"/",ListIndiv[i],"_Col_",VcfSuffix,"_",Chromosom[j],".vcf",sep="")) 
        setnames(colM1DT, names(colM1DT) ,c("Chrom","ColPos","Ref","DP","colQSUM","colMAPQSUM"))
        stopifnot(file.exists(paste(TetradPath,"/",ListIndiv[i],"/",ListIndiv[i],"_Col_",VcfSuffix,"_",Chromosom[j],".vcf",sep="")))
        
        lerPosDT <- fread(paste(LerPositionFilePathPref,"_",Chromosom[j],".txt",sep=""))
        setnames(lerPosDT,names(lerPosDT),c("Chrom","LerPos","ID"))
        stopifnot(file.exists(paste(LerPositionFilePathPref,"_",Chromosom[j],".txt",sep="")))

        lerM1DT <- fread(paste(TetradPath,"/",ListIndiv[i],"/",ListIndiv[i],"_Ler_",VcfSuffix,"_",Chromosom[j],".vcf",sep=""))
        setnames(lerM1DT, names(lerM1DT), c("Chrom","LerPos","Ref","DP","lerQSUM","lerMAPQSUM"))
        stopifnot(file.exists(paste(TetradPath,"/",ListIndiv[i],"/",ListIndiv[i],"_Ler_",VcfSuffix,"_",Chromosom[j],".vcf",sep="")))
        
        j1DT <- merge(colPosDT,colM1DT,by=c("Chrom","ColPos"),all.x=TRUE)
        j2DT <- merge(lerPosDT,lerM1DT,by=c("Chrom","LerPos"),all.x=TRUE)
        
        setkey(j1DT,"ID")
        setkey(j2DT,"ID")
        
        j1DTMean <- j1DT[, lapply(.SD, function(x){ round(mean(x, na.rm=TRUE),0)} ), by=ID, .SDcols=c("DP","colQSUM","colMAPQSUM") ] 
        j2DTMean <- j2DT[, lapply(.SD, function(x){ round(mean(x, na.rm=TRUE),0)} ), by=ID, .SDcols=c("DP","lerQSUM","lerMAPQSUM") ] 
        
        sm <- merge(j1DTMean, j2DTMean, by=c("ID"),all.x=TRUE) 
        setnames(sm,names(sm),c("ID","colDP","colQsum","colMAPQsum","lerDP","lerQsum","lerMAPQsum"))
        
        
        resultat <- merge(sdiDT,sm,all.y=TRUE)[,c(1:5,7:8,6,12,15,13,16,14,17),with=FALSE]

        # Attention: Remplacer les NA par 0 pour les positions non couvertes
        resultat$lerDP[which(is.na(resultat$lerDP))] <- 0
        resultat$colDP[which(is.na(resultat$colDP))] <- 0

        if(j==1){
            write.table(resultat[order(resultat$Chrom,resultat$colPos),],file=paste(TetradPath,"/","Analysis","/",basename(GoldVariantPath),"/","SM","/",ListIndiv[i],"_","SM.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
        }
        else{
            write.table(resultat[order(resultat$Chrom,resultat$colPos),],file=paste(TetradPath,"/","Analysis","/",basename(GoldVariantPath),"/","SM","/",ListIndiv[i],"_","SM.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
         }
        
    }
}

#############################@
# Analyse de la Tetrade
#############################@


## 1) Convert individuals from SM to VCF format

registerDoParallel(cores=8)

VCFL <- foreach(i=1:4) %dopar% new("VCFlikeParser",Format="SM",
          File=paste(TetradPath,"/","Analysis","/",basename(GoldVariantPath),"/","SM","/",paste("M",i,sep=""),"_","SM.txt",sep=""),
          RefGenome="TAIR10",Sample=paste("M",i,sep=""))

## 1bis) Write Input .Rdata file for hmm 
if(! dir.exists(paste(TetradPath,"Analysis",basename(GoldVariantPath),"InputHMM",sep="/"))){
	stopifnot(dir.create(paste(TetradPath,"Analysis",basename(GoldVariantPath),"InputHMM",sep="/"), showWarnings = TRUE, recursive = FALSE, mode = "0770"))
}
SM2InputHMM(VCFL=VCFL, OutFile=paste(TetradPath,"/","Analysis","/",basename(GoldVariantPath),"/","InputHMM","/",basename(TetradPath),sep="")) 

## 2) Genotype individuals

GVCFL <- lapply(VCFL, function(x){
  x@VCFlike <- Genotype(x@VCFlike,0.003)
  return(x)
})

## 3) Filter individuals

# Get DPmin et DPmax

M1F <- FilterDP(info(GVCFL[[1]]@VCFlike@VCF)$LER_DP[which(geno(GVCFL[[1]]@VCFlike@VCF)$GT == "0/1")],distr=c("Normal"),
         side=("left"),filter=c(0.025))
M2F <- FilterDP(info(GVCFL[[2]]@VCFlike@VCF)$LER_DP[which(geno(GVCFL[[2]]@VCFlike@VCF)$GT == "0/1")],distr=c("Normal"),
                side=("left"),filter=c(0.025))
M3F <- FilterDP(info(GVCFL[[3]]@VCFlike@VCF)$LER_DP[which(geno(GVCFL[[3]]@VCFlike@VCF)$GT == "0/1")],distr=c("Normal"),
                side=("left"),filter=c(0.025))
M4F <- FilterDP(info(GVCFL[[4]]@VCFlike@VCF)$LER_DP[which(geno(GVCFL[[4]]@VCFlike@VCF)$GT == "0/1")],distr=c("Normal"),
                side=("left"),filter=c(0.025))

## Apply the filter

GVCFL[[1]]@VCFlike@VCF <- AnnotateVCF(GVCFL[[1]]@VCFlike@VCF, 
                                      vcfField="LER_DP", 
                                      vcfFieldType=c("info"), 
                                      FilterType=c("numeric"), 
                                      FilterName="DP_M1_Ler_Het_Min", 
                                      FilterValue= M1F$values["min"], 
                                      FilterDescription="DP M1 Het Ler", 
                                      FilterSign="less",  
                                      FilterNA=FALSE,
                                      ConditionalFilterVector=geno(GVCFL[[1]]@VCFlike@VCF)$GT, 
                                      ConditionalFilterValue="0/1",
                                      append=TRUE)

GVCFL[[2]]@VCFlike@VCF <- AnnotateVCF(GVCFL[[2]]@VCFlike@VCF, 
                                      vcfField="LER_DP", 
                                      vcfFieldType=c("info"), 
                                      FilterType=c("numeric"), 
                                      FilterName="DP_M2_Ler_Het_Min", 
                                      FilterValue= M2F$values["min"], 
                                      FilterDescription="DP M2 Het Ler", 
                                      FilterSign="less",  
                                      FilterNA=FALSE,
                                      ConditionalFilterVector=geno(GVCFL[[2]]@VCFlike@VCF)$GT, 
                                      ConditionalFilterValue="0/1",
                                      append=TRUE)


GVCFL[[3]]@VCFlike@VCF <- AnnotateVCF(GVCFL[[3]]@VCFlike@VCF, 
                                      vcfField="LER_DP", 
                                      vcfFieldType=c("info"), 
                                      FilterType=c("numeric"), 
                                      FilterName="DP_M3_Ler_Het_Min", 
                                      FilterValue= M3F$values["min"], 
                                      FilterDescription="DP M3 Het Ler", 
                                      FilterSign="less",  
                                      FilterNA=FALSE,
                                      ConditionalFilterVector=geno(GVCFL[[3]]@VCFlike@VCF)$GT, 
                                      ConditionalFilterValue="0/1",
                                      append=TRUE)

GVCFL[[4]]@VCFlike@VCF <- AnnotateVCF(GVCFL[[4]]@VCFlike@VCF, 
                                      vcfField="LER_DP", 
                                      vcfFieldType=c("info"), 
                                      FilterType=c("numeric"), 
                                      FilterName="DP_M4_Ler_Het_Min", 
                                      FilterValue= M4F$values["min"], 
                                      FilterDescription="DP M4 Het Ler", 
                                      FilterSign="less",  
                                      FilterNA=FALSE,
                                      ConditionalFilterVector=geno(GVCFL[[4]]@VCFlike@VCF)$GT, 
                                      ConditionalFilterValue="0/1",
                                      append=TRUE)

# Export Individuals in  VCF format
if(! file.exists(paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,sep="/"))){
    stopifnot(dir.create(paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,sep="/"), showWarnings = TRUE, recursive = FALSE, mode = "0770"))
}


if(! file.exists(paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,"VCF_indiv",sep="/"))){
    stopifnot(dir.create(paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,"VCF_indiv",sep="/"), showWarnings = TRUE, recursive = FALSE, mode = "0770"))
}

for(i in 1:4){
    writeVcf(GVCFL[[i]]@VCFlike@VCF, file= paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,"VCF_indiv",paste("M",i,".vcf",sep=""),sep="/"))
}

## 4) Instanciate the Tetrad

VCFTetrad2 <- new("Tetrad",GVCFL)  

## 5) Find NCO

VCFTetrad2 <- FindNCO(VCFTetrad2)

## 6) Find CO

InputFileWoExt <- paste( basename(TetradPath), "_" , basename(GoldVariantPath), "_", "genotype",sep="")
TetradNameWithCoParam <- paste( basename(TetradPath),"20k5k5k",sep="_")
InputDirPathCO <- paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,"CO","InputCrossOver",sep="/")
OutputDirPathCO <- paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,"CO","OutputCrossOver",sep="/")

if(! dir.exists(paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,"CO",sep="/" ))){
    stopifnot(dir.create(paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,"CO",sep="/" ), showWarnings = TRUE, recursive = FALSE, mode = "0770"))
}

if(! dir.exists(InputDirPathCO )){
    stopifnot(dir.create(InputDirPathCO , showWarnings = TRUE, recursive = FALSE, mode = "0770"))
}

if(! dir.exists(OutputDirPathCO )){
    stopifnot(dir.create(OutputDirPathCO , showWarnings = TRUE, recursive = FALSE, mode = "0770"))
}

VCFTetrad2 <- FindCO(VCFTetrad2,"InputFileWoExt"=InputFileWoExt,"TetradNameWithCoParam"=TetradNameWithCoParam, "InputDirPath"=InputDirPathCO, "OutputDirPath"=OutputDirPathCO)

## Export the results: Tetrad in vcf format, CO table, NCO table and report table.

if(! dir.exists(paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,"VCF_tetrad",sep="/"))){
    stopifnot(dir.create(paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,"VCF_tetrad",sep="/"), showWarnings = TRUE, recursive = FALSE, mode = "0770"))
}

if(! dir.exists(paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,"NCO",sep="/" ))){
    stopifnot(dir.create(paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,"NCO",sep="/" ), showWarnings = TRUE, recursive = FALSE, mode = "0770"))
}

if(! dir.exists(paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,"REPORT",sep="/"))){
    stopifnot(dir.create(paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,"REPORT",sep="/"), showWarnings = TRUE, recursive = FALSE, mode = "0770"))
}

writeTetrad(VCFTetrad2,
            VCFtetradFile=paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,"VCF_tetrad",paste(basename(TetradPath),".vcf",sep=""),sep="/"),
            CoFile = paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,"CO",paste(basename(TetradPath),"_",basename(GoldVariantPath),"_20k5k5k_ListOfCO.txt",sep=""),sep="/"),
            NcoFile = paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,"NCO",paste(basename(TetradPath),"_",basename(GoldVariantPath),"_ListOfNCO.txt",sep=""),sep="/"),
            ReportFile = paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,"REPORT",paste(basename(TetradPath),"_",basename(GoldVariantPath),"_report.txt",sep=""),sep="/"),
	    sdiFile=sdiFilePath)


## Draw the plot

if(! dir.exists(paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,"Plot",sep="/"))){
    stopifnot(dir.create(paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,"Plot",sep="/"), showWarnings = TRUE, recursive = FALSE, mode = "0770"))
}

pdf(paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,"Plot",paste( basename(TetradPath),basename(GoldVariantPath),"_Genotype.pdf",sep=""),sep="/"))
plot(VCFTetrad2)
dev.off()

pdf(paste(TetradPath,"Analysis",basename(GoldVariantPath),AnalysisName,"Plot",paste( basename(TetradPath),basename(GoldVariantPath),"_ProbaOfGenotype_Distribution.pdf",sep=""),sep="/"))
plotGL(VCFTetrad2)
dev.off()

