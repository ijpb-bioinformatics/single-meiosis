TetradPath <-
ListIndiv <- c("M1","M2","M3","M4")

VcfSuffix <- "sorted_reformated"
Chromosom <- paste("Chr",1:5,sep="")

GoldVariantPath <- "/projects/Single_Meiosis/Ressources/Ler_LUHQ01/SI_Schneeberger/Variant_All" 
sdiFilePath <- "/projects/Single_Meiosis/Ressources/Ler_LUHQ01/SI_Schneeberger/Variant_All/newsdi.txt"
ColPositionFilePathPref <- "/projects/Single_Meiosis/Ressources/Ler_LUHQ01/SI_Schneeberger/Variant_All/dfColPos"
LerPositionFilePathPref <- "/projects/Single_Meiosis/Ressources/Ler_LUHQ01/SI_Schneeberger/Variant_All/dfLerPos"



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
            write.table(resultat[order(resultat$Chrom,resultat$colPos),],file=paste(TetradPath,"/","Analysis","/",basename(GoldVariantPath),"/","SM","/",ListIndiv[i],"_","SM.txt",sep=""),quote=FALSE,s
ep="\t",row.names=FALSE,col.names=TRUE)
        }
        else{
            write.table(resultat[order(resultat$Chrom,resultat$colPos),],file=paste(TetradPath,"/","Analysis","/",basename(GoldVariantPath),"/","SM","/",ListIndiv[i],"_","SM.txt",sep=""),quote=FALSE,s
ep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
         }
        
    }
}

