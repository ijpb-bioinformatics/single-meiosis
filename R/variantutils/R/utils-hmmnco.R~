
# CONVERT SM VCFlike Format to input HMM and export .Rdata 

SM2InputHMM <- function(VCFL, OutFile){

	dfVCFL <- lapply( 1:4, function(x){
	# UNIQ_ID CHROMOSOM POSITION LER_POS COL_REF LER_REF LENGTH COL_DP LER_DP COL_QSUM LER_QSUM COL_MAPQSUM LER_MAPQSUM
	# ID Chrom Pos R Y Ind
	 	data.frame(
		  "ID"=VCFL[[x]]@VCFlike@Df$UNIQ_ID,
		  "Chrom"=VCFL[[x]]@VCFlike@Df$CHROMOSOM,
		  "Pos"=VCFL[[x]]@VCFlike@Df$POSITION,
		  "R"=VCFL[[x]]@VCFlike@Df$COL_DP+VCFL[[x]]@VCFlike@Df$LER_DP,
		  "Y"=VCFL[[x]]@VCFlike@Df$COL_DP,
		  "Ind"=rep(paste("M",x,sep=""),dim(VCFL[[x]]@VCFlike@Df)[1]))
	})
	names(dfVCFL)=paste("M",1:4,sep="")

# Remove non informativ SNP
	toberemoved <- unique(c(which(dfVCFL$M1$R == 0),which( dfVCFL$M2$R == 0),which( dfVCFL$M3$R == 0),which (dfVCFL$M4$R == 0)))  	

	if(length(toberemoved) != 0){
	dfVCFL<-lapply(dfVCFL, function(x){
		x[-toberemoved,]
	})
	}
	assign(basename(OutFile),dfVCFL)
	save(list=basename(OutFile),file=paste(OutFile,".Rdata",sep=""))
}
	
