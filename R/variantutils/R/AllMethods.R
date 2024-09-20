#############
## VCFlikeParser
############

#' Constructor method of the [\code{\link{VCFlikeParser}}] Class.
#' @param File A string giving the path to the VCF like file.
#' @param Format A string giving the format of the VCF like file which is also the name of the VCFlike subclass.
#' @param RefGenome A string giving the reference genome
#' @param Sample A string giving the name of the sample  
#' @return An object of class [\code{\link{VCFlikeParser}}]
#' @name VCFlikeParser
#' @rdname VCFlikeParser-class

setMethod("initialize",
          signature(.Object="VCFlikeParser"),
          function(.Object, File, Format, RefGenome,Sample){ 
            .Object@VCFlike = switch(Format,
                                     "sdi-9" = new("sdi-9",File=File,RefGenome=RefGenome,Sample=Sample),
                                     "SM" =  new("SM",File=File,RefGenome=RefGenome,Sample=Sample),
                                     "RefPol" = new("RefPol",File=File,RefGenome=RefGenome,Sample=Sample),
                                     "JGI-aa" = new("JGI-aa",File=File,RefGenome=RefGenome,Sample=Sample),
                                     "JGI-indel" = new("JGI-indel",File=File,RefGenome=RefGenome,Sample=Sample),
                                     "sdi-7" = new("sdi-7",File=File,RefGenome=RefGenome,Sample=Sample),
                                     "VAST" = new("VAST",File=File,RefGenome=RefGenome,Sample=Sample),
                                     "SALK" = new("SALK",File=File,RefGenome=RefGenome,Sample=Sample)
            )
            .Object@VCFlike@VCF = ConvertToVcf(.Object@VCFlike)
            .Object@Format=Format
            return(.Object)
          }
)


###########
### sdi-9
##########

#' Constructor method fof the [\code{\link{sdi-9-like}}] Class.
#' @param File A string giving the path to the VCF like file.
#' @param RefGenome A string giving the reference genome
#' @param Sample A string giving the name of the sample
#' @seealso VCFlike-class, VCFlikeParser-class
#' @name sdi-9-class
#' @rdname sdi-9-class

setMethod("initialize",
          signature(.Object="sdi-9"),
          function(.Object, File, RefGenome,Sample){
            .Object@RefGenome=RefGenome
            .Object@Sample=Sample
            .Object@File=File
            .Object@Df=read.table(File,sep="\t",fill=TRUE,na.strings=c("NA",""))
            
            # Add SVTYPE
            .Object@Df$SVTYPE <- ifelse(.Object@Df[,3]==0,"SUB",ifelse(.Object@Df[,3]<0,"DEL","INS"))
                       
            dfNames <- c(
              "CHROMOSOME",
              "POSITION",
              "LENGTH",
              "REFERENCE_BASE",
              "CONSENSUS_BASE",
              "IMR_DENOM_SCORE",
              "PHRED_SCORE",
              "HMQ", 
              "HMQ_CONS",
              "SVTYPE"
              )
            
            names(.Object@Df) <- dfNames
            
            ContigTmp <- as.factor(.Object@Df$CHROMOSOME)
            tmp <- TestChromFormat(RefGenome = RefGenomes[[RefGenome]] , ChromosomeNames = levels(ContigTmp)) 
          
            if( ! is.null(tmp)){
              if(any(is.na(tmp$UserFileChromName))){
                levels(ContigTmp)<- tmp$RefGenomChromName[na.omit(match(levels(ContigTmp),tmp$UserFileChromName))]
                
              }
              else{
                levels(ContigTmp)<- tmp$RefGenomChromName[match(levels(ContigTmp),tmp$UserFileChromName)]
              }
            }
    
            .Object@Df$CHROMOSOME <- as.vector(ContigTmp)
            
            .Object@MetaDf <- DataFrame(
              Number=rep(".",10),
              Type=c("String","Integer","Double","String","String","Integer","Integer","Double","String","String"),
              Description=c(
                "chromosome: The same name as the chromosome id in reference fasta file",
                "position: 1-based leftmost position",
                "length: the length difference of the changed sequence against reference (0 for SNPs, negative for deletions, positive for insertions)",
                "reference base: [-A-Z]+ (regular expression range)",
                "consensus base: [-A-Z]+ (regular expression range), IUPAC code is used for heterozygous sites",
                "detection score: (1 detected by both IMR and DENOM, 2 detected by IMR as heterozygous SNP but detected by DENOM as homozygous SNP, 4: detected by IMR only, 5: detected by DENOM only)",
                "Phred score: (for SNPs only)",
                "HMQ: high mapping quality score",
                "HMQ consensus base: the consensus base when considering reads with high mapping quality only",
                "SVTYPE: Type of structural variant"
              ))
            
            dimnames(.Object@MetaDf)[[1]] <- dfNames
          
            return(.Object)
          }
)

#' ConvertToVcf method for the sdi-9 class
#' @details 
#' This function is called by the class VCFlikeParser to set the slot VCF of an object of class sdi-9. 
#' @param VCFlike A VCFlike object
#' @return VCF
#' @seealso VCFlike-class, VCFlikeParser-class 
#' @name sdi-9-class
#' @rdname sdi-9-class

setMethod(
  f ="ConvertToVcf",
  signature ="sdi-9",
  definition = function( object ){
    
  
    infodf = DataFrame(
      "LENGTH"=as.character(object@Df$LENGTH),
      "IMR_DENOM_SCORE"=as.character(object@Df$IMR_DENOM_SCORE),
      "SVTYPE"=as.character(object@Df$SVTYPE)
      )
    
    # Convert IUPAC code
    ALT=as.character(object@Df$CONSENSUS_BASE)
    ALT[which(ALT %in% "-")] <- ""
    
    REF <- as.character(object@Df$REFERENCE_BASE)
    REF[which(REF %in% "-")] <- ""
    
    INDELindex <- which(object@Df$LENGTH != 0)
    INDELBASES = rep(NA,length(INDELindex))
    
    INDELBASES <- FetchBaseBeforEvent(object, INDELindex,"CHROMOSOME","POSITION")
    
    REF[INDELindex] <- apply(cbind(INDELBASES,REF[INDELindex]),1,paste,collapse="")
    ALT[INDELindex] <- apply(cbind(INDELBASES,unlist(ALT[INDELindex])),1,paste,collapse="")
    
    ALTlist=as.list(ALT)
    ALTlist[which(ALT %in%  names(IUPAC_CODE_MAP)[-c(1:4)])] <- lapply(IUPAC_CODE_MAP[unlist(ALT[which(ALT %in%  names(IUPAC_CODE_MAP)[-c(1:4)])])],s2c)
    ALTlistForID <-  unlist(lapply(ALTlist,function(x){
      paste(x,collapse=",")
    }))
    
    # CHANGE POSITION FOR INDEL, give POSITION-1
    POSITION <- object@Df$POSITION
    POSITION[INDELindex] <- POSITION[INDELindex] - 1
    
    fixeddf=DataFrame(
      "REF"=DNAStringSet(as.character(REF)),
      "ALT"=DNAStringSetList(ALTlist), 
      "QUAL"= as.numeric(object@Df$PHRED_SCORE),
      "FILTER"=as.character(rep("PASS",dim(object@Df)[1])))
    
    width=object@Df$LENGTH
    width[which(width<0)]<- c(-1)
    
    header=VCFHeader(reference= object@RefGenome, samples= object@Sample, header=DataFrameList("INFO"=object@MetaDf[c(10,3,6),]))
    
    VCF=VCF(
      rowRanges=GRanges(seqnames = Rle(object@Df$CHROMOSOME),IRanges(start=POSITION,width=as.numeric(width+1),names=paste(object@Df$CHROMOSOME,"_",POSITION,"_",REF,"_",ALTlistForID ,sep=""))),
      fixed=fixeddf,
      info=infodf,
      exptData=SimpleList(header=header),
      collapsed=TRUE)
    
    genome(seqinfo(rowRanges(VCF))) <- object@RefGenome
    
    
    return(VCF)
  }
)

##########
## SM
#########

#' Constructor method fof the [\code{\link{SM}}] Class.
#' @param File A string giving the path to the VCF like file.
#' @param RefGenome A string giving the reference genome
#' @param Sample A string giving the name of the sample
#' @name SM-class
#' @rdname SM-class

setMethod("initialize",
          signature(.Object="SM"),
          function(.Object, File, RefGenome,Sample){
            .Object@RefGenome=RefGenome
            .Object@Sample=Sample
            .Object@File=File
            .Object@Df=read.table(File,h=TRUE)[,-3]
            dfNames <- c(
              "UNIQ_ID",
              "CHROMOSOM",
              "POSITION",
              "LER_POS",
              "COL_REF",
              "LER_REF",
              "LENGTH",
              "COL_DP",
              "LER_DP",
              "COL_QSUM",
              "LER_QSUM",
              "COL_MAPQSUM",
              "LER_MAPQSUM"
            )
            names(.Object@Df) <- dfNames
            
            .Object@MetaDf <- DataFrame(
              Number=rep(".",13),
              Type=c("String","String","Integer","Integer","String","String","Integer","Integer","Integer","Integer","Integer","Integer","Integer"),
              Description=c(
                "UNIQ_ID: Unique id created by pasting the COL_CHROM and COL_POS",
                "CHROMOSOM: Chromosom localisation",
                "POSITION: Position in Col",
                "LER_POS: Syntenic Position in Ler, (-1) for deleted position in Ler",
                "COL_REF: Reference base: [-A-Z]+ (regular expression range)",
                "LER_REF: Reference base: [-A-Z]+ (regular expression range)",
                "LENGTH: The length difference of the changed sequence against reference (0 for SNPs, positive for insertions and deletions)",
                "COL_DP: Depth on Col genome",
                "LER_DP: Depth on Ler genome",
                "COL_QSUM: Sum of base's quality value of the reads mapping on Col genome",
                "LER_QSUM: Sum of base's quality value of the reads mapping on Ler genome",
                "COL_MAPQSUM: Sum of read's mapping quality value of the reads mapping on Col genome",
                "LER_MAPQSUM: Sum of read's mapping quality value of the reads mapping on Ler genome"
              ))
            
            dimnames(.Object@MetaDf)[[1]] <- dfNames
           
            
            warning(paste(length(which(duplicated(.Object@Df))==TRUE),"duplicated ID have been removed"))
            
            # Attention aux lignes dupliquées
            .Object@Df <- .Object@Df[! duplicated(.Object@Df),]
            
            #.Object@VCF <- NULL
            
            return(.Object)
          }  
)

#' @title ConvertToVcf Method for the class  [\code{\link{SM}}]
#' @details This function is called by the class VCFlikeParser to set the slot VCF of an object of class  [\code{\link{SM}}]. 
#' @param VCFlike A VCFlike object
#' @return VCF
#' @name SM-class
#' @rdname SM-class

setMethod(
  f ="ConvertToVcf",
  signature ="SM",
  definition = function( object ){
    
    infodf = DataFrame(
      "LER_POS"=as.character(object@Df$LER_POS),
      "LENGTH"=as.character(object@Df$LENGTH),
      "COL_DP"=as.integer(object@Df$COL_DP),
      "LER_DP"=as.integer(object@Df$LER_DP),
      "COL_QSUM"=as.integer(object@Df$COL_QSUM),
      "LER_QSUM"=as.integer(object@Df$LER_QSUM),
      "COL_MAPQSUM"=as.integer(object@Df$COL_MAPQSUM),
      "LER_MAPQSUM"=as.integer(object@Df$LER_MAPQSUM),
      "UNIQ_ID"=as.character(object@Df$UNIQ_ID)
    )
    
    
    # Convert IUPAC code
    ALT=as.character(object@Df$LER_REF)
    ALT[which(ALT %in% "-")] <- ""
    
    REF <- as.character(object@Df$COL_REF)
    REF[which(REF %in% "-")] <- ""
    
    INDELindex <- which(object@Df$LENGTH != 0)
    INDELBASES = rep(NA,length(INDELindex))
    
    INDELBASES <- FetchBaseBeforEvent(object, INDELindex,"CHROMOSOM","POSITION")
    
    REF[INDELindex] <- apply(cbind(INDELBASES,REF[INDELindex]),1,paste,collapse="")
    ALT[INDELindex] <- apply(cbind(INDELBASES,unlist(ALT[INDELindex])),1,paste,collapse="")
    
    ALTlist=as.list(ALT)
    ALTlist[which(ALT %in%  names(IUPAC_CODE_MAP)[-c(1:4)])] <- lapply(IUPAC_CODE_MAP[unlist(ALT[which(ALT %in%  names(IUPAC_CODE_MAP)[-c(1:4)])])],s2c)
    ALTlistForID <-  unlist(lapply(ALTlist,function(x){
      paste(x,collapse=",")
    }))
    
    
    fixeddf=DataFrame(
      "REF"=DNAStringSet(as.character(REF)),
      "ALT"=DNAStringSetList(ALTlist), 
      "QUAL"= as.numeric(rep(NA,dim(object@Df)[1])),
      "FILTER"=as.character(rep("PASS",dim(object@Df)[1])))
    
    width=object@Df$LENGTH
    width[which(width<0)]<- c(-1)
  

      
    header=VCFHeader(reference= object@RefGenome, samples= object@Sample, header=DataFrameList("INFO"=object@MetaDf[c(4,7,8:13,1),]))
    
    VCF=VCF(
      rowRanges=GRanges(seqnames = Rle(object@Df$CHROMOSOM),IRanges(start=object@Df$POSITION,width=as.numeric(width+1),names=paste(object@Df$CHROMOSOM,"_",object@Df$POSITION,"_",REF,"_",ALTlistForID ,sep=""))),
      fixed=fixeddf,
      info=infodf,
      exptData=SimpleList(header=header),
      #change 12 mai
      #geno= SimpleList(lapply(vcf$GENO, `dimnames<-`, NULL)),
      colData=DataFrame(Samples=seq_along(object@Sample),row.names=object@Sample),
      collapsed=TRUE)
    
    genome(seqinfo(rowRanges(VCF))) <- object@RefGenome
    
    return(VCF)
  }
)


##########
## RefPol
#########

#' Constructor method fof the [\code{\link{NOEL}}] Class.
#' @param File A string giving the path to the VCF like file.
#' @param RefGenome A string giving the reference genome
#' @param Sample A string giving the name of the sample
#' @name RefPol-class
#' @rdname RefPol-class

setMethod("initialize",
          signature(.Object="RefPol"),
          function(.Object, File, RefGenome,Sample){
            .Object@RefGenome=RefGenome
            .Object@Sample=Sample
            .Object@File=File
            .Object@Df=read.table(File,h=TRUE)[,-3]
            
            dfNames <- c(
              "UNIQ_ID",
              "CHROMOSOM",
              "POSITION",
              "LER_POS",
              "LENGTH",
              "COL_REF",
              "LER_REF",
              "IMR_DENOM_SCORE",
              "PHRED_SCORE",
              "HMQ", 
              "IS_IN_TE",
              "IS_HETEROZYGOUS",
              "IS_IN_CNV",
              "COL_DP_F1",
              "LER_DP_F1",
              "COL_QSUM_F1",
              "LER_QSUM_F1",
              "COL_MAPQSUM_F1",
              "LER_MAPQSUM_F1",
              "COL_DP_Parent",
              "LER_DP_Parent",
              "COL_QSUM_Parent",
              "LER_QSUM_Parent",
              "COL_MAPQSUM_Parent",
              "LER_MAPQSUM_Parent"
            )
            names(.Object@Df) <- dfNames
            
            .Object@MetaDf <- DataFrame(
              Number=rep(".",25),
              Type=c("String","String","Integer","Integer","Integer","String","String","Integer","Integer","Integer","String","String","String","Integer","Integer","Integer","Integer","Integer","Integer","Integer","Integer","Integer","Integer","Integer","Integer"),
              Description=c(
                "Unique id created by pasting the CHROMOSOM and POSITION",
                "Chromosom",
                "Position in Col",
                "Syntenic Position in Ler, (-1) for deleted position in Ler",
                "The indel size",
                "COL_REF, Reference base: [-A-Z]+ (regular expression range)",
                "LER_REF, Consensus base: [-A-Z]+ (regular expression range)",
                "detection score: (1 detected by both IMR and DENOM, 2 detected by IMR as heterozygous SNP but detected by DENOM as homozygous SNP, 4: detected by IMR only, 5: detected by DENOM only)",
                "Phred score: (for SNPs only)",
                "HMQ: high mapping quality score", 
                "IS_IN_TE, boolean, TRUE if the polymorphism is in a TE (TAIR10, GFF)",
                "IS_HETEROZYGOUS, boolean, TRUE if the polymorphism is in a TE (TAIR10, GFF)",
                "IS_IN_CNV, boolean, TRUE if the polymorphism is in a CNV (Ler-0, shore)",
                "COL_DP_F1, The sequencing depth on Col genome in the F1  at this position",
                "LER_DP_F1, The sequencing depth on Ler genome in the F1 at this position",
                "COL_QSUM_F1, The sum of the base quality on Col genome in the F1 at this position",
                "LER_QSUM_F1, The sequencing depth on Ler genome in the F1 at this position",
                "COL_MAPQSUM_F1, The sum of the mapping quality of the reads covering this position on Col genome in the F1",
                "LER_MAPQSUM_F1, The sum of the mapping quality of the reads covering this position on Ler genome in the F1",
                "COL_DP_Parent, The sequencing depth on Col genome in the Col Parent at this position ",
                "LER_DP_Parent, The sequencing depth on Col genome in the Col Parent at this position ",
                "COL_QSUM_Parent, The sum of the base quality on Col genome in the Col Parent at this position",
                "LER_QSUM_Parent, The sum of the base quality on Ler genome in the Ler Parent at this position",
                "COL_MAPQSUM_Parent, The sum of the mapping quality of the reads covering this position on Col genome in the Col Parent",
                "LER_MAPQSUM_Parent, The sum of the mapping quality of the reads covering this position on Ler genome in the Ler Parent"
              ))
            
            dimnames(.Object@MetaDf)[[1]] <- dfNames
            
            return(.Object)
          }
)     
            

#' @title ConvertToVcf Methof for the class  [\code{\link{RefPol}}]
#' @details This function is called by the class VCFlikeParser to set the slot VCF of an object of class  [\code{\link{RefPol}}]. 
#' @param VCFlike A VCFlike object
#' @return VCF
#' @name RefPol-class
#' @rdname RefPol-class

setMethod(
  f ="ConvertToVcf",
  signature ="RefPol",
  definition = function( object ){
    
    infodf = DataFrame(
      "UNIQ_ID"=as.character(object@Df$UNIQ_ID),
      "LER_POS"=as.character(object@Df$LER_POS),
      "LENGTH"=as.character(object@Df$LENGTH),
      "IMR_DENOM_SCORE"=as.character(object@Df$IMR_DENOM_SCORE),
      "HMQ"=as.integer(object@Df$HMQ), 
      "IS_IN_TE"=as.character(object@Df$IS_IN_TE),
      "IS_HETEROZYGOUS"=as.character(object@Df$IS_HETEROZYGOUS),
      "IS_IN_CNV"=as.character(object@Df$IS_IN_CNV),
      "COL_DP_F1"=as.integer(object@Df$COL_DP_F1),
      "LER_DP_F1"=as.integer(object@Df$LER_DP_F1),
      "COL_QSUM_F1"=as.integer(object@Df$COL_QSUM_F1),
      "LER_QSUM_F1"=as.integer(object@Df$LER_QSUM_F1),
      "COL_MAPQSUM_F1"=as.integer(object@Df$COL_MAPQSUM_F1),
      "LER_MAPQSUM_F1"=as.integer(object@Df$LER_MAPQSUM_F1),
      "COL_DP_Parent"=as.integer(object@Df$COL_DP_Parent),
      "LER_DP_Parent"=as.integer(object@Df$LER_DP_Parent),
      "COL_QSUM_Parent"=as.integer(object@Df$COL_QSUM_Parent),
      "LER_QSUM_Parent"=as.integer(object@Df$LER_QSUM_Parent),
      "COL_MAPQSUM_Parent"=as.integer(object@Df$COL_MAPQSUM_Parent),
      "LER_MAPQSUM_Parent"=as.integer(object@Df$LER_MAPQSUM_Parent)
    )
    
    
    # Convert IUPAC code
    ALT=as.character(object@Df$LER_REF)
    ALT[which(ALT %in% "-")] <- ""
    
    REF <- as.character(object@Df$COL_REF)
    REF[which(REF %in% "-")] <- ""
    
    INDELindex <- which(object@Df$LENGTH != 0)
    INDELBASES = rep(NA,length(INDELindex))
    
    INDELBASES <- FetchBaseBeforEvent(object, INDELindex,"CHROMOSOM","POSITION")
    
    REF[INDELindex] <- apply(cbind(INDELBASES,REF[INDELindex]),1,paste,collapse="")
    ALT[INDELindex] <- apply(cbind(INDELBASES,unlist(ALT[INDELindex])),1,paste,collapse="")
    
    ALTlist=as.list(ALT)
    ALTlist[which(ALT %in%  names(IUPAC_CODE_MAP)[-c(1:4)])] <- lapply(IUPAC_CODE_MAP[unlist(ALT[which(ALT %in%  names(IUPAC_CODE_MAP)[-c(1:4)])])],s2c)
    ALTlistForID <-  unlist(lapply(ALTlist,function(x){
      paste(x,collapse=",")
    }))
    
    
    fixeddf=DataFrame(
      "REF"=DNAStringSet(as.character(REF)),
      "ALT"=DNAStringSetList(ALTlist), 
      "QUAL"= as.integer(object@Df$PHRED_SCORE),
      "FILTER"=as.character(rep("PASS",dim(object@Df)[1])))
    
    width=object@Df$LENGTH
    width[which(width<0)]<- c(-1)
    
    
    header=VCFHeader(reference= object@RefGenome, samples= object@Sample, header=DataFrameList("INFO"=object@MetaDf[c(1,4,5,8,10:25),]))
    
    VCF=VCF(
      rowRanges=GRanges(seqnames = Rle(object@Df$CHROMOSOM),IRanges(start=object@Df$POSITION,width=as.numeric(width+1),names=paste(object@Df$UNIQ_ID,"_",REF,"_",ALTlistForID ,sep=""))),
      fixed=fixeddf,
      info=infodf,
      exptData=SimpleList(header=header),
      collapsed=TRUE)
    
    genome(seqinfo(rowRanges(VCF))) <- object@RefGenome
    
    return(VCF)
  }
)

            

#' @title Genotype Method for the class [\code{\link{SM}}]
#' @details We followed the VCF specification to add genotype information. A call to this function
#' add the required Meta-information relativ to the genotype to the header of the VCF object. It also
#' add the slot [\code{geno}] to the VCF thanks to the [\code{\link{geno}}] accessors.
#' At each variant site, we computed the three probabilities to be: 
#' \itemize{
#' \item (i) heterozygous Col-/Ler by considering that the proportion of either allele should follow a binomial distribution with a probability of success of 0.5
#' \item (ii) homozygous Col by considering that the proportion of the reference allele should follow a binomial distribution with a probability of success for 
#' the reference allele of 1 minus the error sequencing rate 
#' \item (iii) homozygous Ler (idem to ii).
#' }
#' The Genotype at variant site is given by the highest probability. 
#' The genotype encoding is: 
#' \itemize{ 
#' \item 0/0 for the homozygous Col genotype
#' \item 0/1 for the heterozygous Col/Ler genotype.
#' \item 1/1 for the homozygous Ler genotype
#' }
#' @param object An object of class SM
#' @param SequencingErrorRate The estimated sequencing error rate 
#' @return An object of class [\code{\link{SM}}]
#' @aliases Genotype, SM-method
#' @rdname SM-class


setMethod(
  f ="Genotype",
  signature ="SM",
  definition = function( object , SequencingErrorRate){

      #Change 12 mai
    #colData(object@VCF)<-DataFrame(Samples=as.integer(1),row.names=object@Sample)
    
    colDP <- info(object@VCF)$COL_DP
    lerDP <- info(object@VCF)$LER_DP

 
    
    hypcol = mapply(function(x,y){
      if(!is.na(x) & !is.na(y)){ 
        tmp <- dbinom(y,x+y,SequencingErrorRate)  
      }
      else{
        NA
      }
    }, colDP, lerDP)
    
    
    hyphet = mapply(function(x,y){
      if(!is.na(x) & !is.na(y)){ 
         dbinom(y,x+y,1/2)   
      }
      else{
        NA
      }
    }, colDP, lerDP)
    
    hypler =  mapply(function(x,y){
      if(!is.na(x) & !is.na(y)){ 
         dbinom(x,x+y,SequencingErrorRate)    
      }
      else{
        NA
      }
    }, colDP, lerDP)
    
    Genotype=rep("./.",length(colDP))
    
    Genotype[which(hypcol > hypler & hypcol > hyphet)] <- "0/0"
    Genotype[which(hypler > hyphet & hypler > hypcol)] <- "1/1"
    Genotype[which(hyphet > hypcol & hyphet > hypler)] <- "0/1"
    #Genotype[which(Genotype=="0/0" & (lerDP >= 1))] <- "./."
    #Genotype[which((log10(hypcol)-log10(hypler)) >=3 & (log10(hypcol)-log10(hyphet))>=3)] <- "0/0"
    #Genotype[which((log10(hypler)-log10(hyphet)) >=3 & (log10(hypler)-log10(hypcol))>=3)] <- "1/1"
    #Genotype[which((log10(hyphet)-log10(hypcol)) >=3 & (log10(hyphet)-log10(hypler))>=3)] <- "0/1" 
    
    GenotypeLik=paste(round(hypcol,6),round(hyphet,12),round(hypler,12),sep=",")
    
    GenotypeAsMatrix=as.matrix(Genotype)
    #dimnames(GenotypeAsMatrix) <- list(info(object@VCF)$UNIQ_ID,object@Sample)
    dimnames(GenotypeAsMatrix) <- dimnames(object@VCF)
    
    GenotypeLikAsMatrix=as.matrix(GenotypeLik)
    #dimnames(GenotypeLikAsMatrix) <- list(info(object@VCF)$UNIQ_ID,object@Sample)
      dimnames(GenotypeLikAsMatrix) <- dimnames(object@VCF)
    
    geno(object@VCF) <- SimpleList(GT=GenotypeAsMatrix,GL=GenotypeLikAsMatrix)
      
    genodf <- DataFrame(
      Number=rep("1",2),
      Type=c("String","String"),
      Description=c("Genotype", "Genotype Likelihoods")
      )
    row.names(genodf)=c("GT","GL")
    
    geno(header(object@VCF)) <- genodf
    #colData(object@VCF)<-DataFrame(Samples=as.integer(1),row.names=object@Sample)
    return(object)
  }
)

#' @title Filter Method for the class [\code{\link{SM}}]
#' @details We followed the VCF specification to add filter information. A call to this function
#' add the required Meta-information relativ to the Filter field to the header of the VCF object. It also
#' adds the slot [\code{filter}] to the VCF thanks to the [\code{\link{filer}}] accessors.
#' The encoding of the filter in the filter field (of the fixed slot in the VCF) will be: FilterName_FilterSign_FilterValue.
#' @param Object an object of class SM
#' @param FilterName The name of the field that has been filter.  
#' @param FilterDescription The description of the filter in the VCF header
#' @param FilterValue The threshold value
#' @param FilterSign The sens of the filtering: more or less   
#' @param append Logical: TRUE if the filter has to replace the existing one, FALSE if it has to be append to the existing one.
#' @return An object of class [\code{\linkS4class{SM}}].
#' @exportMethod Filter
#' @aliases Filter, SM-method
#' @rdname SM-class

setMethod(
  f ="Filter",
  signature ="SM",
  # Code for filter: FieldSignValue (ex:COL_DP_more_5)
  # Description: "Normal low; (2.5%,5%); bilateral"
  
  # fixed(header(vcf)) # DataFrameList
  # fixed(header(vcf))$FILTER # DataFrame
  
  definition <- function(object, FilterName=c("COL_DP","LER_DP","COL_MAPQSUM","LER_MAPQSum","COL_QSUM","LER_QSUM"), 
                         FilterDescription=character(), FilterValue, FilterSign=c("more","less"), append=logical()){
    
    # check arguments
    if(! FilterName %in%  c("COL_DP","LER_DP","COL_MAPQSUM","LER_MAPQSum","COL_QSUM","LER_QSUM"))
      stop("Invalid FilterName")
    
    if(! (is.character(FilterDescription) | nchar(FilterDescription !=0)))
      stop("Invalid FilterDescription")
    
    if(! (is.numeric(FilterValue) & length(FilterValue)==1))
      stop("Invalid FilterValue")
    
    if(! (FilterSign %in%  c("more","less")))
      stop("Invalid FilterSign")
    
    if(! is.logical(append))
      stop("Invalid append")
    
    ValuesToFilter <- switch(FilterName,
                             "COL_DP"=info(object@VCF)$COL_DP,
                             "LER_DP"=info(object@VCF)$LER_DP,
                             "COL_MAPQSUM"=info(object@VCF)$COL_MAPQSUM,
                             "LER_MAPQSUM"=info(object@VCF)$LER_MAPQSUM,
                             "COL_QSUM"=info(object@VCF)$COL_DP_QSUM,
                             "LER_QSUM"=info(object@VCF)$LER_DP_QSUM
    )  
    
    if(FilterSign == "more"){
      ValuesToFilterStatus <-  ValuesToFilter >=  FilterValue
    }
    else if(ValuesToFilter == "less"){
      ValuesToFilterStatus <- ValuesToFilter <=  FilterValue
    }
    # get the DataFrame corresponding to filters in the VCF header 
    
    #FilterDF <- fixed(header(object@VCF))$FILTER
    FilterDF <- attr(fixed(header(object@VCF)),"header")$FILTER
    
    if(append == TRUE & ! is.null(FilterDF)){  
      # check if this filter already exists
      check <- grep(FilterName, row.names(FilterDF)) 
      if(length(check) != 0){
        # update this filter description 
        UpdatedFilterDF <- FilterDF
        UpdatedFilterDF[check,] <- paste(FilterName, FilterSign, "than", FilterValue, sep=" ")
        row.names(UpdatedFilterDF[check]) <- paste(FilterName, FilterSign, FilterValue, sep="_")
        # update the field filter
        FilterField <- ifelse(filt(object@VCF)=="PASS" & ValuesToFilterStatus,
                              "PASS",
                              paste(FilterName, FilterSign, FilterValue, sep="_"))
        
        if(filt(object@VCF) != "PASS" & ValuesToFilterStatus==FALSE){
          FilterField <- gsub(paste(paste(FilterName, FilterSign, "[0-9]*", sep="_"),paste(FilterName, FilterSign, FilterValue, sep="_"),filt(object@VCF)))
        }
        else   if(filt(object@VCF) != "PASS" & ValuesToFilterStatus==TRUE){
          tmp <- gsub(paste(paste(FilterName, FilterSign, "[0-9]*", sep="_"),"",filt(object@VCF)))
          FilterField <- ifelse(tmp==";","PASS","")
        }         
      }
      else{
        # append this filter description to existing one
        UpdatedFilterDF <- rbind(FilterDF,DataFrame("Description"= paste(FilterDescription, FilterSign, "than", FilterValue, sep=" ")))
        
        row.names(UpdatedFilterDF) <- c(row.names(FilterDF),paste(FilterName, FilterSign, FilterValue, sep="_"))
        
        # append this filter to existing one                                 
        # if PASS and TRUE => PASS, else filter
        FilterField <-  ifelse(filt(object@VCF) =="PASS" & ValuesToFilterStatus,
                               "PASS",
                               paste(FilterName, FilterSign, FilterValue, sep="_")) 
        # if != PASS and TRUE => nothing, else merge
        FilterField <- ifelse(filt(object@VCF) != "PASS" & ValuesToFilterStatus, 
                             filt(object@VCF) , 
                              paste(filt(object@VCF),paste(FilterName, FilterSign, FilterValue, sep="_"),sep=";"))                           
      }
    }
    if(append == FALSE | is.null(FilterDF)){ 
      # create a DataFrame object for filter
      UpdatedFilterDF <- DataFrame("Description"= paste(FilterDescription, FilterSign, "than", FilterValue, sep=" "))
      row.names(UpdatedFilterDF)= paste(FilterName, FilterSign, FilterValue, sep="_") 
      
      # Set the filter field
      FilterField <- ifelse(ValuesToFilterStatus==TRUE,
                            "PASS",
                            paste(FilterName, FilterSign, FilterValue, sep="_"))
    }
    #print(UpdatedFilterDF)
    #print(FilterField)
    
    #
    # FIXED replacement method for VCFHeader class is not yet avalaibale
    #
    
    #fixed(header(object@VCF))[["FILTER"]] <- UpdatedFilterDF
    filt(object@VCF) <- FilterField  
    
    return(object)
  }
)

#

setMethod(
  f ="AnnotateVCF",
  signature ="CollapsedVCF",
  definition <- function(object, vcfField, vcfFieldType=c("fixed","info","geno"), FilterType=c("character","numeric"), 
                   FilterName, FilterValue, FilterDescription, 
                   FilterSign=c("more","less"), FilterNA=FALSE, 
                   ConditionalFilterVector=NA,  ConditionalFilterValue=NA , 
                   append=TRUE){
 

    # Terms below or different
  
    if(FilterType=="character"){
      term <- "different from"
  }else if(FilterSign == "less"){     
    term <- "below"
  }
  else{
    term <- "up to"
  }
  
  # EXCEPTION for PV4
  
  if(vcfField == "PV4"){
    ValuesToFilterStatus <- unlist(lapply(info(object)[[vcfField]] <= FilterValue,function(x){length(which(x==TRUE))}))>0
  }
  
  else{
    
    # Get the value to filter
    
    if(vcfFieldType == "fixed"){
      ValuesToFilter <- fixed(object)[[vcfField]]
    }
    else if(vcfFieldType == "info"){
      ValuesToFilter <- unlist(info(object)[[vcfField]])
    }
    else{
      ValuesToFilter <- unlist(geno(object)[[vcfField]])
    }
    
    # print( ValuesToFilter)
    
    # Mark the value to filter by TRUE
    
    if( FilterType == "character"){    
      ValuesToFilterStatus <- ValuesToFilter == FilterValue    
    }
    
    else
      if( FilterType == "numeric"){
        if(FilterSign == "less"){
        ValuesToFilterStatus <- ValuesToFilter <= FilterValue      
        }
        else if(FilterSign == "more"){
          ValuesToFilterStatus <- ValuesToFilter >= FilterValue   
        }
      }
  }
  
  if(FilterNA==TRUE){
  ValuesToFilterStatus[which(is.na(ValuesToFilter))] <- TRUE
  }
  
  # print(ValuesToFilterStatus)
  
  # Change the ValuesToFilterStatus to FALSE 
  # if a ConditionalFilter does exist
  
  if(length(ConditionalFilterVector)>1){
    ValuesToFilterStatus[which(ConditionalFilterVector != ConditionalFilterValue)] <- FALSE
  }
  
  # get the filter field
  FilterField <- rowRanges(object)$FILTER
  
  # get the DataFrame corresponding to filters in the VCF header

  FilterDF <- attr(fixed(header(object)),"header")$FILTER
  
  # print(FilterField)
  # print(FilterDF)
  
  # Case: UPDATE FILTER, 
  
  if(append == TRUE & ! is.null(FilterDF)){
    
    # check if this filter already exists
    check <- grep(FilterName, row.names(FilterDF)) 
    if(length(check) != 0){
      
      FilterFieldTMP <- FilterField
      
      # update this filter description 
      UpdatedFilterDF <- FilterDF
      UpdatedFilterDF[check,] <- paste(FilterDescription, term , FilterValue, sep=" ")
      row.names(UpdatedFilterDF[check]) <- paste(FilterName, FilterValue, sep="")
      
      # update the field filter
      
      FilterFieldTMP[which(FilterField == "PASS" & ValuesToFilterStatus == TRUE)] <-
        paste(FilterName, FilterValue, sep="")
      
      if(FilterType == "numeric"){
        
        FilterFieldTMP[which(FilterField != "PASS" & ValuesToFilterStatus == TRUE)] <-
          gsub(paste(paste(FilterName, "[0-9]+", sep=""),paste(FilterName, FilterValue, sep="")),FilterField[which(FilterField != "PASS" & ValuesToFilterStatus == TRUE)])
      }
      else  FilterFieldTMP[which(FilterField != "PASS" & ValuesToFilterStatus == TRUE)] <-
        gsub(paste(paste(FilterName, "[a-z]+", sep=""),paste(FilterName, FilterValue, sep="")),FilterField[which(FilterField != "PASS" & ValuesToFilterStatus == TRUE)])
      
      if(FilterType == "numeric"){                          
        tmp <- gsub(paste(paste(FilterName,"[0-9]+", sep=""),"",FilterField[which(FilterField != "PASS" & ValuesToFilterStatus == FALSE)]))
      }
      else{
        tmp <- gsub(paste(paste(FilterName,"[a-z]+", sep=""),"",FilterField[which(FilterField != "PASS" & ValuesToFilterStatus == FALSE)]))
      }
      tmp <- ifelse(tmp==";","PASS","")
      FilterFieldTMP[which(FilterField != "PASS" & ValuesToFilterStatus == FALSE)] <- tmp
      
      FilterField <- FilterFieldTMP  
    }
    
    
    else{
      
      FilterFieldTMP <- FilterField
      
      # append this filter description to existing one
      UpdatedFilterDF <- rbind(FilterDF,DataFrame("Description"= paste(FilterDescription, term , FilterValue, sep=" ")))
      row.names(UpdatedFilterDF) <- c(row.names(FilterDF),paste(FilterName, FilterValue, sep=""))
      
      
      # append this filter description to existing one    
      FilterFieldTMP[which(FilterField == "PASS" & ValuesToFilterStatus==TRUE)] <-
        paste(FilterName, FilterValue, sep="")
      FilterFieldTMP[which(FilterField != "PASS" & ValuesToFilterStatus==TRUE)] <-
        paste(FilterField[which(FilterField != "PASS" & ValuesToFilterStatus == TRUE)],paste(FilterName,FilterValue, sep=""),sep=";")
      
      FilterField <- FilterFieldTMP                 
    }
    
  }

   
  
  if(append == FALSE | is.null(FilterDF)){ 
    # create a DataFrame object for filter
    UpdatedFilterDF <- DataFrame("Description"= paste(FilterDescription, term , FilterValue, sep=" "))
    row.names(UpdatedFilterDF)= paste(FilterName, FilterValue, sep="") 
    # print(UpdatedFilterDF)
    # Set the filter field
    FilterField[which(ValuesToFilterStatus==TRUE)] <- paste(FilterName, FilterValue, sep="")
    # print(FilterField)
  }     
  
  #print(FilterField )
  #print(UpdatedFilterDF)
   
   fixed(header(object))$FILTER <- UpdatedFilterDF
   filt(object) <- FilterField 
  
  return(object)
}
)


############
## Tetrad
############

setMethod("initialize",
          signature(.Object="Tetrad"),
          function(.Object, ListOfVCFLP){
            
            # Validity Check 
            
            if(length(ListOfVCFLP) != 4 ){
              stop("Invalid number of Sample, must be equal to 4")
            }
            checkGenotype <- lapply(ListOfVCFLP, function(x){
              is.null(geno(x@VCFlike@VCF)$GT)
            })
            if(any(unlist(checkGenotype))){
              stop("At least one sample had not yet been genotyped")
              # TODO: add sample name 
            }
            checkSMclass <- lapply(ListOfVCFLP, function(x){
              class(x@VCFlike) == "SM"
            })
            if( ! any(unlist(checkSMclass))){
              stop("At least one VCFlike is not of class SM")
              # TODO: add sample name 
            }
            
            # Extract VCF
            
            ListOfVCF <- lapply(ListOfVCFLP, function(x){
              x@VCFlike@VCF 
            })
            
            
            # Intesect Granges
            
            ListOfGRanges <- lapply(ListOfVCF,function(x){
            rowRanges(x)[which(filt(x) == "PASS")]
             })
           
            intersect <- names(ListOfGRanges[[1]]) %in% names(ListOfGRanges[[2]]) & 
              names(ListOfGRanges[[1]]) %in% names(ListOfGRanges[[3]]) & 
              names(ListOfGRanges[[1]]) %in% names(ListOfGRanges[[4]])
              
            fixedTetrad <- ListOfGRanges[[1]][intersect]
            
            #fixedTetrad <- Reduce(subsetByOverlaps,ListOfGRanges )
            
            # Intersect info
            
            ListOfInfo <- lapply(ListOfVCF,function(x){
              info(x)[which(filt(x) == "PASS"),]
            })
            
            ListOfSamples <- lapply(ListOfVCF,function(x){
              dimnames(colData(x))[[1]]
            })
            
            TetradInfo <- IntersectDataFrame(ListOfInfo,c("COL_DP", "LER_DP",  "COL_QSUM",  "LER_QSUM", "COL_MAPQSUM", "LER_MAPQSUM"),unlist(ListOfSamples))
            TetradInfo$NS <- as.integer(rep(4,dim(TetradInfo)[1]))
          
            # Intersect geno
            
            ListOfGT <- lapply(ListOfVCF,function(x){
              geno(x)$GT[which(filt(x) == "PASS"),]
            })
            
            ListOfGL <- lapply(ListOfVCF,function(x){
              geno(x)$GL[which(filt(x) == "PASS"),]
            })
            
            TetradGT <- as.matrix(IntersectGeno(ListOfGT))
            dimnames(TetradGT)[[2]] <- unlist(ListOfSamples)
            TetradGL <- as.matrix(IntersectGeno(ListOfGL))
            dimnames(TetradGL)[[2]] <- unlist(ListOfSamples)         
            
            TetradGeno <- SimpleList(GT=TetradGT,GL=TetradGL)
            
            # Construct TetradHeader 
            
            TMPH <-  as.data.frame(t(as.data.frame(lapply(info(header(ListOfVCF[[1]]))[3:8,]$Description,function(x){
              unlist(strsplit(x,":",perl=TRUE))
            }))))
            
            TetradDescription=c(paste(TMPH$V1,".",ListOfSamples[[1]],":",TMPH$V2," in sample ",ListOfSamples[[1]],sep=""),
                                paste(TMPH$V1,".",ListOfSamples[[2]],":",TMPH$V2," in sample ",ListOfSamples[[2]],sep=""),
                                paste(TMPH$V1,".",ListOfSamples[[3]],":",TMPH$V2," in sample ",ListOfSamples[[3]],sep=""),
                                paste(TMPH$V1,".",ListOfSamples[[4]],":",TMPH$V2," in sample ",ListOfSamples[[4]],sep=""))
            TetradHeaderInfo <- DataFrame(
              "Number" <- c(rep(info(header(ListOfVCF[[1]]))$Number[3:8],4),"."),
              "Type" <- c(rep(info(header(ListOfVCF[[1]]))$Type[3:8],4),"Integer"),
              "Description" <- c(TetradDescription,"Number of Samples")
            )
            names(TetradHeaderInfo)=c("Number","Type","Description")
            TMPH$V1 <- as.vector(TMPH$V1)
            dimnames(TetradHeaderInfo)[[1]]=c(paste(c(TMPH$V1),".",ListOfSamples[[1]],sep=""),
                                              paste(c(TMPH$V1),".",ListOfSamples[[2]],sep=""),
                                              paste(c(TMPH$V1),".",ListOfSamples[[3]],sep=""),
                                              paste(c(TMPH$V1),".",ListOfSamples[[4]],sep=""),
                                              "NS")
            TetradHeader=VCFHeader(reference= header(ListOfVCF[[1]])@reference, samples=unlist(ListOfSamples)) 
            info(TetradHeader)=TetradHeaderInfo
            geno(TetradHeader)=geno(header(ListOfVCF[[1]]))
            
            # Construct Tetrad VCF
            
            TetradVCF=VCF(
              #rowRanges=GRanges(seqnames = Rle(TetradChrom),TetradIRanges),
              #fixed=TetradFixed,
              rowRanges=fixedTetrad,
              fixed=fixedTetrad@elementMetadata,
              colData=DataFrame(Samples=rep(1,4),row.names=unlist(ListOfSamples)),
              info=TetradInfo,
              geno=TetradGeno,
              exptData=SimpleList(header=TetradHeader),
              collapsed=TRUE)
            
            .Object@VCF=TetradVCF
            .Object@DfOfNCO=data.frame()
            .Object@DfOfCO=data.frame()            
            
            return(.Object)
          }
)

#' Wrapper function Tetrad
#' @details This method allow to merge the VCF files from the 4 individus of a tetrad. Merging is done by 
#' doing a strict intersection based on the variant position.
#' @param ListOfVCFLP A list of VCFlikeParser containing VCFlike of sublclass SM
#' @return An object of class [\code{\link{Tetrad}}]
#' @name Tetrad
#' @rdname Tetrad-class
#' @export
Tetrad <- function(ListOfVCFLP) new("Tetrad",ListOfVCFLP=ListOfVCFLP)


#' @aliases setDfOfNCO, Tetrad-method
#' @param value a data.frame
#' @rdname Tetrad-class
#' @name Tetrad-class

setReplaceMethod(f="setDfOfNCO", 
                   signature="Tetrad", 
                   definition=function(object,value){ 
                     object@DfOfNCO <- value
                     return(object) 
                   } ) 

#' @title FindNCO Method for the class [\code{\link{Tetrad}}]
#' @details
#' To find NCO events, the following encoding has been used:
#' \describe{
#' \item{1}{for the homozygous Col genotype}
#' \item{0}{for the heterozygous Col/Ler genotype.}
#' }
#' At tetrad variant site, we sum this numbers and got:
#' \describe{
#' \item{4}{if the 4 individuals are heterozygous}
#' \item{0}{if the 4 individuals are homozygous}
#' These two cases don’t have to occur.
#' \item{3}{if one individual is homozygous and the 3 others are heterozygous}
#' \item{1}{if one individual is heterozygous and the 3 others are homozygous}
#' This occurs for gene conversion events. We respectively called these events for being of type 3:1 or of type 1:3.
#' \item{2}{if two individuals are homozygous and the 2 others are heterozygous. Most of cases, it is a 2:2 segregation.
#' We searched for gene conversion events by looking for 3* and 1* in the string representing the tetrad genotype along the variant. 
#' We reported the positions of the first and last converted variant sites together 
#' with the leftmost and rightmost non-converted sites. We reported the two tracks length.}
#' }
#' @param object An object of class Tetrad
#' @return An object of class [\code{\linkS4class{Tetrad}}] with the DfOfNCO slots.
#' @aliases FindNCO, Tetrad-method
#' @rdname Tetrad-class
#' @name Tetrad-class
#' @exportMethod FindNCO

setMethod(f="FindNCO",
          signature(object="Tetrad"),        
          definition = function(object){
            
            # Changer l'encodage de la matrice des génotypes
            GENOM <- geno(object@VCF)$GT
            GENOM[which(GENOM=="0/1")]<-1
            GENOM[which(GENOM=="0/0")]<-0
            GENOM[which(GENOM=="1/1")]<-100
            GENOM[which(GENOM=="./.")]<-1000
            GENOM=data.frame(apply(GENOM,2,function(x){as.numeric(as.vector(x))}))
            GENO <- data.frame("Chrom"= as.vector(rowRanges(object@VCF)@seqnames),"colPos"=start(rowRanges(object@VCF)@ranges),"M1"=as.numeric(GENOM[,1]),"M2"=as.numeric(GENOM[,2]),"M3"=as.numeric(GENOM[,3]),"M4"=as.numeric(GENOM[,4]))
            #GENO <- GENO[order(GENO$Chrom,GENO$colPos),]
            
            TetradGenotypeByChrom <- split(GENO,as.factor(GENO$Chrom))
  
            FindNCOByChrom<-function(TetradGenotypeByChrom){  
              rS <- rowSums(TetradGenotypeByChrom[,c(3:6)])
              
              rSc <- as.character(rS)
              rSc[which(rS>=100 & rS <1000)]<-"F"
              rSc[which(rS>=1000)]<-"N"
              
              bigS = paste(rSc,collapse="")
          
              # Type 1:3
              res1=str_locate_all(bigS,"11*")
              if(sum(unlist(lapply(res1,length)))!=0){
                a1=apply(res1[[1]],1,function(x){
                  substr(bigS,x[1]-10,x[2]+10)
                })
                t1=c(res1[[1]][,2]-res1[[1]][,1])+1
                inter1Conv=TetradGenotypeByChrom$colPos[res1[[1]][,2]]-TetradGenotypeByChrom$colPos[res1[[1]][,1]]
                inter1Flanq = TetradGenotypeByChrom$colPos[res1[[1]][,2]+1]-TetradGenotypeByChrom$colPos[res1[[1]][,1]-1]
                
                df1=data.frame(rep("NCO1",length(t1)),rep(TetradGenotypeByChrom[1,1],length(t1)),TetradGenotypeByChrom$colPos[res1[[1]][,1]],TetradGenotypeByChrom$colPos[res1[[1]][,1]-1],TetradGenotypeByChrom$colPos[res1[[1]][,2]],TetradGenotypeByChrom$colPos[res1[[1]][,2]+1],t1,inter1Conv,inter1Flanq,a1)
                names(df1)=c("Type","Chrom","S (conv)", "S (flanq)", "E (conv)", "E (flanq)","# markers","L (conv)", "L (flanq)","Context (+/- 10 markers)")
              }
              else{
                warning(paste("No NCO of type 1:3 detected for chromosom",TetradGenotypeByChrom[1,1],sep=""))
                df1=data.frame()
              }
              
              # Type 3:1
              res3=str_locate_all(bigS,"33*")
              if(sum(unlist(lapply(res3,length)))!=0){
                a3=apply(res3[[1]],1,function(x){
                  substr(bigS,x[1]-10,x[2]+10)
                })
                t3=c(res3[[1]][,2]-res3[[1]][,1])+1
                inter3Conv=TetradGenotypeByChrom$colPos[res3[[1]][,2]]-TetradGenotypeByChrom$colPos[res3[[1]][,1]]
                inter3Flanq = TetradGenotypeByChrom$colPos[res3[[1]][,2]+1]-TetradGenotypeByChrom$colPos[res3[[1]][,1]-1]
                
                df3=data.frame(rep("NCO3",length(t3)),rep(TetradGenotypeByChrom[1,1],length(t3)),TetradGenotypeByChrom$colPos[res3[[1]][,1]],TetradGenotypeByChrom$colPos[res3[[1]][,1]-1],TetradGenotypeByChrom$colPos[res3[[1]][,2]],TetradGenotypeByChrom$colPos[res3[[1]][,2]+1],t3,inter3Conv,inter3Flanq,a3)
                names(df3)=c("Type","Chrom","S (conv)", "S (flanq)", "E (conv)", "E (flanq)","# markers","L (conv)", "L (flanq)","Context (+/- 10 markers)")
              }
              else{
                warning(paste("No NCO of type 3:1 detected for chromosom",TetradGenotypeByChrom[1,1],sep=""))
                df3=data.frame()
              }
              #print(df3)
              return(list("NCO1"=df1, "NCO3"=df3))
            }
            
            NCO <- lapply(TetradGenotypeByChrom,FindNCOByChrom)
            
            dfNCO <- do.call("rbind",lapply(NCO,function(x){
              do.call("rbind",x)
            }))
            
            setDfOfNCO(object)<-dfNCO
            return(object)
          }
          
)  


#' @aliases setDfOfCO, Tetrad-method
#' @param value a data.frame
#' @rdname Tetrad-class
#' @name Tetrad-class

setReplaceMethod(f="setDfOfCO", 
                   signature="Tetrad", 
                   definition=function(object,value){ 
                     object@DfOfCO <- value
                     return(object) 
                   } )



#' @title FindCO Method for the class [\code{\link{Tetrad}}]
#' @details
#' This method does a system call to the python program CrossOver.py
#' which aims to find CO. Parameters 20K5K5K
#' \describe{
#' \item{1}{Export the tetrad genotype in a inputFile}
#' \item{2}{Find and write CO in outputfiles}
#' }
#' @param object An object of class Tetrad
#' @return An object of class [\code{\linkS4class{Tetrad}}] with the DfOfCO slots.
#' @aliases FindCO, Tetrad-method
#' @rdname Tetrad-class
#' @name Tetrad-class
#' @exportMethod FindCO

setMethod(f="FindCO",
          signature(object="Tetrad"),        
          definition = function(object, InputFileWoExt, TetradNameWithCoParam, InputDirPath, OutputDirPath){

                                        # Write InputFile
             genotype <- geno(object@VCF)$GT

                                        # Undetermined genotype
             UG <- as.vector(sort(c(which(genotype[,1] %in% c("./.","1/1")),which(genotype[,2] %in% c("./.","1/1")),which(genotype[,3] %in% c("./.","1/1") ),which(genotype[,4] %in% c("./.","1/1")))))
             
                                        # Construct the dataframe to export as Input to CrossOver
             genotype <- genotype[-UG,]
             genotype[which(genotype=="0/1")] <- 1
             genotype[which(genotype=="0/0")] <- 0
             Pos <- ranges(rowRanges(object@VCF))@start[-UG]
             Chrom <- as.factor(seqnames(rowRanges(object@VCF)))[-UG]
             levels(Chrom)=1:5
             
                                        # Export the dataframe
             InputToCO <- data.frame("Chrom"=Chrom,"Pos"=Pos,"M1"=genotype[,1],"M2"=genotype[,2],"M3"=genotype[,3],"M4"=genotype[,4]) 

             if(! file.exists(InputDirPath)){
                 stopifnot(dir.create(InputDirPath, showWarnings = TRUE, recursive = FALSE, mode = "0770"))
             }
             tmp<-paste(InputDirPath,"/",InputFileWoExt,"_tmp.txt",sep="")
             write.table(InputToCO,sep="\t",file=tmp,quote=FALSE,row.names=FALSE,col.names=FALSE)
             system(paste("awk \'{print $1 \"\t\" $2 \"\t\t\" $3 \"\t\" $4 \"\t\" $5 \"\t\" $6 \"\r\"}\'",tmp,">",paste(InputDirPath,"/",InputFileWoExt,".txt",sep="")))
                    
                                        # Verify that OutputDirPath does exist, create dir if needed
             
              if(! file.exists(OutputDirPath)){
                 stopifnot(dir.create(OutputDirPath, showWarnings = TRUE, recursive = FALSE, mode = "0770"))
             }
              if(! file.exists(paste(OutputDirPath,TetradNameWithCoParam,sep="/"))){
                 stopifnot(dir.create(paste(OutputDirPath,TetradNameWithCoParam,sep="/"), showWarnings = TRUE, recursive = FALSE, mode = "0770"))
             }
             
                                        # Run Crossovers
             #status <- system(paste("python /usr/local/bin/crossOver_adaptation.py", InputFileWoExt ,  TetradNameWithCoParam , InputDirPath , OutputDirPath, sep= " "), intern=FALSE)
	     status <- system(paste("python /save/project/ijpb/bioinfo-code/src/Recombine_v2.1/CrossOver_v6.2/crossOver_adaptation.py", InputFileWoExt , TetradNameWithCoParam , InputDirPath , OutputDirPath, sep= " "), intern=FALSE)

                                        # Test the output
             stopifnot(status == 0)

             COList <- read.table(paste(OutputDirPath,"/",TetradNameWithCoParam,"/CoList_",TetradNameWithCoParam,".txt",sep=""),skip=1,sep="\t")[,-c(1,8,9)]
             names(COList) <- c("Chr","pos(bp)","co_type","chromatids","GCasso_tract_len", "GCasso_marker")
             COList$Chr <- paste("Chr",COList$Chr,sep="")
             
             setDfOfCO(object)<- COList[which(COList$co_type %in% c(0,1)),]
            return(object) 
         }
        )


#' plot Method for the class [\code{\link{Tetrad}}]
#' @aliases plot, Tetrad-method
#' @rdname Tetrad-class

setMethod(f="plot",
          signature(x="Tetrad"),        
          definition = function(x, Chromosom=c(1:5),addNCO=FALSE, addCO=FALSE, addGenomicFeatures=TRUE, ...){
            
            #load("/projects/Single_Meiosis/Ressources/Ler0/variantutils/data/ColGenomeFeatures.Rdata")
	     load("/save/dcharif/bioinfo-dev/r-dev/variantutils/data/ColGenomeFeatures.Rdata")     

            VariantPos <- start(rowRanges(x@VCF)@ranges)
            VariantChrom <-  as.vector(rowRanges(x@VCF)@seqnames)
            VariantGenotype <- geno(x@VCF)$GT
            centro <- ColGenomeFeatures[[2]]$centromere
            hetero <- ColGenomeFeatures[[2]]$heterochromatine
            
            par(mfrow=c(3,2))
            # Chromosom iteration
            for(j in Chromosom){
              Pos <- VariantPos[which(VariantChrom==paste("Chr",j,sep=""))]
              Geno <- VariantGenotype[which(VariantChrom==paste("Chr",j,sep="")),]
              plot(c(1, ColGenomeFeatures[[1]][[j]]) ,c(1,6),type="n",xlab="Positions",ylab="Individu",axes=FALSE,main=names(ColGenomeFeatures[[1]][j]))
              axis(2, at = 1:4, labels = c("M1","M2","M3","M4"),cex.axis=0.8)
              axis(1)
              for(i in 1:4){
                points(Pos[which(Geno[,i]=="0/1")],rep(i,length(Pos[which(Geno[,i]=="0/1")])),col="red",cex=0.9,pch="|")
                points(Pos[which(Geno[,i]=="0/0")],rep(i,length(Pos[which(Geno[,i]=="0/0")])),col="black",cex=0.3,pch="|")
              }
              rect(centro[which(centro$V1==names(ColGenomeFeatures[[1]][j])),2],0.5,centro[which(centro$V1==names(ColGenomeFeatures[[1]][j])),3],4.5,col="yellow",density=80,border=NA)
              rect(hetero[which(hetero$V1==names(ColGenomeFeatures[[1]][j])),2],0.5,hetero[which(hetero$V1==names(ColGenomeFeatures[[1]][j])),3],4.5,col="grey",density=30,border=NA)
              
              legend(x=0.75*max(Pos),y=6,c("Col","Het"),col=c(1,2),pch="|")
            }
          })

#' Write Method for the class [\code{\link{Tetrad}}]
#' @aliases writeTetrad, Tetrad-method
#' @rdname Tetrad-class

setMethod(f="writeTetrad",
          signature(object="Tetrad"),        
          definition = function(object, VCFtetradFile, CoFile, NcoFile, ReportFile, sdiFile, ...){

	      # test if sdiFile exist
	      stopifnot(file.exists(sdiFile))

              sdiLer0 = read.table(sdiFile,h=TRUE)

              # Export the tetrad
              writeVcf(object@VCF, file=VCFtetradFile)

              # Export the DfOfNCO
              # Merge Annotation
              NCO=merge(object@DfOfNCO,sdiLer0,by.x=c("Chrom","S (conv)"),by.y=c("Chrom","colPosNew"))
              # Merge Genotype
              id <-pmatch(paste(NCO[,1],"_",NCO[,2],"_",sep=""),row.names(geno(object@VCF)$GT)) 
              
              FileNCOGL <- cbind(NCO,geno(object@VCF)$GT[id,])

              tmp <- geno(object@VCF)$GL[id,]

              f <- function(x){
                  tmp <- lapply(x,function(x){
                      unlist(strsplit(as.character(x),split=c(",")))
                  })
                  names(tmp)=NULL
                  df <- as.data.frame(t(as.data.table(tmp)))[,1:2]
                  names(df) <- c("pCol","pHet")
                  df
              }

              M <- lapply(1:4, function(i){ f(tmp[,i]) })
              names(M) <- paste("M",1:4,sep="")

              tmpM <- as.data.frame(M)
              
              FileNCOGT <- cbind(FileNCOGL,tmpM)
              
              FileNCO <- cbind(FileNCOGT,info(object@VCF)[id,grep("DP",names(info(object@VCF)))])
              
              write.table(FileNCO,file=NcoFile,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
              write.csv(FileNCO,file=paste(NcoFile,".csv",sep=""))
              
              # Export the DfOfCO
              write.table(object@DfOfCO,file=CoFile,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

              # Export informations in the report file
 
              tmp <-as.data.frame(rbind( table(geno(VCFTetrad2@VCF)$GT[,1]), table(geno(VCFTetrad2@VCF)$GT[,2]), table(geno(VCFTetrad2@VCF)$GT[,3]) ,table(geno(VCFTetrad2@VCF)$GT[,4])))
             row.names(tmp) <- paste("M",1:4,sep="")
 
              cat(paste("Nombre de sites pour la tétrade: ",dim(geno(object@VCF)$GT)[1],sep=" "),"\n",file=ReportFile)
              write.table(tmp,file=ReportFile,append=TRUE)
          })



###############
### BamCounter
###############

setMethod("initialize",
          signature(.Object="BamCounter"),
          function(.Object, file, param){ 
            .Object@file = BamFile(file=file) 
            .Object@param = param
            .Object@res = scanBam(file=.Object@file, param=.Object@param)[[1]] 
            validObject(.Object) 
			return(.Object)
          }
)

setMethod(
		  f = "countMismatches",
		  signature = "BamCounter",
		  definition = function(object, mismatchTag="NM",tags){
			 if (!(is.null(tags))) {
				 if (!(class(mismatchTag)=="character") && !(nchar(mismatchTag)==2))
					stop("provided 'mismatchTag' value is not a 2-character string")
				 if (!(is(tags,"character") && any(nchar(tags)>0)))
					stop("tags must be a character vector and each tag element must be non-empty string")

				# cross NM tag with given other tags
				df <- as.data.frame(xtabs(as.formula(paste(paste("~",mismatchTag,sep=""),tags, sep="+",collapse="+")), data=object@res[["tag"]]))
				return(df)
			 } else {
				# only NM tag
				df <- as.data.frame(xtabs(paste("~",mismatchTag,sep=""), data=object@res[["tag"]]))
				return(df)
			 }
		  })

setMethod(
		  f = "countPrimaryTag",
		  signature = "BamCounter",
		  definition = function(object, primaryTag="IE", tags){
			if (!(is.null(tags))) {
                 if (!(class(primaryTag)=="character") && !(nchar(primaryTag)==2))
                    stop("provided 'primaryTag' value is not a 2-character string")
                 if (!(is(tags,"character") && any(nchar(tags)>0)))
                    stop("tags must be a character vector and each tag element must be non-empty string")
				 if (is.null(object@res[["tag"]][[primaryTag]]))
					stop(paste(primaryTag, "tag does not exist in BamCounter object", sep=" "))

				 # cross primaryTag with given other tags
				 df <- as.data.frame(xtabs(as.formula(paste(paste("~",primaryTag,sep=""),tags, sep="+",collapse="+")), data=object@res[["tag"]]))
				 df
			} else {
				# only primaryTag
				df <- as.data.frame(xtabs(paste("~",primaryTag,sep=""), data=object@res[["tag"]]))
				df
			}
		  })

setReplaceMethod(
				 f = "setCounts",
				 signature = "BamCounter",
				 definition = function(object,value){
					 object@counts = value
					 validObject(object)
					 return(object)
				 }
				)




###########
### JGI-aa
##########

#' Constructor method fof the [\code{\link{JGI-aa-like}}] Class.
#' @param File A string giving the path to the VCF like file.
#' @param RefGenome A string giving the reference genome
#' @param Sample A string giving the name of the sample
#' @seealso VCFlike-class, VCFlikeParser-class
#' @name JGI-aa-class
#' @rdname JGI-aa-class

setMethod("initialize",
          signature(.Object="JGI-aa"),
          function(.Object, File, RefGenome,Sample){
            .Object@RefGenome=RefGenome
            .Object@Sample=Sample
            .Object@File=File
            .Object@Df=read.table(File)
            
            dfNames <- c(
              "contig",
              "pos",
              "type",
              "name",
              "strand",
              "ref_nt",
              "cds_nt_pos",
              "ref_aa:codon", 
              "cds_aa_pos",
              "cons_qual",
              "depth",
              "avg_hits",
              "cons_base",
              "1st_base_count",
              "2nd_base_count",
              "cons_aa:codon"
              )
            
            names(.Object@Df) <- dfNames
            
            ContigTmp <- as.factor(.Object@Df$contig)
            tmp <- TestChromFormat(RefGenome = RefGenomes[[RefGenome]] , ChromosomeNames = levels(ContigTmp)) 
            
            if( ! is.null(tmp)){
              if(any(is.na(tmp$UserFileChromName))){
              levels(ContigTmp)<- tmp$RefGenomChromName[na.omit(match(levels(ContigTmp),tmp$UserFileChromName))]
                
              }
              else{
                levels(ContigTmp)<- tmp$RefGenomChromName[match(levels(ContigTmp),tmp$UserFileChromName)]
              }
            }
            
            .Object@Df$contig <- as.vector(ContigTmp)
            
            
            .Object@MetaDf <- DataFrame(
              Number=rep(".",16),
              Type=c("String","Integer","String","String","String","String","String","String","String","String","Integer","Double","String","Integer","Integer","String"),
              Description=c(
                "contig =  contig/chromosome/scaffold name",
                "pos = genomic position in contig/chromosome/scaffold",
                "type = effect",
                "name = gene name", 
                "strand = strand",
                "ref_nt = reference base",
                "cds_nt_pos = cds nucleotide position",
                "ref_aa:codon = reference amino acid:reference codon",
                "cds_aa_pos = cds codon position",
                "cons_qual = phred-like consensus quality",
                "depth = read depth",
                "avg_hits = the average number of hits in reference",
                "cons_base = consensus base",
                "1st_base:count = most common allele count",
                "2nd_base:count = 2nd most common allele count",
                "cons_aa:codon = consensus amino acid:consensus codon"
              ))
            
            dimnames(.Object@MetaDf)[[1]] <- dfNames
            
            #.Object@VCF = NULL
            
            return(.Object)
          }
)

#' ConvertToVcf method for the JGI-aa class
#' @details This function is called by the class VCFlikeParser to set the slot VCF of an object of class JGI-aa. 
#' @param VCFlike A VCFlike object
#' @return VCF
#' @seealso VCFlike-class, VCFlikeParser-class 
#' @name JGI-aa-class
#' @rdname JGI-aa-class

setMethod(
  f ="ConvertToVcf",
  signature ="JGI-aa",
  definition = function( object ){
    
    infodf = DataFrame(
      "DP"=as.character(object@Df$depth)
    )
    
    CONS_BASE_TMP <- lapply( IUPAC_CODE_MAP[as.character(object@Df$cons_base)],s2c)
    names(CONS_BASE_TMP) <- NULL
    CONS_BASE_TMP_FOR_ID <- unlist(lapply(CONS_BASE_TMP,function(x){
      paste(x,collapse=",")
    }))
    
    fixeddf=DataFrame(
      "REF"=DNAStringSet(as.character(gsub("-","",object@Df$ref_nt))),
      "ALT"=DNAStringSetList(CONS_BASE_TMP),
      "QUAL"= as.numeric(object@Df$cons_qual),
      "FILTER"=as.character(rep("PASS",dim(object@Df)[1])))

    
    header=VCFHeader(reference= object@RefGenome, samples= object@Sample, header=DataFrameList("INFO"=object@MetaDf[c(11),]))
    
    VCF=VCF(
      rowRanges=GRanges(seqnames = Rle(object@Df$contig),IRanges(start=object@Df$pos,width=rep(1,dim(object@Df)[1]),names=paste(object@Df$contig,"_",object@Df$pos,"_",object@Df$ref_nt,"_",CONS_BASE_TMP_FOR_ID,sep=""))),
      fixed=fixeddf,
      info=infodf,
      exptData=SimpleList(header=header),
      collapsed=TRUE)
    
    genome(seqinfo(rowRanges(VCF))) <- object@RefGenome
    
    
    return(VCF)
  }
)




###########
### JGI-indel
##########

#' Constructor method fof the [\code{\link{JGI-indel-like}}] Class.
#' @param File A string giving the path to the VCF like file.
#' @param RefGenome A string giving the reference genome
#' @param Sample A string giving the name of the sample
#' @seealso VCFlike-class, VCFlikeParser-class
#' @name JGI-indel-class
#' @rdname JGI-indel-class

setMethod("initialize",
          signature(.Object="JGI-indel"),
          function(.Object, File, RefGenome,Sample){
            .Object@RefGenome=RefGenome
            .Object@Sample=Sample
            .Object@File=File
            .Object@Df=read.table(File)
  
            
            # Split indel type and allele
            tmp <- as.data.frame(t(as.data.frame(lapply(as.character(.Object@Df$V7), strsplit,split=":"))))
            .Object@Df$indel_size <- as.numeric(as.vector(tmp$V1))
            .Object@Df$indel_seq <- as.vector(tmp$V2)
            .Object@Df$SVTYPE <- ifelse(.Object@Df$indel_size<0,"DEL","INS")
            
            dfNames <- c(
              "type",
              "name",
              "contig",
              "pos",
              "evidence",
              "depth",
              "size_seq",
              "nb_read_fwd", 
              "nb_read_rev",
              "up_flank",
              "down_flank",
              "nb_read_wo_indel",
              "indel_size",
              "indel_seq",
              "SVTYPE"
            )   
            names(.Object@Df) <- dfNames
            
            ContigTmp <- as.factor(.Object@Df$contig)
            tmp <- TestChromFormat(RefGenome = RefGenomes[[RefGenome]] , ChromosomeNames = levels(ContigTmp)) 
           
            if( ! is.null(tmp)){
              if(any(is.na(tmp$UserFileChromName))){
                levels(ContigTmp)<- tmp$RefGenomChromName[na.omit(match(levels(ContigTmp),tmp$UserFileChromName))]
                
              }
              else{
                levels(ContigTmp)<- tmp$RefGenomChromName[match(levels(ContigTmp),tmp$UserFileChromName)]
              }
            }
            
            .Object@Df$contig <- as.vector(ContigTmp)
            
            .Object@MetaDf <- DataFrame(
              Number=rep(".",15),
              Type=c("String","String","String","Integer","String","Integer","String","Integer","Integer","String","String","Integer","String","String","String"),
              Description=c(
                "type = effect",
                "name = gene name",
                "contig =  contig/chromosome/scaffold name",
                "pos = genomic position in contig/chromosome/scaffold",
                "evidence = evidence for indel",
                "depth = read depth",
                "size_seq  = size-of-indel:sequence-of-indel",
                "nb_read_fwd = number of forward reads",
                "nb_read_rev = number of reverse reads",
                "up_flank   = upstream flanking sequence",
                "down_flank = downstream flanking sequence",
                "nb_read_wo_indel = number of reads crossing this site without indel",
                "indel_size = size-of-indel",
                "indel_seq = sequence-of-indel",
                "SVTYPE = Type of structural variant"  
              ))
            
            dimnames(.Object@MetaDf)[[1]] <- dfNames
            
            return(.Object)
          }
)

#' ConvertToVcf method for the JGI-indel class
#' @details This function is called by the class VCFlikeParser to set the slot VCF of an object of class JGI-indel. 
#' @param VCFlike A VCFlike object
#' @return VCF
#' @seealso VCFlike-class, VCFlikeParser-class 
#' @name JGI-indel-class
#' @rdname JGI-indel-class

setMethod(
  f ="ConvertToVcf",
  signature ="JGI-indel",
  definition = function( object ){
    
    LENGTH<-abs(object@Df$indel_size) 
    REF<-rep(NA,length(LENGTH))
    ALT<-rep(NA,length(LENGTH))
    BASES <- rep(NA,length(LENGTH))
    #Fetch the base before the event
    
    BASES <- FetchBaseBeforEvent(object,c(1:length(BASES)),"contig","pos")
  
    REF[which(object@Df$indel_size < 0)] <- apply(data.frame(BASES[which(object@Df$indel_size < 0)], as.character(object@Df$indel_seq[which(object@Df$indel_size < 0)])),1,paste,collapse="")
    ALT[which(object@Df$indel_size < 0)] <- BASES[which(object@Df$indel_size < 0)]
    
    REF[which(object@Df$indel_size > 0)] <- BASES[which(object@Df$indel_size > 0)]
    ALT[which(object@Df$indel_size > 0)] <- apply(data.frame(BASES[which(object@Df$indel_size > 0)], as.character(object@Df$indel_seq[which(object@Df$indel_size > 0)])),1,paste,collapse="")
     
    infodf = DataFrame(
      "DP"=as.character(object@Df$depth),
      "SVTYPE"=object@Df$SVTYPE
    )
 
    
    fixeddf=DataFrame(
      "REF"=DNAStringSet(as.character(REF)),
      "ALT"=VariantAnnotation:::.formatALT(as.list(as.character(ALT))), 
      "QUAL"= as.numeric(rep(NA,dim(object@Df)[1])),
      "FILTER"=as.character(rep("PASS",dim(object@Df)[1])))
    
    header=VCFHeader(reference= object@RefGenome, samples= object@Sample, header=DataFrameList("INFO"=object@MetaDf[c(6,15),]))
    
    VCF=VCF(
      rowRanges=GRanges(seqnames = Rle(object@Df$contig),IRanges(start=object@Df$pos-1,width=nchar(gsub("-","",REF))+1,names=paste(object@Df$contig,"_",object@Df$pos-1,"_",REF,"_",ALT,sep=""))),
      fixed=fixeddf,
      info=infodf,
      exptData=SimpleList(header=header),
      collapsed=TRUE)
    
    genome(seqinfo(rowRanges(VCF))) <- object@RefGenome
    
    
    return(VCF)
  }
)



###########
### sdi-7
##########

#' Constructor method fof the [\code{\link{sdi-7-like}}] Class.
#' @param File A string giving the path to the VCF like file.
#' @param RefGenome A string giving the reference genome
#' @param Sample A string giving the name of the sample
#' @seealso VCFlike-class, VCFlikeParser-class
#' @name sdi-7-class
#' @rdname sdi-7-class

setMethod("initialize",
          signature(.Object="sdi-7"),
          function(.Object, File, RefGenome,Sample){
            .Object@RefGenome=RefGenome
            .Object@Sample=Sample
            .Object@File=File
            .Object@Df=read.table(File,sep="\t",fill=TRUE,na.strings=c("NA",""))
            
            # Add SVTYPE
            .Object@Df$SVTYPE <- ifelse(.Object@Df[,3]==0,"SUB",ifelse(.Object@Df[,3]<0,"DEL","INS"))
            
            dfNames <- c(
              "CHROMOSOME",
              "POSITION",
              "LENGTH",
              "REFERENCE_BASE",
              "CONSENSUS_BASE",
              "QUALITY_VALUE",
              "PERCENTAGE_VALUE",
              "SVTYPE"
            )
            
            names(.Object@Df) <- dfNames
            
            ContigTmp <- as.factor(.Object@Df$CHROMOSOME)
            tmp <- TestChromFormat(RefGenome = RefGenomes[[RefGenome]] , ChromosomeNames = levels(ContigTmp)) 
            
 
            if( ! is.null(tmp)){
              if(any(is.na(tmp$UserFileChromName))){
                levels(ContigTmp)<- tmp$RefGenomChromName[na.omit(match(levels(ContigTmp),tmp$UserFileChromName))]
                
              }
              else{
                levels(ContigTmp)<- tmp$RefGenomChromName[match(levels(ContigTmp),tmp$UserFileChromName)]
              }
            }
            
            .Object@Df$CHROMOSOME <- as.vector(ContigTmp)
            
            .Object@MetaDf <- DataFrame(
              Number=rep(".",10),
              Type=c("String","Integer","Double","String","String","Integer","Integer","Double","String","String"),
              Description=c(
                "chromosome: The same name as the chromosome id in reference fasta file",
                "position: 1-based leftmost position",
                "length: the length difference of the changed sequence against reference (0 for SNPs, negative for deletions, positive for insertions)",
                "reference base: [-A-Z]+ (regular expression range)",
                "consensus base: [-A-Z]+ (regular expression range), IUPAC code is used for heterozygous sites",
                "quality_value: Phred quality score for SNP only",
                "percentage_value: Concordance",
                "SVTYPE: Type of structural variant"
              ))
            
            dimnames(.Object@MetaDf)[[1]] <- dfNames
            
            return(.Object)
          }
)

#' ConvertToVcf method for the sdi-7 class
#' @details This function is called by the class VCFlikeParser to set the slot VCF of an object of class sdi-7. 
#' @param VCFlike A VCFlike object
#' @return VCF
#' @seealso VCFlike-class, VCFlikeParser-class 
#' @name sdi-7-class
#' @rdname sdi-7-class

setMethod(
  f ="ConvertToVcf",
  signature ="sdi-7",
  definition = function( object ){
    
    
    infodf = DataFrame(
      "LENGTH"=as.character(object@Df$LENGTH),
      "SVTYPE"=as.character(object@Df$SVTYPE)
    )
    
    # Convert IUPAC code
    ALT=as.character(object@Df$CONSENSUS_BASE)
    ALTlist=as.list(ALT)
    ALTlist[which(ALT %in%  names(IUPAC_CODE_MAP))] <- lapply(IUPAC_CODE_MAP[ALT[which(ALT %in%  names(IUPAC_CODE_MAP))]],s2c)
    ALTlistForID <-  unlist(lapply(ALTlist,function(x){
      paste(x,collapse=",")
    }))
    ALTlist[which(ALT %in% "-")] <- ""
    
    fixeddf=DataFrame(
      "REF"=DNAStringSet(as.character(gsub("-","",object@Df$REFERENCE_BASE))),
      "ALT"=DNAStringSetList(ALTlist), 
      "QUAL"= as.numeric(object@Df$QUALITY_VALUE),
      "FILTER"=as.character(rep("PASS",dim(object@Df)[1])))
    
    width=object@Df$LENGTH
    width[which(width<0)]<- c(-1)
    
    header=VCFHeader(reference= object@RefGenome, samples= object@Sample, header=DataFrameList("INFO"=object@MetaDf[c(3,10),]))
    
    VCF=VCF(
      rowRanges=GRanges(seqnames = Rle(object@Df$CHROMOSOME),IRanges(start=object@Df$POSITION,width=as.numeric(width+1),names=paste(object@Df$CHROMOSOME,"_",object@Df$POSITION,"_",object@Df$REFERENCE_BASE,"_",ALTlistForID,sep=""))),
      fixed=fixeddf,
      info=infodf,
      exptData=SimpleList(header=header),
      collapsed=TRUE)
    
    genome(seqinfo(rowRanges(VCF))) <- object@RefGenome
    
    
    return(VCF)
  }
)


###########
### VAST
##########

#' Constructor method fof the [\code{\link{VAST-like}}] Class.
#' @param File A string giving the path to the VCF like file.
#' @param RefGenome A string giving the reference genome
#' @param Sample A string giving the name of the sample
#' @seealso VCFlike-class, VCFlikeParser-class
#' @name VAST-class
#' @rdname VAST-class

setMethod("initialize",
          signature(.Object="VAST"),
          function(.Object, File, RefGenome,Sample){
            .Object@RefGenome=RefGenome
            .Object@Sample=Sample
            .Object@File=File
            .Object@Df=read.table(File,sep="\t",fill=TRUE,na.strings=c("NA",""))
            
            .Object@Df[,3] <- as.character(.Object@Df[,3] )
            .Object@Df[,4] <- as.character(.Object@Df[,4] )
            
            # Add SVTYPE
            .Object@Df$SVTYPE <- ifelse(nchar(.Object@Df[,3])==nchar(.Object@Df[,4]),"SUB",ifelse(nchar(.Object@Df[,3]) < nchar(.Object@Df[,4]),"INS","DEL"))
            
            dfNames <- c(
              "CHROM",
              "POS",
              "REF",
              "ALT",
              "QUAL",
              "SVTYPE"
            )
            
            names(.Object@Df) <- dfNames
            
            ContigTmp <- as.factor(.Object@Df$CHROM)
            tmp <- TestChromFormat(RefGenome = RefGenomes[[RefGenome]] , ChromosomeNames = levels(ContigTmp)) 
            
            if( ! is.null(tmp)){
              if(any(is.na(tmp$UserFileChromName))){
                levels(ContigTmp)<- tmp$RefGenomChromName[na.omit(match(levels(ContigTmp),tmp$UserFileChromName))]
                
              }
              else{
                levels(ContigTmp)<- tmp$RefGenomChromName[match(levels(ContigTmp),tmp$UserFileChromName)]
              }
            }
            
            .Object@Df$CHROM <- as.vector(ContigTmp)
            
            .Object@MetaDf <- DataFrame(
              Number=rep(".",6),
              Type=c("String","Integer","String","String","Double","String"),
              Description=c(
                "CHROM: The same name as the chromosome id in reference fasta file",
                "POS: 1-based leftmost position",
                "REF: the length difference of the changed sequence against reference (0 for SNPs, negative for deletions, positive for insertions)",
                "ALT: [-A-Z]+ (regular expression range)",
                "QUAL: [-A-Z]+ (regular expression range)",
                "SVTYPE: Type of structural variant"
              ))
            
            dimnames(.Object@MetaDf)[[1]] <- dfNames
            
            return(.Object)
          }
)

#' ConvertToVcf method for the VAST class
#' @details This function is called by the class VCFlikeParser to set the slot VCF of an object of class VAST. 
#' @param VCFlike A VCFlike object
#' @return VCF
#' @seealso VCFlike-class, VCFlikeParser-class 
#' @name VAST-class
#' @rdname VAST-class

setMethod(
  f ="ConvertToVcf",
  signature ="VAST",
  definition = function( object ){
    
    
    infodf = DataFrame(
      "SVTYPE"=as.character(object@Df$SVTYPE)
    )
    
    fixeddf=DataFrame(
      "REF"=DNAStringSet(as.character(gsub("-","",object@Df$REF))),
      "ALT"=DNAStringSetList(as.list(object@Df$ALT)), 
      "QUAL"= as.numeric(rep(NA,dim(object@Df)[1])),
      "FILTER"=as.character(rep("PASS",dim(object@Df)[1])))
    
    header=VCFHeader(reference= object@RefGenome, samples= object@Sample, header=DataFrameList("INFO"=object@MetaDf[c(6),]))
    
    VCF=VCF(
      rowRanges=GRanges(seqnames = Rle(object@Df$CHROM),IRanges(start=object@Df$POS,width=nchar(object@Df$ALT),names=paste(object@Df$CHROM,"_",object@Df$POS,"_",object@Df$REF,"_",object@Df$ALT,sep=""))),
      fixed=fixeddf,
      info=infodf,
      exptData=SimpleList(header=header),
      collapsed=TRUE)
    
    genome(seqinfo(rowRanges(VCF))) <- object@RefGenome
    
    
    return(VCF)
  }
)




###########
### SALK
##########

#' Constructor method fof the [\code{\link{SALK-like}}] Class.
#' @param File A string giving the path to the VCF like file.
#' @param RefGenome A string giving the reference genome
#' @param Sample A string giving the name of the sample
#' @seealso VCFlike-class, VCFlikeParser-class
#' @name SALK-class
#' @rdname SALK-class

setMethod("initialize",
          signature(.Object="SALK"),
          function(.Object, File, RefGenome,Sample){
            .Object@RefGenome=RefGenome
            .Object@Sample=Sample
            .Object@File=File
            .Object@Df=read.table(File,sep="\t",fill=TRUE)
            
            # Add SVTYPE
            #.Object@Df$SVTYPE <- ifelse(.Object@Df[,3]==0,"SUB",ifelse(.Object@Df[,3]<0,"DEL","INS"))
              
            
            dfNames <- c(
              "NAME",
              "CHROMOSOME",
              "POSITION",
              "REFERENCE_BASE",
              "SUBSTITUTION_BASE",
              "QUALITY",
              "NUMBER_OF_NONREPETITIVE_READS",
              "CONCORDANCE",
              "AVERAGE_HITS"
            )
            
            names(.Object@Df) <- dfNames
            
            ContigTmp <- as.factor(.Object@Df$CHROMOSOME)
            tmp <- TestChromFormat(RefGenome = RefGenomes[[RefGenome]] , ChromosomeNames = levels(ContigTmp)) 
            
            if( ! is.null(tmp)){
              if(any(is.na(tmp$UserFileChromName))){
                levels(ContigTmp)<- tmp$RefGenomChromName[na.omit(match(levels(ContigTmp),tmp$UserFileChromName))]
                
              }
              else{
                levels(ContigTmp)<- tmp$RefGenomChromName[match(levels(ContigTmp),tmp$UserFileChromName)]
              }
            }
            
            .Object@Df$CHROMOSOME <- as.vector(ContigTmp)
            
            .Object@MetaDf <- DataFrame(
              Number=rep(".",9),
              Type=c("String","String","Integer","String","String","Integer","Integer","Double","Double"),
              Description=c(
                "Name: The name of the strain",
                "Chromosome: The same name as the chromosome id in reference fasta file",
                "Position: The position",
                "Reference_base: The reference base",
                "Substitution_base: The alternative base [ACGT] for substitutions, [-] for 1-bp deletions.",
                "Quality: Phred quality score",
                "Nb_of_NonRepetitive_reads: Number of nonrepetitive reads supporting substitution.",
                "Concordance: Ratio of reads supporting a predicted feature to total coverage (excluding quality masked bases).",
                "Avg_Hits: Average number of alignments of all reads covering this genomic position."
              ))
            
            dimnames(.Object@MetaDf)[[1]] <- dfNames
            
            return(.Object)
          }
)

#' @title ConvertToVcf Methof for the class  [\code{\link{SALK}}]
#' @details This function is called by the class VCFlikeParser to set the slot VCF of an object of class  [\code{\link{SALK}}]. 
#' @param VCFlike A VCFlike object
#' @return VCF
#' @name SALK-class
#' @rdname SALK-class

setMethod(
  f ="ConvertToVcf",
  signature ="SALK",
  definition = function( object ){
    
    infodf = DataFrame(
      "NUMBER_OF_NONREPETITIVE_READS"=as.character(object@Df$NUMBER_OF_NONREPETITIVE_READS),
      "CONCORDANCE"=as.character(object@Df$CONCORDANCE),
      "AVERAGE_HITS"=as.character(object@Df$AVERAGE_HITS)
    )
    
    REF <- as.character(object@Df$REFERENCE_BASE)
      
    ALT <- as.character(object@Df$SUBSTITUTION_BASE)
    ALT[which(ALT %in% "-")] <- ""
    
    INDELindex <- which(object@Df$SUBSTITUTION_BASE == "-")
    INDELBASES = rep(NA,length(INDELindex))
    
    INDELBASES <- FetchBaseBeforEvent(object, INDELindex,"CHROMOSOME","POSITION")
    
    REF[INDELindex] <- apply(cbind(INDELBASES,REF[INDELindex]),1,paste,collapse="")
    ALT[INDELindex] <- apply(cbind(INDELBASES,unlist(ALT[INDELindex])),1,paste,collapse="")
    
    POSITION <- object@Df$POSITION
    POSITION[INDELindex] <- POSITION[INDELindex] - 1
    
    fixeddf=DataFrame(
      "REF"=DNAStringSet(as.character(REF)),
      "ALT"=DNAStringSetList(as.list(ALT)), 
      "QUAL"= as.numeric(object@Df$QUALITY),
      "FILTER"=as.character(rep("PASS",dim(object@Df)[1])))
       
    header=VCFHeader(reference= object@RefGenome, samples= object@Sample, header=DataFrameList("INFO"=object@MetaDf[c(7,8,9),]))
    
    VCF=VCF(
      rowRanges=GRanges(seqnames = Rle(object@Df$CHROMOSOME),IRanges(start=POSITION,width=rep(1,dim(object@Df)[1]),names=paste(object@Df$CHROMOSOME,"_",POSITION,"_",REF,"_",ALT,sep=""))),
      fixed=fixeddf,
      info=infodf,
      exptData=SimpleList(header=header),
      collapsed=TRUE)
    
    genome(seqinfo(rowRanges(VCF))) <- object@RefGenome
    
    return(VCF)
  }
)

