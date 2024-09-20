
FilterDP <- function(x , distr=c("Normal","Negative Binomial"), side=c("right","left","bilateral"), filter=c(0.025,0.025)){
  
  
  if(distr == "Normal"){
    nor <-  c("mean"=mean(x,na.rm=TRUE),"sd"=sd(x,na.rm=TRUE))
    resFilter <- switch(side,
                        right = c("max"=qnorm(1-filter[1], nor["mean"], nor["sd"])),
                        left =  c("min"=qnorm(filter[1], nor["mean"], nor["sd"])),
                        bilateral = c("min"=qnorm(filter[1], nor["mean"], nor["sd"]),"max"=qnorm(1-filter[2], nor["mean"], nor["sd"]))
    )
  }
  
  else if(distr == "Negative Binomial"){
    nbi <- fitdistr(x, "Negative Binomial")$estimate
    resFilter <- switch(side,
                        right = c("max"=qnbinom(1-filter[1],size=nbi["size"],mu=nbi["mu"])),
                        left = c("min"=qnbinom(filter[1],size=nbi["size"],mu=nbi["mu"])),
                        bilateral = c("min"=qnbinom(filter[1],size=nbi["size"],mu=nbi["mu"]),"max"=qnbinom(1-filter[2],size=nbi["size"],mu=nbi["mu"]))
    )  
  }  
  
  resFilter <- round(resFilter)
  # Attention au DP à 0, passer à 1.
  resFilter[which(resFilter==0)] <- 1
  
  
  VecOfRes <- rep(FALSE,length(x))
  if(side=="right"){
    VecOfRes[which(x >= resFilter["max"])] = TRUE
  } else if(side=="left"){
    VecOfRes[which(x <= resFilter["min"])] = TRUE
  } else if(side=="bilateral"){
    VecOfRes[which(x <= resFilter["min"] | x >= resFilter["max"])] = TRUE
  }   
  return(list("values"=resFilter, "filter"=VecOfRes))
}

ComputeCoverageThreshold <- function(SiteCoverage , Distribution=c("Normal","Negative Binomial"), DistributionTail=c("right","left"), Prob=c(0.025)){
  
  
  if(Distribution == "Normal"){
    nor <-  c("mean"=mean(SiteCoverage,na.rm=TRUE),"sd"=sd(SiteCoverage,na.rm=TRUE))
    CoverageCutOff <- switch(DistributionTail,
                        right = c("max"=qnorm(1-Prob[1], nor["mean"], nor["sd"])),
                        left =  c("min"=qnorm(Prob[1], nor["mean"], nor["sd"]))                   
    )
  }
  
  else if(Distribution == "Negative Binomial"){
    nbi <- fitdistr(SiteCoverage, "Negative Binomial")$estimate
    CoverageCutOff <- switch(DistributionTail,
                        right = c("max"=qnbinom(1-Prob[1],size=nbi["size"],mu=nbi["mu"])),
                        left = c("min"=qnbinom(Prob[1],size=nbi["size"],mu=nbi["mu"]))
    )  
  }  
  
  CoverageCutOff <- round(CoverageCutOff,0)
 
  FilterDescription <- paste(Distribution,"low;",Prob, DistributionTail,"tail;",sep=" ")
  FilterValue <- as.numeric(as.vector(CoverageCutOff))
  
  return(list("FilterDescription"=FilterDescription,"FilterValue"=FilterValue))    
}
 

#  = c("ID","Chrom","colPos","lerPos","colRef","lerRef","size") 


 IntersectDataFrame <- function(ListOfDf , VectorOfColToKeep, VectorOfSuffixes){
   
  # Check for consistenties
   if( (length(ListOfDf) != length(VectorOfSuffixes)))
     stop("The number of dataframes must be equal to the number of suffixes")
   
   # Check that the column to keep are in all data frame
   CheckByNames <- lapply(ListOfDf, function(x){
     if( length(which( (names(x) %in% VectorOfColToKeep) == TRUE))  != length(VectorOfColToKeep) ){
      stop("One of the dataframe don't have all the column")
   }
     else{
       return(which(names(x) %in% VectorOfColToKeep))
     }
   })
   
  ListOfDfMod <- lapply(1:length(VectorOfSuffixes),function(x){
    dimnames(ListOfDf[[x]])[[2]][CheckByNames[[x]]]<-paste(names(ListOfDf[[x]])[CheckByNames[[x]]],".",VectorOfSuffixes[[x]],sep="")
    return(ListOfDf[[x]][,CheckByNames[[x]]])
    })    
   
   MergedDf <- Reduce(function(x, y){
     tmp = merge(x, y, by=0,sort=FALSE)
     row.names(tmp)=tmp$Row.names
     tmp[,-1]
   }
                   , ListOfDfMod , accumulate=F)
   return(MergedDf)
 }


IntersectIRanges <- function(ListOfIRanges){
     
  MergedIRanges <- Reduce(function(x, y){
    merge(x, y, by=c("start","end","width","names"),sort=FALSE)
  }
                     , ListOfIRanges , accumulate=F)
  IR <- IRanges(start=MergedIRanges$start,end=MergedIRanges$end,width=MergedIRanges$width,names=MergedIRanges$names)
  
  return(IR)
}

IntersectGeno <- function(ListOfGeno){
  
  MergedGeno <- Reduce(function(x, y){
    tmp = merge(x, y, by=0,sort=FALSE)
    row.names(tmp)=tmp$Row.names
    tmp[,-1]
  }
                          , lapply(ListOfGeno,as.data.frame) , accumulate=F)
  return(MergedGeno)
}


TestChromFormat <- function(RefGenome,ChromosomeNames){
 
   tmp <- ChromosomeNames %in% names(RefGenome)
   
   if( any( ! tmp) ){
     cat("ChromosomeNames\n")
     cat(ChromosomeNames)
     cat("\n")
     cat("dont fit the names of the RefGenome, please enter the names of the Chromosome in the following order:\n")
     cat(sort(names(RefGenome)))
     cat("\n")
     cat("if some chromosome are missing, please enter NA\n")
     Names <- scan(file="",what="character",n=7) 
     if( any(! Names %in% c(ChromosomeNames,NA)) | ! (length(unique(na.omit(Names)))==length(ChromosomeNames)) ){
       stop("ChromosomeNames dont fit the name of the given values")
     }
     return(data.frame("RefGenomChromName" = sort(names(RefGenome)), "UserFileChromName" = Names))
   } 
}


#FetchBaseBeforEvent <- function(RefGenome, Chromosome, Position){
#   toupper(RefGenome[[Chromosome]][Position])
#}

FetchBaseBeforEvent <- function(object, Index, ChromosomCol, PositionCol){
  tmp <- object@Df[Index,]
  tmp.split <- split(tmp,f=as.factor(tmp[,ChromosomCol]))
  tmp.res <- list()
  for(i in 1:length(tmp.split)){
    tmp.res[[i]] <- toupper(RefGenomes[[object@RefGenome]][[names(tmp.split)[i]]][tmp.split[[i]][,PositionCol]-1])
  }
  names(tmp.res)<-levels(as.factor(tmp[,ChromosomCol])) 
  unlist(tmp.res[unique(tmp[,ChromosomCol])])
}



plotGL <- function(Tetrad){

    stopifnot(class(Tetrad)=="Tetrad")

    f <- function(x){
        tmp <- lapply(x,function(x){
            unlist(strsplit(as.character(x),split=c(",")))
        })
        names(tmp)=NULL
        as.data.table(t(as.data.table(tmp)))
    }

    registerDoParallel(cores=8)
    M <- foreach(i=1:4) %dopar% f(geno(Tetrad@VCF)$GL[,i])
  
    df<-rbind(M[[1]],M[[2]],M[[3]],M[[4]])                        
    df$V1=as.numeric(df$V1)
    df$V2=as.numeric(df$V2)
    df$indiv=rep(c("M1","M2","M3","M4"),each=dim(M[[1]])[1])
    df$geno=c(geno(Tetrad@VCF)$GT[,1],
        geno(Tetrad@VCF)$GT[,2],
        geno(Tetrad@VCF)$GT[,3],
        geno(Tetrad@VCF)$GT[,4])
    
   # Elimination des 1/1
    dfF <- df[-which(df$geno=="1/1" | df$geno=="./."),]
    dfF$Proba<-ifelse(dfF$geno=="0/0",as.numeric(as.vector(dfF$V1)),as.numeric(as.vector(dfF$V2)))
    
    dfF <- as.data.frame(dfF)
    dfF.split <- split(dfF,dfF$geno)

    quantile1 <- quantile(dfF.split[[1]]$Proba,probs=c(0.005,0.01,0.025,0.05))
    data1=data.frame("values"=as.vector(quantile1),"quantile"=names(quantile1))
    data1$quantile=paste("q",data1$quantile,"=",round(data1$values,12),sep="")

    quantile2 <- quantile(dfF.split[[2]]$Proba,probs=c(0.005,0.01,0.025,0.05))
    data2=data.frame("values"=as.vector(quantile2),"quantile"=names(quantile2))
    data2$quantile=paste("q",data2$quantile,"=",round(data2$values,12),sep="")
    
    p1 <-  ggplot(data=dfF.split[[1]])+geom_histogram(aes(x=Proba)) +
        ggtitle("For hom, Distribution of p(hom)") +
            theme(plot.title = element_text(color="grey", size=10)) +
                geom_vline(data=data1,aes(xintercept=values,linetype=quantile),show_guide=TRUE)
    p2 <-  ggplot(data=dfF.split[[2]])+geom_histogram(aes(x=Proba)) +
        ggtitle("For het, Distribution of p(het)") +
            theme(plot.title = element_text(color="grey", size=10)) +
                geom_vline(data=data2,aes(xintercept=values,linetype=quantile),show_guide=TRUE)

    multiplot(p1, p2, cols=1)
}




#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
