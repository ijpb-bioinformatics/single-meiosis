### ======================================
### Utils functions
### ======================================

clusterMapCountMismatches <- function(bamcounter_list,mismatchTag,crosstags_list) {
	# check validity of arguments
	bclen <- length(bamcounter_list)
	if (bclen==0L)
		stop("bamcounter list length must be greater than 0")
	if (bclen != length(crosstags_list))
		stop("bamcounter object and bam tags lists must be of same length")
	if (!(any(llply(bamcounter_list,class)=="BamCounter"))) 
		stop("all elements in bamcounter list must be of class BamCounter")
	if (!(class(mismatchTag)=="character") && !(nchar(mismatchTag)==2))
		                stop("provided 'mismatchTag' value is not a 2-character string")
	if (!(any(llply(crosstags_list, function(t){
											  if (is.null(class(t))) return(TRUE) 
											  ifelse(class(t)!="character", FALSE, TRUE)}
					)==TRUE)))
		stop("all elements in crosstags list must be NULL or of class character")

	# parallel clusterMap
	limit_cores=detectCores()/6
	cs <- ifelse(bclen<=limit_cores, bclen, limit_cores)
	cl <- makeCluster(cs, type="FORK")
	bamcounter_list <- BiocGenerics::clusterApplyLB(cl, 1:bclen, function(i){
																			bci <- bamcounter_list[[i]]
																			xtagi <- crosstags_list[[i]]
																			setCounts(bci) <- countMismatches(bci,mismatchTag,xtagi)
																			bci
																			})
	stopCluster(cl)
	bamcounter_list
}

clusterMapCountPrimaryTag <- function(bamcounter_list,primaryTag,crosstags_list) {
	# check validity of arguments
	bclen <- length(bamcounter_list)
	if (bclen==0L)
		stop("bamcounter list length must be greater than 0")
	if (bclen != length(crosstags_list))
		stop("bamcounter object and bam tags lists must be of same length")
	if (!(any(llply(bamcounter_list,class)=="BamCounter"))) 
		stop("all elements in bamcounter list must be of class BamCounter")
	if (!(class(primaryTag)=="character") && !(nchar(primaryTag)==2))
		                stop("provided 'primaryTag' value is not a 2-character string")
	if (!(any(llply(crosstags_list, function(t){
											  if (is.null(class(t))) return(TRUE) 
											  ifelse(class(t)!="character", FALSE, TRUE)}
					)==TRUE)))
		stop("all elements in crosstags list must be NULL or of class character")

	# parallel clusterMap
	limit_cores=detectCores()/6
	cs <- ifelse(bclen<=limit_cores, bclen, limit_cores)
	cl <- makeCluster(cs, type="FORK")
	bamcounter_list <- BiocGenerics::clusterApplyLB(cl, 1:bclen, function(i){
																			bci <- bamcounter_list[[i]]
																			xtagi <- crosstags_list[[i]]
																			setCounts(bci) <- countPrimaryTag(bci,primaryTag,xtagi)
																			bci
																			})
	stopCluster(cl)
	bamcounter_list
}

clusterMapCalculateIndependentEvents <- function(bamcounter_list,tagIE,snpIs,clusterUse) {
	# check validity of arguments
    bclen <- length(bamcounter_list)
    if (bclen==0L)
        stop("bamcounter list length must be greater than 0")
	if (!(any(llply(bamcounter_list,class)=="BamCounter")))
        stop("all elements in bamcounter list must be of class BamCounter")
	if (!(class(tagIE)=="character") && !(nchar(tagIE)==2))
		stop("provided 'tagIE' value is not a 2-character string")
	if (!(is.logical(snpIs)))
        stop("provided 'snpIs' value is not logical/boolean")
    if (!(is.logical(clusterUse)))
        stop("provided 'clusterUse' value is not logical/boolean")

	# parallel clusterMap
    limit_cores=detectCores()/6
    csie <- ifelse(bclen<=limit_cores, bclen, limit_cores)
    clie <- makeCluster(csie, type="FORK")
	bamcounter_list <- BiocGenerics::clusterApplyLB(clie, seq_len(bclen), function(i){
                                                                            bci <- bamcounter_list[[i]]
                                                                            out <- calculateIndependentEvents(bci, ieTag=tagIE, isSnp=snpIs, useCluster=clusterUse)
																			bci@res$tag=c(bci@res$tag,out)
																			bci
                                                                            })
    stopCluster(clie)
	bamcounter_list 
}

joinCounts <- function(bamcounter_list,freq_labels,by=NULL,type="left",match="first") {
	# check validity of arguments
	bclen <- length(bamcounter_list)
	if (bclen==0L)
		stop("bamcounter list length must be greater than 0")
	if (bclen != length(freq_labels))
		stop("bamcounter object and frequency labels lists must be of same length")
	if (!(any(llply(freq_labels, function(t){
					ifelse(class(t)!="character", FALSE, TRUE)}
	)==TRUE)))
		stop("all elements in frequency labels list must be of class character")
	if (!(any(llply(seq_len(bclen), function(i){
					ifelse(dim(bamcounter_list[[i]]@counts)[1]==0L,FALSE,TRUE)}
	)==TRUE)))
		stop("all dataframes in bamcounter objects must not be void")

	# rename "Freq" column name in each dataframe with the right label
	dfs <- llply(seq_len(bclen), function(i){
											  bci <- bamcounter_list[[i]]
											  colnames(bci@counts) <- replace(colnames(bci@counts), which(colnames(bci@counts) == "Freq"), freq_labels[[i]])
											  bci@counts
											}
	)
																			
	# join dataframes
	dfj <- join_all(dfs, by=by, type=type, match=match)
	dfj
}

filterTag <- function(bamcounter, tag, value, skip=FALSE) {
    # check validity of arguments
    if (!(class(bamcounter)=="BamCounter"))
        stop("provided bamcounter value is not of class BamCounter")
    bamlist<-bamcounter@res
    bamlistlen <- length(bamlist)
    if (bamlistlen==0L)
        stop("provided bam list length must be greater than 0")
    if (!(any(llply(list(tag, value), function(t){
                          if (is.null(class(t))) return(TRUE)
                          ifelse(!(class(t) %in% c("character","numeric")), FALSE, TRUE)}

        )==TRUE)))
        stop("provided tag, value pair must be NULL or of class character or numeric")
    if (!(is.logical(skip)))
        stop("provided skip value is not logical/boolean")

    # filter for tag/value
    if (!(is.null(class(tag))) && !(is.null(class(value)))) {
        if (!(is.null(bamlist$tag[[tag]]))) {
			# mandatory fields 
            fields<-names(bamlist)[which(!names(bamlist) == "tag")]

			# parallel clusterMap
			limit_cores=detectCores()/6
			flen=length(fields)
			cs <- ifelse(flen<=limit_cores, flen, limit_cores)
			cl <- makeCluster(cs, type="FORK")
            
			# output only alignments with the given tag/value pair
            if (!(skip==TRUE)) {
                # idx       
                idx<-which(bamlist$tag[[tag]]==value)

            # skip alignments with the given tag/value pair
            } else {
                # idx       
                idx<-which(bamlist$tag[[tag]]!=value)
            }
            # filter fields
            bamlistfilt <- BiocGenerics::clusterApplyLB(cl,bamlist[fields], function(v) {
                                        if (!(is.vector(v)))
                                            stop("should be a vector")
                                        v[idx]
                })
			names(bamlistfilt)<-names(bamlist[fields])
			stopCluster(cl)

            # filter tag list
			tlen=length(bamlist$tag)
			cs <- ifelse(tlen<=limit_cores, tlen, limit_cores)
            cl <- makeCluster(cs, type="FORK")
            bamlistfilt$tag <- BiocGenerics::clusterApplyLB(cl,bamlist$tag, function(v) {
                                          if (!(is.vector(v)))
                                              stop("should be a vector")
                                          v[idx]
                })
			names(bamlistfilt$tag)<-names(bamlist$tag)
            stopCluster(cl)
			bamlistfilt
        }
    }
}

countMismatchesAndJoin <- function(bam, mismatchTag="NM", withOnlyFirstHit=TRUE, by="NH"){
	# check validity of arguments
    if (!(file.exists(bam)))
        stop("provided Bam file does not exist. Check file name/path.")
	if (!(class(mismatchTag)=="character") && !(nchar(mismatchTag)==2))
		        stop("provided 'mismatchTag' value is not a 2-character string")
	if (!(is.logical(withOnlyFirstHit)))
		stop("provided withOnlyFirstHit value is not logical/boolean")
	if (!(class(by)=="character") && !(nchar(by)==2))
		stop("provided 'by' value is not a 2-character string")
 
    # count and join
	HI=ifelse(withOnlyFirstHit,"HI",NULL)
	taglist=c(mismatchTag, by, HI)
    p1=ScanBamParam(tag=taglist, what="flag")
    p2=ScanBamParam(tag=taglist, what="flag", flag=scanBamFlag(isProperPair=TRUE))
    bc1 = BamCounter(file=bam, param=p1)
    bc2 = BamCounter(file=bam, param=p2)
	if (withOnlyFirstHit) {
		bc1filt <- filterTag(bc1, "HI", 1)
    	bc2filt <- filterTag(bc2, "HI", 1)
		bc1@res<-bc1filt
		bc2@res<-bc2filt
		countColnames=c("AllAln_WithOnlyFirstHit_Freq","ProperPairAln_WithOnlyFirstHit_Freq")
		outFile=paste(dirname(path.expand(bam)),"/",basename(bam),"_countMM_",mismatchTag,"By",by,"_withOnlyFirstHit.tab",sep="")
    } else {
		countColnames=c("AllAln_WithAllHits_Freq","ProperPairAln_WithAllhHits_Freq")
		outFile=paste(dirname(path.expand(bam)),"/",basename(bam),"_countMM_",mismatchTag,"By",by,"_withAllHits.tab",sep="")
	}
	bcl=list(bc1, bc2)
	xtags=c(by)
    bcl <- clusterMapCountMismatches(bcl,mismatchTag,list(xtags,xtags))
    dfj <- joinCounts(bcl,countColnames,by=c(mismatchTag,by),type="left",match="first")

    write.table(dfj,file=outFile,sep="\t",row.names = FALSE)
}

#
# Independent events
#

countIndels <- function(cigar) {
	# check validity of arguments
	if (!(class(cigar)=="character") && grepl(pattern="([0-9]+[MIDNSHP])+|*",cigar))
		stop("provided 'cigar' value is not a character string or do not match the expected pattern")

	matches=gregexpr(pattern="([0-9]+[ID])+", text=cigar)
	indels=unlist(regmatches(cigar, m=matches))
	indelsSizeSum=sum(as.numeric(substr(indels,1,nchar(indels)-1)))
	indelsCountRef=c(length(indels),indelsSizeSum)
	indelsCountRef
}


calculateIndependentEvents <- function(bamcounter, ieTag="IE", isSnp=FALSE, useCluster=FALSE) {
    # check validity of arguments
	if (!(class(ieTag)=="character") && !(nchar(ieTag)==2))
		stop("provided 'ieTag' value is not a 2-character string")
	if (!(is.logical(isSnp)))
        stop("provided isSnp value is not logical/boolean")
	if (!(is.logical(useCluster)))
        stop("provided useCluster value is not logical/boolean")
	if (!(class(bamcounter)=="BamCounter"))
		stop("provided bamcounter value is not of class BamCounter")	
	bamlist<-bamcounter@res
    bamlistlen <- length(bamlist)
    if (bamlistlen==0L)
        stop("provided bam list length must be greater than 0")

	if (!(is.null(class(ieTag)))) {
        if (!(is.null(bamlist$cigar))) {
			# init output list
			out=list()

			if (useCluster) {
				# parallel clusterMap
            	limit_cores=detectCores()/6
				clen=length(bamlist$cigar)
				cs <- ifelse(clen<=limit_cores, clen, limit_cores)
				cl <- makeCluster(cs, type="FORK")

				idl <- BiocGenerics::clusterApplyLB(cl,bamlist$cigar, function(c) {
                                        countIndels(c)
                })                                                            			
				stopCluster(cl)                                                   			
			# count indels without clusterMap
			} else {
				idl <- llply(seq_len(length(bamlist$cigar)), function(i) {
																		countIndels(bamlist$cigar[[i]])
				})
			}	
			m=matrix(unlist(idl),ncol=2, byrow=TRUE)
			out$ID=m[,1]
			out$IS=m[,2]

			# count mismatches and independent events
			if (isSnp) {
				# check for XW tag
				if (is.null(bamlist$tag$XW))
					stop("XW tag is absent in provided alignments. Please check your bam file.")
				# calculate IE=XW+ID
				out$IE <- bamlist$tag$XW+out$ID
			} else
			{
				# mismatches: IM=NM-IS
				if (is.null(bamlist$tag$NM))
					stop("NM tag is absent in provided alignments. Please check your bam file.")

				#out$IM <- unlist(llply(seq_len(length(bamlist$tag$NM)), function(i) {
				#																	bamlist$tag$NM[i]-out$IS[i]
				#}))
				out$IM=bamlist$tag$NM-out$IS
				# independent events: IE=IM+ID	
				out$IE <- out$IM+out$ID
			}

			# output
			out
		}
	}	
}

countIndependentEventsAndJoin <- function(bam, ieTag="IE", withOnlyFirstHit=TRUE, by="NH", mmTag="XW", vrTag="XV", tagsCrossing=c("IM"),isSnp=FALSE, useCluster=FALSE, outFilePrefix=NULL){
	# check validity of arguments
    if (!(file.exists(bam)))
        stop("provided Bam file does not exist. Check file name/path.")
	if (!(class(ieTag)=="character") && !(nchar(ieTag)==2))
		        stop("provided 'ieTag' value is not a 2-character string")
	if (!(is.logical(withOnlyFirstHit)))
		stop("provided 'withOnlyFirstHit' value is not logical/boolean")
	if (!(class(by)=="character") && !(nchar(by)==2))
		stop("provided 'by' value is not a 2-character string")
	if (!(class(mmTag)=="character") && !(nchar(mmTag)==2))
        stop("provided 'mmTag' value is not a 2-character string")
 	if (!(class(vrTag)=="character") && !(nchar(vrTag)==2))
        stop("provided 'vrTag' value is not a 2-character string")
	if (!(any(llply(tagsCrossing, function(t){
                          ifelse(class(t)=="character" && nchar(t)==2, TRUE, FALSE)}

        )==TRUE)))
        stop("provided tags crossing list elements must be of class character and 2-character long")
	if (!(is.logical(isSnp)))
        stop("provided 'isSnp' value is not logical/boolean")
	if (!(is.logical(useCluster)))
        stop("provided 'useCluster' value is not logical/boolean")

    # load bam 
	HI=ifelse(withOnlyFirstHit,"HI",NULL)
	if (isSnp==FALSE) {
		taglist=c("NM", by, HI)
	} else {
		taglist=c("NM", by, mmTag, vrTag, HI)
    }
	fields=c("cigar","flag")
	p1=ScanBamParam(tag=taglist, what=fields)
    p2=ScanBamParam(tag=taglist, what=fields, flag=scanBamFlag(isProperPair=TRUE))
    bc1 = BamCounter(file=bam, param=p1)
    bc2 = BamCounter(file=bam, param=p2)

	print(taglist	)

	xtags=c(tagsCrossing, by)
	xtagsCollapsed=paste(xtags,collapse="-")
	# filter tag
	if (withOnlyFirstHit) {
		bc1@res <- filterTag(bc1, "HI", 1)
    	bc2@res <- filterTag(bc2, "HI", 1)
#		bc1@res<-bc1filt
#		bc2@res<-bc2filt
		countColnames=c("AllAln_WithOnlyFirstHit_Freq","ProperPairAln_WithOnlyFirstHit_Freq")
		outFile=ifelse(is.null(outFilePrefix),paste(dirname(path.expand(bam)),"/",basename(bam),"_count",ieTag,"By",xtagsCollapsed,"_withOnlyFirstHit.tab",sep=""),
									paste(outFilePrefix,"_count",ieTag,"By",xtagsCollapsed,"_withOnlyFirstHit.tab",sep=""))
    } else {
		countColnames=c("AllAln_WithAllHits_Freq","ProperPairAln_WithAllhHits_Freq")
		outFile=ifelse(is.null(outFilePrefix),paste(dirname(path.expand(bam)),"/",basename(bam),"_count",ieTag,"By",xtagsCollapsed,"_withAllHits.tab",sep=""),
									paste(outFilePrefix,"_count",ieTag,"By",xtagsCollapsed,"_withAllHits.tab",sep=""))
	}
	# calculating IE
	bcl=list(bc1, bc2)
	bcl <- clusterMapCalculateIndependentEvents(bcl, tagIE=ieTag,snpIs=isSnp,clusterUse=useCluster)
	# cross counting and join
    bcl <- clusterMapCountPrimaryTag(bcl,ieTag,list(xtags,xtags))
    dfj <- joinCounts(bcl,countColnames,by=c(ieTag,xtags),type="left",match="first")

    write.table(dfj,file=outFile,sep="\t",row.names = FALSE)
}

countIndependentEvents <- function(bam, ieTag="IE", withOnlyFirstHit=TRUE, by="NH", mmTag="XW", vrTag="XV", tagsCrossing=c("IM"),isSnp=FALSE, useCluster=FALSE, outFilePrefix=NULL){
	# check validity of arguments
    if (!(file.exists(bam)))
        stop("provided Bam file does not exist. Check file name/path.")
	if (!(class(ieTag)=="character") && !(nchar(ieTag)==2))
		        stop("provided 'ieTag' value is not a 2-character string")
	if (!(is.logical(withOnlyFirstHit)))
		stop("provided 'withOnlyFirstHit' value is not logical/boolean")
	if (!(class(by)=="character") && !(nchar(by)==2))
		stop("provided 'by' value is not a 2-character string")
	if (!(class(mmTag)=="character") && !(nchar(mmTag)==2))
        stop("provided 'mmTag' value is not a 2-character string")
 	if (!(class(vrTag)=="character") && !(nchar(vrTag)==2))
        stop("provided 'vrTag' value is not a 2-character string")
	if (!(any(llply(tagsCrossing, function(t){
                          ifelse(class(t)=="character" && nchar(t)==2, TRUE, FALSE)}

        )==TRUE)))
        stop("provided tags crossing list elements must be of class character and 2-character long")
	if (!(is.logical(isSnp)))
        stop("provided 'isSnp' value is not logical/boolean")
	if (!(is.logical(useCluster)))
        stop("provided 'useCluster' value is not logical/boolean")

    # load bam 
	HI=ifelse(withOnlyFirstHit,"HI",NULL)
	if (isSnp==FALSE) {
		taglist=c("NM", by, HI)
	} else {
		taglist=c("NM", by, mmTag, vrTag, HI)
    }
	fields=c("cigar","flag")
	p1=ScanBamParam(tag=taglist, what=fields)
    bc1 = BamCounter(file=bam, param=p1)

	print(taglist	)

	xtags=c(tagsCrossing, by)
	xtagsCollapsed=paste(xtags,collapse="-")
	# filter tag
	if (withOnlyFirstHit) {
		bc1@res <- filterTag(bc1, "HI", 1)
#		bc1@res<-bc1filt
		countColnames=c("AllAln_WithOnlyFirstHit_Freq")
		outFile=ifelse(is.null(outFilePrefix),paste(dirname(path.expand(bam)),"/",basename(bam),"_count",ieTag,"By",xtagsCollapsed,"_withOnlyFirstHit.tab",sep=""),
									paste(outFilePrefix,"_count",ieTag,"By",xtagsCollapsed,"_withOnlyFirstHit.tab",sep=""))
    } else {
		countColnames=c("AllAln_WithAllHits_Freq")
		outFile=ifelse(is.null(outFilePrefix),paste(dirname(path.expand(bam)),"/",basename(bam),"_count",ieTag,"By",xtagsCollapsed,"_withAllHits.tab",sep=""),
									paste(outFilePrefix,"_count",ieTag,"By",xtagsCollapsed,"_withAllHits.tab",sep=""))
	}
	# calculating IE
	bcl=list(bc1)
	bcl <- clusterMapCalculateIndependentEvents(bcl, tagIE=ieTag,snpIs=isSnp,clusterUse=useCluster)
	# cross counting and join
    bcl <- clusterMapCountPrimaryTag(bcl,ieTag,list(xtags))

	# rename "Freq" column name in each dataframe with the right label
    dfl <- llply(seq_len(length(bcl)), function(i){ 
                                              bci <- bcl[[i]]
                                              colnames(bci@counts) <- replace(colnames(bci@counts), which(colnames(bci@counts) == "Freq"), countColnames[[i]])
                                              bci@counts
                                            })

    write.table(dfl[[1]],file=outFile,sep="\t",row.names = FALSE)
}
