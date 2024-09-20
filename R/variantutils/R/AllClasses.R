#' @docType class
#' @title The VCFlike class
#' @description
#' Containers for holding data from Variant Call Format like files.
#' @details 
#' The VCFlike class is a virtual class extended from \code{\link{CollapsedVCF}}. 
#' This class can not be instanciated directly.  
#' @slot VCF An object of class collapsedVCF
#' @slot File A character giving the path of the VCFlike file
#' @slot Df A data.frame storing the data of the VCFlike
#' @slot MetaDf A DataFrame object describing all the columns of the data.frame storing the VCFlike data
#' @slot RefGenome A character giving the reference genome
#' @slot Sample A character giving the name of the sample
#' @seealso \code{\link{VCFlikeParser}} for the existing subclasses.
#' @name VCFlike-class
#' @rdname VCFlike-class
#' @exportClass VCFlike

setClass("VCFlike",
         contains=c("VIRTUAL"),
         representation(
           VCF="CollapsedVCF",
           File="character",
           Df="data.frame",
           MetaDf="DataFrame",
           RefGenome="character",
           Sample="character"
         )
)


#' @docType class
#' @title The VCFlikeParser class 
#' @details 
#' The VCFlikeParser allows to instanciate several subclasses of VCFlike.
#' The subclasses (or VCF like input format) implemented are: 
#' \describe{
#' \item{sdi-9}{The \code{\link{sdi-9}} format corresponding to the class of the same name}
#' \item{SM}{The \code{\link{SM}} format corresponding to the class of the same name}
#' }
#' @slot VCFlike An object of class \code{\link{VCFlike}}
#' @slot Format The format of the input VCF like file which is also the name of the VCFlike subclass
#' @name VCFlikeParser-class
#' @rdname VCFlikeParser-class
#' @exportClass VCFlikeParser

setClass("VCFlikeParser",
         representation(
           VCFlike="VCFlike",
           Format="character"
         )
)

#' @docType class
#' @title The sdi-9 Class 
#' @description 
#' This class is a subclass of the \code{\link{VCFlike}} class and it don't have any additionaly slot. 
#' @details 
#' The sdi format is the format used by the IMR-DENOM software to store sequence variations.
#' We call it sdi-9 because it specified the format of a text file with 9 columns which are:
#'  \describe{ 
#'    \item{\emph{chromosome}:}{The name of the chromosome}
#'    \item{\emph{position}:}{The 1-base leftmost position}
#'    \item{\emph{length}:}{The length difference of the changed sequence against reference
#'      \describe{
#'        \item{\code{0}}{SNPs}
#'        \item{\code{negative}}{deletions}
#'        \item{\code{positive}}{insertions}
#'      }
#'    }
#'    \item{\emph{reference_base}:}{reference base [-A-Z]+ (regular expression range)}
#'    \item{\emph{consensus_base}:}{consensus base [-A-Z]+ (regular expression range), IUPAC code is used for heterozygous sites}
#'    \item{\emph{Detection_score}:}{Detection score 
#'      \describe{
#'        \item{ \code{1:}}{ detected by both IMR and DENOM }
#'        \item{ \code{2:}}{ detected by IMR as heterozygous SNP but detected by DENOM as homozygous SNP }
#'        \item{ \code{4:}}{ detected by IMR only }
#'        \item{ \code{5:}}{ detected by DENOM only }
#'      }
#'    }
#'    \item{\emph{Phred_score}:} {Phred quality score for SNP only.}
#'    \item{\emph{HMQ}:} {HMQ score for SNP only.}
#'    \item{\emph{Consensus_base}:} {Consensus base when considering reads with HMQ only.}
#' }
#' @references Multiple reference genomes and transcriptomes for Arabidopsis thaliana (Nature 2011).
#' \url{http://mus.well.ox.ac.uk/19genomes/variants.SDI/readme.txt}
#' @name sdi-9-class
#' @rdname sdi-9-class
#' @exportClass sdi-9

setClass("sdi-9",
         contains="VCFlike"
)

#' @docType class
#' @title The SM Class 
#' @description This class is a subclass of the \code{\link{VCFlike}} class and it don't have any additionaly slot. 
#' @details 
#' The SM format is the format used by the pipeline devoted to the analysis of the Single Meiosis data project to store sequence variations.
#' As the particularity of the pipeline is to consider two references genome at particular synteniq position the colum are:
#'  \describe{
#'    \item{\emph{ID}:}{The ID }
#'    \item{\emph{Chrom}:}{The name of the chromosom}
#'    \item{\emph{colPos}:}{Position on Col genome as in sdi-9}
#'    \item{\emph{colPosNew}:}{Position on Col genome, VCF formated}  
#'    \item{\emph{lerPos}:}{Position on Ler genome}
#'    \item{\emph{Ref}:}{The reference base on the Col genome}
#'    \item{\emph{alt}:}{The reference base on the Ler genome} 
#'    \item{\emph{size}:}{The size of the indel}
#'    \item{\emph{colDP}:}{The sequencing depth at this position}
#'    \item{\emph{lerDP}:}{The sequencing depth at this position}
#'    \item{\emph{colQsum}:}{The sum of the base quality at this position}
#'    \item{\emph{lerQsum}:}{The sum of the base quality at this position} 
#'    \item{\emph{colMAPQsum}:}{The sum of the mapping quality of the reads covering this position }
#'    \item{\emph{lerMAPQsum}:}{The sum of the mapping quality of the reads covering this position}
#' }
#' @seealso VCFlike-class, VCFlikeParser-class 
#' @name SM-class
#' @rdname SM-class
#' @exportClass SM

setClass("SM",
         contains="VCFlike"
)

#' @docType class
#' @title The RefPol Class 
#' @description This class is a subclass of the \code{\link{VCFlike}} class and it don't have any additionaly slot. 
#' @details 
#' The RefPol format is the format used by the pipeline devoted to the analysis of the Single Meiosis data project 
#' to store sequence variations of the reference list of Polymoprhism.
#'  \describe{
#'    \item{\emph{ID}:}{The ID of the polymporphism}
#'    \item{\emph{Chrom}:}{The name of the chromosome}
#'    \item{\emph{colPos}:}{Position on Col genome as in sdi-9}
#'    \item{\emph{colPosNew}:}{Position on Col genome, VCF formated}  
#'    \item{\emph{lerPos}:}{Position on Ler genome}
#'     \item{\emph{size}:}{The size of the indel}
#'    \item{\emph{Ref}:}{The reference base on the Col genome}
#'    \item{\emph{alt}:}{The reference base on the Ler genome} 
#'    \item{\emph{Detection_score}:}{Detection score 
#'      \describe{
#'        \item{ \code{1:}}{ detected by both IMR and DENOM }
#'        \item{ \code{2:}}{ detected by IMR as heterozygous SNP but detected by DENOM as homozygous SNP }
#'        \item{ \code{4:}}{ detected by IMR only }
#'        \item{ \code{5:}}{ detected by DENOM only }
#'      }
#'    }
#'    \item{\emph{Phred_score}:} {Phred quality score for SNP only.}
#'    \item{\emph{HMQ}:} {HMQ score for SNP only.}
#'    \item{\emph{IsInTE}:} {boolean, TRUE if the polymorphism is in a TE (TAIR10, GFF)}
#'    \item{\emph{IsHeterozygous}:} {boolean, TRUE if the alt or ref allele are heterozygous}
#'     \item{\emph{IsInCNV}:} {boolean, TRUE if the polymorphism is in a CNV (Ler-0, shore)}
#'    \item{\emph{size}:}{The size of the indel}
#'    \item{\emph{colDP.F1}:}{The sequencing depth on Col genome in the F1  at this position}
#'    \item{\emph{lerDP.F1}:}{The sequencing depth on Ler genome in the F1 at this position}
#'    \item{\emph{colQsum.F1}:}{The sum of the base quality on Col genome in the F1 at this position}
#'    \item{\emph{lerQsum.F1}:}{The sum of the base quality on Ler genome in the F1 at this position} 
#'    \item{\emph{colMAPQsum.F1}:}{The sum of the mapping quality of the reads covering this position on Col genome in the F1}
#'    \item{\emph{lerMAPQsum.F1}:}{The sum of the mapping quality of the reads covering this position on Ler genome in the F1}
#'    \item{\emph{colDP.Parent}:}{The sequencing depth on Col genome in the Col Parent at this position }
#'    \item{\emph{lerDP.Parent}:}{The sequencing depth on Ler genome in the Ler Parent at this position}
#'    \item{\emph{colQsum.Parent}:}{The sum of the base quality on Col genome in the Col Parent at this position}
#'    \item{\emph{lerQsum.Parent}:}{The sum of the base quality on Ler genome in the Ler Parent at this position} 
#'    \item{\emph{colMAPQsum.Parent}:}{The sum of the mapping quality of the reads covering this position on Col genome in the Col Parent}
#'    \item{\emph{lerMAPQsum.Parent}:}{The sum of the mapping quality of the reads covering this position on Ler genome in the Ler Parent}
#' }
#' @seealso VCFlike-class, VCFlikeParser-class 
#' @name RefPol-class
#' @rdname RefPol-class
#' @exportClass RefPol

setClass("RefPol",
         contains="VCFlike"
)
#' @docType class
#' @title The class Tetrad 
#' @description Class tetrad
#' @details class tetrad
#' @slot VCF An object of class VCF
#' @slot DfOfNCO An object of type data.frame 
#' @slot DfOfCO 
#' @name Tetrad-class
#' @rdname Tetrad-class
#' @exportClass Tetrad

setClass("Tetrad",
         representation(
           VCF="CollapsedVCF",
           DfOfNCO="data.frame",
           DfOfCO="data.frame"
         )
)

setClass("BamCounter",
		representation(
			file="BamFile",
			param="ScanBamParam",
			res="list",
			counts="data.frame"
			),
		 sealed = FALSE
		)

BamCounter<-function(file=NA_character_, param=NA_ScanBamParam_){ 
    new("BamCounter", file=file, param=param)  
}


#' @docType class
#' @title The JGI-aa Class 
#' @description This class is a subclass of the \code{\link{VCFlike}} class and it don't have any additionaly slot. 
#' @details The JGI-aa format is the format used by the JGI to store sequence variations.
#'  \describe{ 
#'    \item{\emph{contig}:}{contig/chromosome/scaffold name}
#'    \item{\emph{pos}:}{genomic position in contig/chromosome/scaffold}
#'    \item{\emph{type}:}{
#'      \describe{
#'      \item{\code{NA}}{not annotated}
#'      \item{\code{NC}}{non-coding}
#'      \item{\code{Syn}}{synonymous}
#'      \item{\code{Non-Syn}}{non-synonymous}
#'      \item{\code{Nonsense}}{introducing or abolishing a stop codon}
#'      \item{\code{3\'}}{3\'UTR}
#'      \item{\code{5\'}}{5\'UTR}
#'      \item{\code{Int}}{Intronic}
#'      \item{\code{Spl}}{splice junction}
#'      }
#'    }
#' \item{\emph{name}:}{gene name or NC=non-coding or NA=not annotated}   
#' \item{\emph{strand}}{
#' \describe{
#' \item{\code{+}}{gene translated from plus strand}
#' \item{\code{-}}{gene translated from minus strand}
#' \item{\code{NC}}{non-coding}
#' \item{\code{NA}}{not annotated}
#'  }
#'  }
#' \item{\emph{ref_nt}}{reference base} 
#' \item{\emph{cds_nt_pos}}{cds nucleotide position or NC=non-coding or NA=not annotated}
#' \item{\emph{ref_aa:codon}}{reference amino acid:reference codon NC=non-coding, NA=not annotated}
#' \item{\emph{cds_aa_pos}}{cds codon position NC=non-coding, NA=not annotated}
#' \item{\emph{cons_qual}}{phred-like consensus quality (255 max value)}
#' \item{\emph{depth  }}{read depth (255 max value)}
#' \item{\emph{avg_hits  }}{the average number of hits in reference for reads covering this
#' position (essentially the \# of times the flanking region appears in the genome. 
#' zero here means there were 3 or more mismatches in close proximity, and it is too complicated for MAQ to calculate.
#' For more information please refer to documentation for Maq)}
#' \item{\emph{cons_base }}{ Consensus base}
#' \item{\emph{1st_base:count }}{ Most common allele count}
#' \item{\emph{2nd_base:count }}{ 2nd most common allele count}
#' \item{\emph{cons_aa:codon  }}{ Consensus amino acid:consensus codon ( NC=non-coding, NA=not annotated )}
#' }
#' @references \url{http://1001genomes.org/data/JGI/JGIHeazlewood2008/releases/current/TAIR10/README_JGIHeazlewood2008.html#Putative_SNV_files}
#' @name JGI-aa-class
#' @rdname JGI-aa-class
#' @exportClass JGI-aa

setClass("JGI-aa",
         contains="VCFlike"
)


#' @docType class
#' @title The JGI-indel Class 
#' @description This class is a subclass of the \code{\link{VCFlike}} class and it don't have any additionaly slot. 
#' @details 
#' The JGI-indel format is the format used by the JGI to store sequence variations.
#'  \describe{ 
#'    \item{\emph{type}:}{NC=non-coding CDS=coding 3\'=3\'UTR 5\'=5\'UTR, Int=intronic NA=not annotated}
#'    \item{\emph{name}:}{gene name, NC=non-coding, NA=not annotated}
#'    \item{\emph{contig}:}{contig/chromosome/scaffold name}
#' \item{\emph{pos}}{genomic position in contig/chromosome/scaffold} 
#' \item{\emph{evidence}}{evidence for indel, highest \"*\" (confirmed by reads from both strands), \"+\" (confirmed by reads at least 2 reads from a single strand)}
#' \item{\emph{depth}}{read depth (255 max value)}
#' \item{\emph{size:seq}}{size_seq size-of-indel:sequence-of-indel. A negative number indicates deletion,a positive number indicates insertion. 
#'  For example -3:TTT is a deletion of TTT}
#' \item{\emph{\#fwd}}{nb_read_fwd number of forward reads}
#' \item{\emph{\#rev}}{nb_read_rev number of reverse reads}
#' \item{\emph{up_flank}}{upstream flanking sequence}
#' \item{\emph{down_flank}}{downstream flanking sequence. If deletion, deleted sequence is included in the flanking sequence}
#' \item{\emph{#w/o_indel}}{nb_read_wo_indel number of reads crossing this site that don't contain the indel}
#' }
#' @references \url{http://1001genomes.org/data/JGI/JGIHeazlewood2008/releases/current/TAIR10/README_JGIHeazlewood2008.html#Putative_SNV_files}
#' @name JGI-indek-class
#' @rdname JGI-indel-class
#' @exportClass JGI-indel

setClass("JGI-indel",
         contains="VCFlike"
)



#' @docType class
#' @title The sdi-7 Class 
#' @description This class is a subclass of the \code{\link{VCFlike}} class and it don't have any additionaly slot. 
#' @details The sdi format is the format used by the IMR-DENOM software to store sequence variations.
#' We call it sdi-7 because it specified the format of a text file with 7 columns which are:
#'  \describe{ 
#'    \item{\emph{chromosome}:}{The name of the chromosome}
#'    \item{\emph{position}:}{The 1-base leftmost position}
#'    \item{\emph{length}:}{The length difference of the changed sequence against reference
#'      \describe{
#'        \item{\code{0}}{SNPs}
#'        \item{\code{negative}}{deletions}
#'        \item{\code{positive}}{insertions}
#'      }
#'    }
#'    \item{\emph{reference_base}:}{reference base [-A-Z]+ (regular expression range)}
#'    \item{\emph{consensus_base}:}{consensus base [-A-Z]+ (regular expression range), IUPAC code is used for heterozygous sites}
#'    \item{\emph{quality_value}:}{Phred quality score for SNP only ???}
#'    \item{\emph{percentage_value}:}{Concordance ?} 
#' }
#' @references Multiple reference genomes and transcriptomes for Arabidopsis thaliana (Nature 2011).
#' @name sdi-7-class
#' @rdname sdi-7-class
#' @exportClass sdi-7

setClass("sdi-7",
         contains="VCFlike"
)


#' @docType class
#' @title The VAST Class 
#' @description This class is a subclass of the \code{\link{VCFlike}} class and it don't have any additionaly slot. 
#' @details The VAST format is the format used by the VAST team to store sequence variations.
#'  \describe{ 
#'    \item{\emph{CHROM}:}{The name of the chromosome}
#'    \item{\emph{POS}:}{The position}
#'    \item{\emph{REF}:}{reference base [-A-Z]+ (regular expression range)}
#'    \item{\emph{ALT}:}{consensus base [-A-Z]+ (regular expression range)}
#'    \item{\emph{QUAL}:}{???}
#' }
#' @name VAST-class
#' @rdname VAST-class
#' @exportClass VAST

setClass("VAST",
         contains="VCFlike"
)


#' @docType class
#' @title The SALK Class 
#' @description This class is a subclass of the \code{\link{VCFlike}} class and it don't have any additionaly slot. 
#' @details The SALK format is a format derived from the SHORE pipeline output to store sequence variations. 
#' This format reports all potential SNPs and 1-bp deletions regardless of the quality of the call.
#'  \describe{ 
#'    \item{\emph{name}:}{The name of the strain}
#'    \item{\emph{Chromosome}:}{The name of the chromosome}
#'    \item{\emph{Position}:}{The position}
#'    \item{\emph{Reference_base}:}{The reference base}
#'    \item{\emph{Substitution_base}:}{The alternative base [ACGT] for substitutions, [-] for 1-bp deletions.}
#'    \item{\emph{Quality}:}{Phred quality score}
#'    \item{\emph{Nb_of_NonRepetitive_reads}:}{Number of nonrepetitive reads supporting substitution}
#'    \item{\emph{Concordance}:}{Ratio of reads supporting a predicted feature to total coverage (excluding quality masked bases)}
#'    \item{\emph{Avg_Hits}:}{Average number of alignments of all reads covering this genomic position}
#' }
#' @references See \url{http://1001genomes.org/data/Salk/releases/current/README.txt} for more details.
#' @name SALK-class
#' @rdname SALK-class
#' @exportClass SALK

setClass("SALK",
         contains="VCFlike"
)

