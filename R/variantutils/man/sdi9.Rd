\name{sdi9}
\alias{sdi9}
\docType{data}
\title{
Data frame of the first twenty sequence variants of Oy-9 relative to
TAIR10.
}
\description{
Data frame of the first twenty sequence variants of Oy-9 relative to
TAIR10.
The IMR-DENOM software was used to call this sequence variants and the
format is in \emph{sdi-9}.
}
\usage{data(sdi9)}
\format{
  A data frame with 20 observations on the following 9 variables.
  \describe{
    \item{\code{V1}}{\emph{chromosome}: a factor with levels \code{Chr1} \code{Chr2}
      \code{Chr3} \code{Chr4} \code{Chr5} giving the name of the chromosome}
    \item{\code{V2}}{\emph{position}: a numeric vector giving the 1-base leftmost position}
    \item{\code{V3}}{\emph{length}: the length difference of the changed sequence against reference
      \describe{
	\item{ \code{0}}{SNPs}
	\item{ \code{negative}}{deletions}
	\item{ \code{positive}}{insertions}
      }
    }
    \item{\code{V4}}{\emph{reference_base}: reference base [-A-Z]+ (regular expression range)}
    \item{\code{V5}}{\emph{consensus_base}: consensus base [-A-Z]+ (regular expression range), IUPAC code is used for heterozygous sites}
    \item{\code{V6}}{\emph{Detection_score}: Detection score 
      \describe{
	\item{ \code{1:}}{ detected by both IMR and DENOM }
	\item{ \code{2:}}{ detected by IMR as heterozygous SNP but detected by
	DENOM as homozygous SNP }
	\item{ \code{4:}}{ detected by IMR only }
	\item{ \code{5:}}{ detected by DENOM only }
      }
    }
    \item{\code{V7}}{\emph{Phred_score}: Phred quality score for SNP
      only.}
    \item{\code{V8}}{\emph{HMQ}: HMQ score for SNP only.}
    \item{\code{V9}}{\emph{Consensus_base}: Consensus base when
      considering reads with HMQ only.}
  }
}

\details{
SNPs data comes from the resources relating to the 19 genomes of Arabidopsis thaliana
\url{http://mus.well.ox.ac.uk/19genomes/variants.SDI/ler_0.v7c.sdi}
  which contribute to the 1001 Arabidopsis genomes project.
}
\references{
  Multiple reference genomes and transcriptomes for Arabidopsis thaliana (Nature 2011).
}
\examples{
data(sdi9)
}

