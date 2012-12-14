\name{transcripts}

\alias{transcripts,OrganismDb-method}
\alias{exons,OrganismDb-method}
\alias{cds,OrganismDb-method}

\alias{transcriptsBy,OrganismDb-method}
\alias{exonsBy,OrganismDb-method}
\alias{cdsBy,OrganismDb-method}


\title{
  Extract genomic features from an object
}

\description{
  Generic functions to extract genomic features from an object.
  This page documents the methods for \link{OrganismDb} objects only.
}

\usage{
transcripts(x, ...)
\S4method{transcripts}{OrganismDb}(x, cols, vals=NULL, columns=c("tx_id", "tx_name"))

exons(x, ...)
\S4method{exons}{OrganismDb}(x, cols, vals=NULL, columns="exon_id")

cds(x, ...)
\S4method{cds}{OrganismDb}(x, cols, vals=NULL, columns="cds_id")


transcriptsBy(x, by=c("gene", "exon", "cds"), ...)
\S4method{transcriptsBy}{OrganismDb}(x, by, cols)

exonsBy(x, by=c("tx", "gene"), ...)
\S4method{exonsBy}{OrganismDb}(x, by, cols)

cdsBy(x, by=c("tx", "gene"), ...)
\S4method{cdsBy}{OrganismDb}(x, by, cols)
}

\arguments{ 
  \item{x}{
    A \link{TranscriptDb} object.
  }
  \item{...}{
    Arguments to be passed to or from methods.
  }
  \item{by}{
    One of \code{"gene"}, \code{"exon"}, \code{"cds"} or \code{"tx"}.
    Determines the grouping.
  }
  \item{cols}{
    The cols or kinds of metadata that can be retrieved from the
    database.  All possible cols are returned by using the \code{cols}
    method.
  }
  \item{vals}{
    Either \code{NULL} or a named list of vectors to be used to
    restrict the output. Valid names for this list are:
    \code{"gene_id"},
    \code{"tx_id"}, \code{"tx_name"},
    \code{"tx_chrom"}, \code{"tx_strand"},
    \code{"exon_id"}, \code{"exon_name"},
    \code{"exon_chrom"}, \code{"exon_strand"},
    \code{"cds_id"}, \code{"cds_name"},
    \code{"cds_chrom"}, \code{"cds_strand"} and \code{"exon_rank"}.
  }
  \item{columns}{
    Columns to include in the output.
    Must be \code{NULL} or a character vector with values in the
    above list of valid names. With the following restrictions:
    \itemize{
      \item \code{"tx_chrom"} and \code{"tx_strand"} are not allowed
      for \code{transcripts}.
      \item \code{"exon_chrom"} and \code{"exon_strand"} are not allowed
      for \code{exons}.
      \item \code{"cds_chrom"} and \code{"cds_strand"} are not allowed
      for \code{cds}.
    }
    If the vector is named, those names are used for the corresponding
    column in the element metadata of the returned object. 
  }
}

\details{
  These are the range based functions for extracting transcript information
  from a \link{OrganismDb} object. 
}

\value{ a GRanges or GRangesList object }

\author{
  M. Carlson
}

\seealso{
  \itemize{
    \item \link[OrganismDbi]{OrganismDb-class} for how to use the
          simple "select" interface to extract information from a
          \code{OrganismDb} object.
    \item \link[GenomicFeatures]{transcripts} for the original
          \code{transcripts} method and related methods.
    \item \link[GenomicFeatures]{transcriptsby} for the original
          \code{transcriptsBy} method and related methods.
	}
}

\examples{
## extracting all transcripts from Homo.sapiens with some extra metadata
library(Homo.sapiens)
cols = c("TXNAME","SYMBOL")
res <- transcripts(Homo.sapiens, cols)

## extracting all transcripts from Homo.sapiens, grouped by gene and
## with extra metadata
res <- transcriptsBy(Homo.sapiens, by="gene", cols)
}

\keyword{methods}
