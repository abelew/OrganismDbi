\name{rangeBasedAccessors}

\alias{transcripts,OrganismDb-method}
\alias{exons,OrganismDb-method}
\alias{cds,OrganismDb-method}
\alias{genes,OrganismDb-method}

\alias{transcriptsBy,OrganismDb-method}
\alias{exonsBy,OrganismDb-method}
\alias{cdsBy,OrganismDb-method}

% \alias{columns}
% \alias{columns,OrganismDb-method}


\title{
  Extract genomic features from an object
}

\description{
  Generic functions to extract genomic features from an object.
  This page documents the methods for \link{OrganismDb} objects only.
}

\usage{
\S4method{transcripts}{OrganismDb}(x, vals=NULL, columns=c("TXID", "TXNAME"))

\S4method{exons}{OrganismDb}(x, vals=NULL, columns="EXONID")

\S4method{cds}{OrganismDb}(x, vals=NULL, columns="CDSID")

\S4method{genes}{OrganismDb}(x, vals=NULL, columns="GENEID")

\S4method{transcriptsBy}{OrganismDb}(x, by, columns)

\S4method{exonsBy}{OrganismDb}(x, by, columns)

\S4method{cdsBy}{OrganismDb}(x, by, columns)

% \S4method{columns}{OrganismDb}(x)
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
  \item{columns}{
    The columns or kinds of metadata that can be retrieved from the
    database.  All possible columns are returned by using the \code{columns}
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
    \item \link[GenomicFeatures]{transcriptsBy} for the original
          \code{transcriptsBy} method and related methods.
	}
}

\examples{
## extracting all transcripts from Homo.sapiens with some extra metadata
library(Homo.sapiens)
cols = c("TXNAME","SYMBOL")
res <- transcripts(Homo.sapiens, columns=cols)

## extracting all transcripts from Homo.sapiens, grouped by gene and
## with extra metadata
res <- transcriptsBy(Homo.sapiens, by="gene", columns=cols)

## list possible values for columns argument:
columns(Homo.sapiens)
}

\keyword{methods}
