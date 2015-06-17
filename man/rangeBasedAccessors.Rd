\name{rangeBasedAccessors}

\alias{transcripts,MultiDb-method}
\alias{exons,MultiDb-method}
\alias{cds,MultiDb-method}
\alias{genes,MultiDb-method}

\alias{transcriptsBy,MultiDb-method}
\alias{exonsBy,MultiDb-method}
\alias{cdsBy,MultiDb-method}

\alias{getTxDbIfAvailable}
\alias{getTxDbIfAvailable,MultiDb-method}

% new stuff:
\alias{asBED,MultiDb-method}
\alias{asGFF,MultiDb-method}

\alias{disjointExons,MultiDb-method}
\alias{microRNAs,MultiDb-method}
\alias{tRNAs,MultiDb-method}
\alias{promoters,MultiDb-method}

\alias{distance,GenomicRanges,MultiDb-method}

\alias{extractTranscriptSeqs,BSgenome-method}

\alias{extractUpstreamSeqs,MultiDb-method}

\alias{intronsByTranscript,MultiDb-method}
\alias{fiveUTRsByTranscript,MultiDb-method}
\alias{threeUTRsByTranscript,MultiDb-method}

\alias{isActiveSeq,MultiDb-method}
\alias{isActiveSeq<-,MultiDb-method}


\title{
  Extract genomic features from an object
}

\description{
  Generic functions to extract genomic features from an object.
  This page documents the methods for \link{OrganismDb} objects only.
}

\usage{
\S4method{transcripts}{MultiDb}(x, vals=NULL, columns=c("TXID", "TXNAME"))

\S4method{exons}{MultiDb}(x, vals=NULL, columns="EXONID")

\S4method{cds}{MultiDb}(x, vals=NULL, columns="CDSID")

\S4method{genes}{MultiDb}(x, vals=NULL, columns="GENEID")

\S4method{transcriptsBy}{MultiDb}(x, by, columns, use.names=FALSE,
                                     outerMcols=FALSE)

\S4method{exonsBy}{MultiDb}(x, by, columns, use.names=FALSE, outerMcols=FALSE)

\S4method{cdsBy}{MultiDb}(x, by, columns, use.names=FALSE, outerMcols=FALSE)

\S4method{getTxDbIfAvailable}{MultiDb}(x, ...)
% \S4method{columns}{MultiDb}(x)

% new stuff: (replace TxDb with MultiDb)
\S4method{asBED}{MultiDb}(x)
\S4method{asGFF}{MultiDb}(x)

\S4method{disjointExons}{MultiDb}(x, aggregateGenes=FALSE, 
              includeTranscripts=TRUE, ...) 
\S4method{microRNAs}{MultiDb}(x)
\S4method{tRNAs}{MultiDb}(x)
\S4method{promoters}{MultiDb}(x, upstream=2000, downstream=200, ...)

\S4method{distance}{GenomicRanges,MultiDb}(x, y, ignore.strand=FALSE,
    ..., id, type=c("gene", "tx", "exon", "cds"))

\S4method{extractTranscriptSeqs}{BSgenome}(x, transcripts, strand = "+")

\S4method{extractUpstreamSeqs}{MultiDb}(x, genes, width=1000, exclude.seqlevels=NULL)

\S4method{intronsByTranscript}{MultiDb}(x, use.names=FALSE)
\S4method{fiveUTRsByTranscript}{MultiDb}(x, use.names=FALSE)
\S4method{threeUTRsByTranscript}{MultiDb}(x, use.names=FALSE)

\S4method{isActiveSeq}{MultiDb}(x)
}

\arguments{ 
  \item{x}{
    A \link{MultiDb} object. Except for the extractTranscriptSeqs
  method.  In that case it's a \link[BSgenome]{BSgenome} object and the
  second argument is an \link{MultiDb} object.
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
  \item{use.names}{Controls how to set the names of the returned
    \link[GenomicRanges]{GRangesList} object.
    These functions return all the features of a given type (e.g.
    all the exons) grouped by another feature type (e.g. grouped by
    transcript) in a \link[GenomicRanges]{GRangesList} object.
    By default (i.e. if \code{use.names} is \code{FALSE}), the
    names of this \link[GenomicRanges]{GRangesList} object
    (aka the group names) are the internal ids of the features
    used for grouping (aka the grouping features), which are
    guaranteed to be unique.
    If \code{use.names} is \code{TRUE}, then the names of the
    grouping features are used instead of their internal ids.
    For example, when grouping by transcript (\code{by="tx"}),
    the default group names are the transcript internal ids
    (\code{"tx_id"}). But, if \code{use.names=TRUE}, the group
    names are the transcript names (\code{"tx_name"}).
    Note that, unlike the feature ids, the feature names are not
    guaranteed to be unique or even defined (they could be all
    \code{NA}s). A warning is issued when this happens.
    See \code{?\link{id2name}} for more information about
    feature internal ids and feature external names and how
    to map the formers to the latters.

    Finally, \code{use.names=TRUE} cannot be used when grouping
    by gene \code{by="gene"}. This is because, unlike for the
    other features, the gene ids are external ids (e.g. Entrez
    Gene or Ensembl ids) so the db doesn't have a \code{"gene_name"}
    column for storing alternate gene names.
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
    \item{upstream}{
    For \code{promoters} : An \code{integer(1)} value indicating the
    number of bases upstream from the transcription start site. For
    additional details see \code{?`promoters,GRanges-method`}.
  }
  \item{downstream}{
    For \code{promoters} : An \code{integer(1)} value indicating the 
    number of bases downstream from the transcription start site. For 
    additional details see \code{?`promoters,GRanges-method`}.
  }
  \item{aggregateGenes}{
    For \code{disjointExons} : A \code{logical}. When \code{FALSE}
    (default) exon fragments that overlap multiple genes are dropped.
    When \code{TRUE}, all fragments are kept and the \code{gene_id}
    metadata column includes all gene ids that overlap the exon fragment.
  }
  \item{includeTranscripts}{
    For \code{disjointExons} : A \code{logical}. When \code{TRUE}
    (default) a \code{tx_name} metadata column is included that 
    lists all transcript names that overlap the exon fragment. 
  }
  \item{y}{For \code{distance}, a \link{MultiDb} instance. The
   \code{id} is used to extract ranges from the \link{MultiDb}
   which are then used to compute the distance from \code{x}.
  }
  \item{id}{A \code{character} vector the same length as \code{x}.
    The \code{id} must be identifiers in the \link{MultiDb} object.
    \code{type} indicates what type of identifier \code{id} is. 
  }
  \item{type}{A \code{character(1)} describing the \code{id}. 
    Must be one of \sQuote{gene}, \sQuote{tx}, \sQuote{exon} or
    \sQuote{cds}.
  }
  \item{ignore.strand}{
    A \code{logical} indicating if the strand of the ranges
    should be ignored. When \code{TRUE}, strand is set to \code{'+'}.
  }
  \item{outerMcols}{
    A \code{logical} indicating if the the 'outer' mcols (metadata
    columns) should be populated for some range based accesors which
    return a GRangesList object.  By default this is FALSE, but if TRUE
    then the outer list object will also have it's metadata columns
    (mcols) populated as well as the mcols for the 'inner' GRanges
    objects.
  }
  \item{transcripts}{
    An object representing the exon ranges of each transcript to extract.
    It must be a \link[GenomicRanges]{GRangesList} or \link{MultiDb}
    object while the \code{x} is a \link[BSgenome]{BSgenome} object.
    Internally, it's turned into a \link[GenomicRanges]{GRangesList}
    object with \code{\link{exonsBy}(transcripts, by="tx", use.names=TRUE)}.
  }
  \item{strand}{
    Only supported when \code{x} is a \link[Biostrings]{DNAString} object.

    Can be an atomic vector, a factor, or an \link[S4Vectors]{Rle} object,
    in which case it indicates the strand of each transcript (i.e. all the
    exons in a transcript are considered to be on the same strand).
    More precisely: it's turned into a factor (or factor-\link[S4Vectors]{Rle})
    that has the "standard strand levels" (this is done by calling the
    \code{\link[BiocGenerics]{strand}} function on it). Then it's recycled
    to the length of \link[IRanges]{RangesList} object \code{transcripts}
    if needed. In the resulting object, the i-th element is interpreted
    as the strand of all the exons in the i-th transcript.

    \code{strand} can also be a list-like object, in which case it indicates
    the strand of each exon, individually. Thus it must have the same
    \emph{shape} as \link[IRanges]{RangesList} object \code{transcripts}
    (i.e. same length plus \code{strand[[i]]} must have the same length
    as \code{transcripts[[i]]} for all \code{i}).

    \code{strand} can only contain \code{"+"} and/or \code{"-"} values.
    \code{"*"} is not allowed.
  }
  \item{genes}{
    An object containing the locations (i.e. chromosome name, start, end, and
    strand) of the genes or transcripts with respect to the reference genome.
    Only \link[GenomicRanges]{GenomicRanges} and \link{MultiDb} objects
    are supported at the moment. If the latter, the gene locations are obtained
    by calling the \code{\link{genes}} function on the \link{MultiDb}
    object internally.
  }
  \item{width}{
    How many bases to extract upstream of each TSS (transcription start site).
  }
  \item{exclude.seqlevels}{
    A character vector containing the chromosome names (a.k.a. sequence levels)
    to exclude when the genes are obtained from a \link{MultiDb} object.
  }

}

\details{
  These are the range based functions for extracting transcript information
  from a \link{MultiDb} object. 
}

\value{ a GRanges or GRangesList object }

\author{
  M. Carlson
}

\seealso{
  \itemize{
    \item \link[OrganismDbi]{MultiDb-class} for how to use the
          simple "select" interface to extract information from a
          \code{MultiDb} object.
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

## Get the TxDb from an MultiDb object (if it's available)
getTxDbIfAvailable(Homo.sapiens)

## Other functions listed above should work in way similar to their TxDb
## counterparts.  So for example:
promoters(Homo.sapiens)
## Should give the same value as:
promoters(getTxDbIfAvailable(Homo.sapiens))


}

\keyword{methods}
