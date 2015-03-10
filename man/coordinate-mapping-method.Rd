\name{mapToTranscripts}

\alias{mapToTranscripts,ANY,OrganismDb-method}

\title{Map range coordinates between transcripts and genome space}

\description{
  Map range coordinates between features in the transcriptome and
  genome (reference) space.

  See \code{?\link[GenomicAlignments]{mapToAlignments}} in the 
  \pkg{GenomicAlignments} package for mapping coordinates between 
  reads (local) and genome (reference) space using a CIGAR alignment.
}

\usage{
\S4method{mapToTranscripts}{ANY,OrganismDb}(x, transcripts,
          ignore.strand = TRUE, extractor.fun = transcripts, ...) 
}

\arguments{
  \item{x}{
    \code{\linkS4class{GenomicRanges}} object of positions to be mapped.
    \code{x} must have names when mapping to the genome.
  }
  \item{transcripts}{
    The \code{OrganismDb} object that will be used to extract features
  using the \code{extractor.fun}.
  }
  \item{ignore.strand}{
    When TRUE, strand is ignored in overlap operations.
  }
  \item{extractor.fun}{
    Function to extract genomic features from a \code{TxDb} object.

    Valid \code{extractor} functions:
    \itemize{
      \item{transcripts}  ## default
      \item{exons}
      \item{cds}
      \item{genes}
      \item{promoters}
      \item{disjointExons}
      \item{microRNAs}
      \item{tRNAs}
      \item{transcriptsBy}
      \item{exonsBy}
      \item{cdsBy}
      \item{intronsByTranscript}
      \item{fiveUTRsByTranscript}
      \item{threeUTRsByTranscript}
    }    
  }
  \item{\dots}{
    Additional arguments passed to \code{extractor.fun} functions. 
  }
}

\details{
  \itemize{
    \item{\code{mapToTranscripts}}{
      The genomic range in \code{x} is mapped to the local position in the
      \code{transcripts} ranges. A successful mapping occurs when \code{x}
      is completely within the \code{transcripts} range, equivalent to: 
      \preformatted{
      findOverlaps(..., type="within")
      }
      Transcriptome-based coordinates start counting at 1 at the beginning
      of the \code{transcripts} range and return positions where \code{x}
      was aligned. The seqlevels of the return object are taken from the
      \code{transcripts} object and should be transcript names. In this
      direction, mapping is attempted between all elements of \code{x} and
      all elements of \code{transcripts}.
    }
  }
}

\value{
  An object the same class as \code{x}.

  Parallel methods return an object the same shape as \code{x}. Ranges that
  cannot be mapped (out of bounds or strand mismatch) are returned as 
  zero-width ranges starting at 0 with a seqname of "UNMAPPED".

  Non-parallel methods return an object that varies in length similar to a
  Hits object. The result only contains mapped records, strand mismatch
  and out of bound ranges are not returned. \code{xHits} and 
  \code{transcriptsHits} metadata columns indicate the elements of \code{x} 
  and \code{transcripts} used in the mapping.

  When present, names from \code{x} are propagated to the output. When
  mapping to transcript coordinates, seqlevels of the output are the names 
  on the \code{transcripts} object; most often these will be transcript
  names. When mapping to the genome, seqlevels of the output are the seqlevels 
  of \code{transcripts} which are usually chromosome names.
}

\seealso{
  \itemize{
    \item \code{?\link[GenomicAlignments]{mapToAlignments}}{in the 
    \pkg{GenomicAlignments} package for methods mapping between
    reads and genome space using a CIGAR alignment.}  
    \item \code{?\link[GenomicFeatures]{mapToTranscripts}}{ in the 
    \pkg{GenomicFeatures} package for orgininal methods and their TxDb uses.}
  }
}

\author{V. Obenchain, M. Lawrence and H. Pages; ported to work with
  OrganismDbi by Marc Carlson}

\examples{

## ---------------------------------------------------------------------
## A. Basic Use
## ---------------------------------------------------------------------

library(Homo.sapiens)
x <- GRanges("chr2", 
             IRanges(c(7500, 8400, 9000), width=200, 
                     names=LETTERS[1:3])) 

## Map to transcript coordinates:
#mapToTranscripts(x, Homo.sapiens)
}

\keyword{methods}
\keyword{utilities}
% maybe I was confused about what this was supposed to do???
% trace(what ='mapToTranscripts',tracer=browser,signature=c('ANY','TxDb'))