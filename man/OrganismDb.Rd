\name{OrganismDb-class}

\alias{OrganismDb-class}
\alias{class:OrganismDb}
\alias{OrganismDb}

% select and select related methods
\alias{cols,OrganismDb-method}
\alias{keytypes,OrganismDb-method}
\alias{keys,OrganismDb-method}
\alias{select,OrganismDb-method}
\alias{metadata,OrganismDb-method}

\alias{seqinfo,OrganismDb-method}
\alias{seqinfo<-,OrganismDb-method}
\alias{seqnameStyle,OrganismDb-method}
\alias{seqnameStyle<-,OrganismDb-method}

\title{OrganismDb objects}

\description{
  The OrganismDb class is a container for storing knowledge
  about existing Annotation packages and the relationships between these
  resources.  The purpose of this object and it's associated methods is
  to provide a means by which users can conveniently query for data from
  several different annotation resources at the same time using a
  familiar interface.

  The supporting methods \code{select}, \code{cols} and \code{keys} are
  used together to extract data from an \code{OrganismDb}
  object in a manner that should be consistent with how these are used
  on the supporting annotation resources.

  The family of \code{seqinfo} style methods (\code{seqinfo},
  \code{seqlevels}, \code{seqlengths}, \code{isCircular}, \code{genome},
  and \code{seqnameStyle}) is also supported for OrganismDb objects
  provided that the object in question has an embedded TranscriptDb
  object.
}

\section{Methods}{
  In the code snippets below, \code{x} is a OrganismDb object. For the
  metadata and show methods, there is also support for FeatureDb objects.

  \describe{
    \item{}{
      \code{keytypes(x)}:
      allows the user to discover which keytypes can be passed in to
      \code{select} or \code{keys} and the \code{keytype} argument.
    }
    \item{}{
      \code{keys(x, keytype)}:
      returns keys for the database contained in the \code{OrganismDb}
      object .  By default it will return the "TXNAME" keys for the database,
      but if used with the \code{keytype} argument, it will return the keys
      from that keytype.
    }
    \item{}{
      \code{cols(x)}:
      shows which kinds of data can be returned for the
      \code{OrganismDb} object.
    }
    \item{}{
      \code{select(x, keys, cols, keytype)}:
      When all the appropriate arguments are specifiedm \code{select}
      will retrieve the matching data as a data.frame based on
      parameters for selected \code{keys} and \code{cols} and
      \code{keytype} arguments.
    }
}
}

\author{Marc Carlson}

\seealso{ 
  \itemize{
    \item \link[OrganismDbi]{makeOrganismPackage} for functions
    used to generate an \code{OrganismDb} based package.
    \item \link[OrganismDbi]{rangeBasedAccessors} for the range based methods
    used in extracting data from a \code{OrganismDb} object.
    \item \link[GenomicRanges]{seqlevels} .
    \item \link[GenomicRanges]{seqlengths} .
    \item \link[GenomicRanges]{isCircular} .
    \item \link[GenomicRanges]{genome} .    
}
}

\examples{
## load a package that creates an OrganismDb
library(Homo.sapiens)
ls(2)
## then the methods can be used on this object.
cols <- cols(Homo.sapiens)[c(7,10,11,12)]
keys <- head(keys(org.Hs.eg.db, "ENTREZID"))
keytype <- "ENTREZID"
res <- select(Homo.sapiens, keys, cols, keytype)
head(res)
}
