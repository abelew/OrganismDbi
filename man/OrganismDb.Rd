\name{OrganismDb-class}

\alias{OrganismDb-class}
\alias{class:OrganismDb}
\alias{OrganismDb}

% select and select related methods
\alias{columns,OrganismDb-method}
\alias{keytypes,OrganismDb-method}
\alias{keys,OrganismDb-method}
\alias{select,OrganismDb-method}
\alias{metadata,OrganismDb-method}
\alias{dbconn,OrganismDb-method}
\alias{dbfile,OrganismDb-method}

\alias{seqinfo,OrganismDb-method}
% \alias{seqinfo<-,OrganismDb-method}
% \alias{seqnameStyle,OrganismDb-method}
% \alias{seqnameStyle<-,OrganismDb-method}

\title{OrganismDb objects}

\description{
  The OrganismDb class is a container for storing knowledge
  about existing Annotation packages and the relationships between these
  resources.  The purpose of this object and it's associated methods is
  to provide a means by which users can conveniently query for data from
  several different annotation resources at the same time using a
  familiar interface.

  The supporting methods \code{select}, \code{columns} and \code{keys} are
  used together to extract data from an \code{OrganismDb}
  object in a manner that should be consistent with how these are used
  on the supporting annotation resources.

  The family of \code{seqinfo} style getters (\code{seqinfo},
  \code{seqlevels}, \code{seqlengths}, \code{isCircular}, \code{genome},
  and \code{seqnameStyle}) is also supported for OrganismDb objects
  provided that the object in question has an embedded TxDb
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
      \code{keys(x, keytype, pattern, column, fuzzy)}: Return keys for
      the database contained in the \link{TxDb} object .

      The \code{keytype} argument specifies the kind of keys that will
      be returned and is always required.

      If \code{keys} is used with \code{pattern}, it will pattern match
      on the \code{keytype}.

      But if the \code{column} argument is also provided along with the
      \code{pattern} argument, then \code{pattern} will be matched
      against the values in \code{column} instead.

      If \code{keys} is called with \code{column} and no \code{pattern}
      argument, then it will return all keys that have corresponding
      values in the \code{column} argument.

      Thus, the behavior of \code{keys} all depends on how many arguments are
      specified.

      Use of the \code{fuzzy} argument will toggle fuzzy matching to
      TRUE or FALSE.  If \code{pattern} is not used, fuzzy is ignored.
    }    
    \item{}{
      \code{columns(x)}:
      shows which kinds of data can be returned for the
      \code{OrganismDb} object.
    }
    \item{}{
      \code{select(x, keys, columns, keytype)}:
      When all the appropriate arguments are specifiedm \code{select}
      will retrieve the matching data as a data.frame based on
      parameters for selected \code{keys} and \code{columns} and
      \code{keytype} arguments.
    }
}
}

\author{Marc Carlson}

\seealso{ 
  \itemize{
    \item \link[AnnotationDbi]{AnnotationDb-class} for more descriptsion
    of methods \code{select},\code{keytypes},\code{keys} and \code{columns}.
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
columns <- columns(Homo.sapiens)[c(7,10,11,12)]
keys <- head(keys(org.Hs.eg.db, "ENTREZID"))
keytype <- "ENTREZID"
res <- select(Homo.sapiens, keys, columns, keytype)
head(res)

## Get the DB connections or DB file paths associated with those for
## each.
dbconn(Homo.sapiens)
dbfile(Homo.sapiens)
}

