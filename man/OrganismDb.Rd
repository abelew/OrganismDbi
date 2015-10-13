\name{MultiDb-class}

\alias{MultiDb-class}
\alias{class:MultiDb}
\alias{MultiDb}

\alias{OrganismDb-class}
\alias{class:OrganismDb}
\alias{OrganismDb}

% select and select related methods
\alias{columns,MultiDb-method}
\alias{keytypes,MultiDb-method}
\alias{keys,MultiDb-method}
\alias{select,MultiDb-method}
\alias{mapIds,MultiDb-method}
% special select enhancers for MultiDb objects
\alias{selectByRanges}
\alias{selectByRanges,MultiDb-method}
\alias{selectRangesById}
\alias{selectRangesById,MultiDb-method}


\alias{metadata,MultiDb-method}
\alias{dbconn,MultiDb-method}
\alias{dbfile,MultiDb-method}
\alias{taxonomyId,MultiDb-method}

\alias{seqinfo,MultiDb-method}

\alias{TxDb}
\alias{TxDb,OrganismDb-method}
\alias{TxDb<-}
\alias{TxDb<-,OrganismDb-method}
\alias{resources}
\alias{resources,MultiDb-method}



\title{MultiDb and OrganismDb objects}

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
  In the code snippets below, \code{x} is a OrganismDb object.

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
    \item{}{
      \code{mapIds(x, keys, columns, keytype, ..., multiVals)}:
      When all the appropriate arguments are specifiedm \code{mapIds}
      will retrieve the matching data as a vector or list based on
      parameters for selected \code{keys} and \code{columns} and
      \code{keytype} arguments.  The multiVals argument can be used to
      choose the format of the values returned.  Possible values for
      multiVals are:
      \describe{
	\item{first:}{This value means that when there are multiple matches only the 1st thing that comes back will be returned. This is the default behavior}
	\item{list:}{This will just returns a list object to the end user}
	\item{filter:}{This will remove all elements that contain multiple matches and will therefore return a shorter vector than what came in whenever some of the keys match more than one value}
	\item{asNA:}{This will return an NA value whenever there are multiple matches}
	\item{CharacterList:}{This just returns a SimpleCharacterList object}
	\item{FUN:}{You can also supply a function to the \code{multiVals} argument for custom behaviors.  The function must take a single argument and return a single value.  This function will be applied to all the elements and will serve a 'rule' that for which thing to keep when there is more than one element.  So for example this example function will always grab the last element in each result: \code{ last <- function(x){x[[length(x)]]} }
	}
      }
    }
    \item{}{
      \code{selectByRanges(x, ranges, columns, overlaps,
      ignore.strand)}: When all the appropriate arguments are specified,
      \code{selectByRanges} will return an annotated GRanges object that
      has been generated based on what you passed in to the ranges
      argument and whether that overlapped with what you specified in
      the overlaps argument. Internally this function will get
      annotation features and overlaps by calling the appropriate
      annotation methods indicated by the overlaps argument.  The value
      for overlaps can be any of: gene, tx, exons, cds, 5utr, introns or
      3utr.  The default value is 'tx' which will return to you, your
      annotated ranges based on whether the overlapped with the
      transcript ranges of any gene in the associated TxDb object based
      on the gene models it contains.  Also: the number of ranges
      returned to you will match the number of genes that your ranges
      argument overlapped for the type of overlap that you specified.
      So if some of your ranges are large and overlap several features
      then you will get many duplicated ranges returned with one for
      each gene that has an overlapping feature.  The columns values
      that you request will be returned in the mcols for the annotated
      GRanges object that is the return value for this function.
      Finally, the ignore.strand argument is provided to indicate
      whether or not findOverlaps should ignore or respect the strand.
    }
    \item{}{
      \code{selectRangesById(x, keys, columns, keytype, feature)}: When
      all the appropriate arguments are specified,
      \code{selectRangesById} will return a GRangesList object that
      correspond to gene models GRanges for the keys that you specify
      with the keys and keytype arguments. The annotation ranges
      retrieved for this will be specified by the feature argument and
      can be: gene, tx, exon or cds.  The default is 'tx' which will
      return the transcript ranges for each gene as a GRanges object in
      the list.  Extra data can also be returned in the mcols values for
      those GRanges by using the columns argument.
    }
    \item{}{
      \code{resources(x)}: shows where the db files are for resources
      that are used to store the data for the \code{OrganismDb} object.
    }
    \item{}{
      \code{TxDb(x)}: Accessor for the TxDb object of a
      \code{OrganismDb} object.  
    }
    \item{}{
      \code{TxDb(x) <- value}: Allows you to swap in an alternative TxDb
      for a given \code{OrganismDb} object.  This is most often useful
      when combined with \code{saveDb(TxDb, file)}, which returns the
      saved TxDb, so that you can save a TxDb to disc and then assign
      the saved version right into your \code{OrganismDb} object.
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
    \item \link[GenomeInfoDb]{seqlevels} .
    \item \link[GenomeInfoDb]{seqlengths} .
    \item \link[GenomeInfoDb]{isCircular} .
    \item \link[GenomeInfoDb]{genome} .    
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
res <-  mapIds(Homo.sapiens, keys=c('1','10'), column='ALIAS',
               keytype='ENTREZID',  multiVals="CharacterList")

## get symbols for ranges in question:
ranges <-  GRanges(seqnames=Rle(c('chr11'), c(2)),
                    IRanges(start=c(107899550, 108025550),
                            end=c(108291889, 108050000)), strand='*',
                    seqinfo=seqinfo(Homo.sapiens))
selectByRanges(Homo.sapiens, ranges, 'SYMBOL')

## Or extract the gene model for the 'A1BG' gene:
selectRangesById(Homo.sapiens, 'A1BG', keytype='SYMBOL')


## Get the DB connections or DB file paths associated with those for
## each.
dbconn(Homo.sapiens)
dbfile(Homo.sapiens)

## extract the taxonomyId
taxonomyId(Homo.sapiens)

##extract the resources
resources(Homo.sapiens)
}

