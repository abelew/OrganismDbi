\name{makeOrganismDbFromTxDb}
\alias{makeOrganismDbFromTxDb}

\title{
  Make an OrganismDb object from an existing TxDb object.
}
\description{
  The \code{makeOrganismDbFromTxDb} function allows the user
  to make a \link{OrganismDb} object from an existing TxDb object.
}
\usage{
makeOrganismDbFromTxDb(txdb, keytype = NULL)
}
\arguments{
  \item{txdb}{a \code{TxDb} object}.
  \item{keytype}{By default, \code{makeOrganismDbFromTxDb} will try to
  guess this information based on the OrgDb object that is inferred to
  go with your TxDb object...  But in some instances, you may need to
  supply an over-ride and that is what this argument is for.  It is the
  column name of the ID type that your OrgDb will use as a foreign key
  when connecting to the data from the associated TxDb.  So for example,
  if you looked at the Homo.sapiens package the keytype for
  \code{org.Hs.eg.db}, would be 'ENTREZID' because that is the kind of
  ID that matches up with it's TxDb GENEID. (Because the GENEID for that
  specific TxDb is from UCSC and uses entrez gene IDs)}.
}
\details{
  \code{makeOrganismDbFromTxDb} is a convenience function that processes
  a \code{TxDb} object and pairs it up with GO.db and an appropriate
  \code{OrgDb} object to make a \code{OrganismDb} object.
  See \code{?\link{makeOrganismDbFromBiomart}} and
  \code{?\link{makeOrganismDbFromUCSC}} for a similar function that
  feeds data from either a BioMart or UCSC.
}

\value{A \link{OrganismDb} object.}

\author{
  M. Carlson and H. Pages
}

\seealso{
  \itemize{
    \item \code{\link{makeOrganismDbFromBiomart}} for convenient ways to make a
          \link{OrganismDb} object from BioMart online resources.

    \item The \link{OrganismDb} class.
  }
}

\examples{
## lets start with a txdb object
transcript_ids <- c(
    "uc009uzf.1",
    "uc009uzg.1",
    "uc009uzh.1",
    "uc009uzi.1",
    "uc009uzj.1"
)
txdbMouse <- makeTxDbFromUCSC(genome="mm9", tablename="knownGene",
                          transcript_ids=transcript_ids)

## Using that, we can call our function to promote it to an OrgDb object:
odb <- makeOrganismDbFromTxDb(txdb=txdbMouse)

columns(odb)
}



