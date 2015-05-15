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
makeOrganismDbFromTxDb(txdb)
}
\arguments{
  \item{txdb}{a \code{TxDb} object}.
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



