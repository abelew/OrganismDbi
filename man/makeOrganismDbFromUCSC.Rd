\name{makeOrganismDbFromUCSC}
\alias{makeOrganismDbFromUCSC}

\title{
  Make a OrganismDb object from annotations available at the
  UCSC Genome Browser
}
\description{
  The \code{makeOrganismDbFromUCSC} function allows the user
  to make a \link{OrganismDb} object from transcript annotations
  available at the UCSC Genome Browser.
}
\usage{
makeOrganismDbFromUCSC(
        genome="hg19",
        tablename="knownGene",
        transcript_ids=NULL,
        circ_seqs=DEFAULT_CIRC_SEQS,
        url="http://genome.ucsc.edu/cgi-bin/",
        goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath",
        miRBaseBuild=NA)
}
\arguments{
  \item{genome}{genome abbreviation used by UCSC and obtained by
    \code{\link[rtracklayer]{ucscGenomes}()[ , "db"]}.
    For example: \code{"hg19"}.}
  \item{tablename}{name of the UCSC table containing the transcript
    annotations to retrieve. Use the \code{supportedUCSCtables} utility
    function to get the list of supported tables. Note that not all
    tables are available for all genomes.}
  \item{transcript_ids}{optionally, only retrieve transcript
    annotation data for the specified set of transcript ids.
    If this is used, then the meta information displayed for the
    resulting \link{OrganismDb} object will say 'Full dataset: no'.
    Otherwise it will say 'Full dataset: yes'.}
  \item{circ_seqs}{a character vector to list out which chromosomes
    should be marked as circular.}
  \item{url,goldenPath_url}{use to specify the location of an
    alternate UCSC Genome Browser.}
  \item{miRBaseBuild}{specify the string for the appropriate build
    Information from mirbase.db to use for microRNAs.  This can be
    learned by calling \code{supportedMiRBaseBuildValues}.  By default,
    this value will be set to \code{NA}, which will inactivate the
    \code{microRNAs} accessor.}
}
\details{
  \code{makeOrganismDbFromUCSC} is a convenience function that feeds
  data from the UCSC source to the lower level \code{\link{OrganismDb}}
  function.
  See \code{?\link{makeOrganismDbFromBiomart}} for a similar function
  that feeds data from a BioMart database.
}

\value{A \link{OrganismDb} object.}

\author{
  M. Carlson and H. Pages
}

\seealso{
  \itemize{
    \item \code{\link{makeOrganismDbFromBiomart}} for convenient ways to make a
          \link{OrganismDb} object from BioMart online resources.

    \item \code{\link[rtracklayer]{ucscGenomes}} in the \pkg{rtracklayer}
          package.

    \item \code{\link{DEFAULT_CIRC_SEQS}}.

    \item The \code{\link{supportedMiRBaseBuildValues}} function for
          listing all the possible values for the \code{miRBaseBuild}
          argument.

    \item The \link{OrganismDb} class.

  }
}

\examples{
## Display the list of genomes available at UCSC:
library(rtracklayer)
ucscGenomes()[ , "db"]

## Display the list of tables supported by makeOrganismDbFromUCSC():
supportedUCSCtables()

\dontrun{
## Retrieving a full transcript dataset for Yeast from UCSC:
txdb1 <- makeOrganismDbFromUCSC(genome="sacCer2", tablename="ensGene")
}

## Retrieving an incomplete transcript dataset for Mouse from UCSC
## (only transcripts linked to Entrez Gene ID 22290):
transcript_ids <- c(
    "uc009uzf.1",
    "uc009uzg.1",
    "uc009uzh.1",
    "uc009uzi.1",
    "uc009uzj.1"
)

txdb2 <- makeOrganismDbFromUCSC(genome="mm9", tablename="knownGene",
                          transcript_ids=transcript_ids)
txdb2
}
