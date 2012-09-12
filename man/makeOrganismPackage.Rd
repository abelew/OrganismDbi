\name{makeOrganismPackage}

\alias{makeOrganismPackage}


\title{
  Making OrganismDb packages from annotation packages.
}
\description{
 \code{makeOrganismPackage} is a method that generates a package
  that will load an appropriate \code{annotationOrganismDb} object that
  will in turn point to existing annotation packages.
}
\usage{
makeOrganismPackage (pkgname,
                     graphData,
                     organism,
                     version,
                     maintainer,
                     author,
                     destDir,
                     license="Artistic-2.0")
}
\arguments{
  \item{pkgname}{What is the desired package name.  Traditionally, this
  should be the genus and species separated by a ".".  So as an example,
  "Homo.sapiens" would be the package name for human}
\item{graphData}{A list of short character vectors.  Each character
  vector in the list is exactly two elements long and represents a join
  relationship between two packages.  The names of these character
  vectors are the package names and the values are the foreign keys that
  should be used to connect each package.  All foreign keys must be
  values that can be returned by the cols method for each package in
  question, and obviously they also must be the same kind of identifier
  as well.}
  \item{organism}{The name of the organism this package represents}
  \item{version}{What is the version number for this package?}
  \item{maintainer}{Who is the package maintainer? (must include email to be 
    valid)}
  \item{author}{Who is the creator of this package?}
  \item{destDir}{A path where the package source should be assembled.}
  \item{license}{What is the license (and it's version)}
}
\details{
  The purpose of this method is to create a special package that will
  depend on existing annotation packages and which will load a special
  \code{annotationOrganismDb} object that will allow proper dispatch of
  special select methods.  These methods will allow the user to easily
  query across multiple annotation resources via information contained
  by the \code{annotationOrganismDb} object. Because the end result will
  be a package that treats all the data mapped together as a single
  source, the user is encouraged to take extra care to ensure that the
  different packages used are from the same build etc.
}

\value{A special package to load an \link{OrganismDb} object.}

\author{
  M. Carlson
}

\seealso{
  \code{\link{OrganismDb}}
}

\examples{
## set up the list with the relevant relationships:
gd <- list(join1 = c(GO.db="GOID", org.Hs.eg.db="GO"),
           join2 = c(org.Hs.eg.db="ENTREZID",
                     TxDb.Hsapiens.UCSC.hg19.knownGene="GENEID"))

## sets up a temporary directory for this example
## (users won't need to do this step)
destination <- tempfile()
dir.create(destination)

## makes an Organism package for human called Homo.sapiens
makeOrganismPackage(pkgname = "Homo.sapiens",
  graphData = gd,
  organism = "Homo sapiens",
  version = "1.0.0",
  maintainer = "Bioconductor Package Maintainer <maintainer@bioconductor.org>",
  author = "Bioconductor Core Team",
  destDir = destination,
  license = "Artistic-2.0")
}





