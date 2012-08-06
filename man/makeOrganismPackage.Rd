\name{makeOrganismPackage}

\alias{makeOrganismPackage}


\title{
  Making AnnotationOrganismDb packages from annotation packages.
}
\description{
 \code{makeAnnotationOrganismDb} is a method that generates a package
  that will load an appropriate \code{annotationOrganismDb} object that
  will in turn point to existing annotation packages.
}
\usage{
makeOrganismPackage (pkgname,
                     OrgPkg,
                     TxDbPkg,
                     version,
                     maintainer,
                     author,
                     destDir=".",
                     license="Artistic-2.0")
}
\arguments{
  \item{pkgname}{What is the desired package name.  Traditionally, this
  should be the genus and species separated by a ".".  So as an example,
  "Homo.sapiens" would be the package name for human}
  \item{OrgPkg}{The name of an existing or generated organism package}
  \item{TxDbPkg}{The name of an existing or generated TranscriptDb package}
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
  query across multiple annotationd resources via information contained
  by the \code{annotationOrganismDb} object.
}

\value{A special package to load an \link{AnnotationOrganismDb} object.}

\author{
  M. Carlson
}

\seealso{
  \code{\link{AnnotationOrganismDb}}
}

\examples{
## set up the data.frame with the relevant relationships:
xDbs <- c("GO.db","org.Hs.eg.db")
yDbs <- c("org.Hs.eg.db","TxDb.Hsapiens.UCSC.hg19.knownGene")
xKeys <- c("GOID","ENTREZID")
yKeys <- c("GO","GENEID")
gd <- data.frame(cbind(xDbs, yDbs, xKeys, yKeys))

## makes an Organism package for human called Homo.sapiens
makeOrganismPackage(pkgname = "Homo.sapiens",
  graphData = gd,
  organism = "Homo sapiens",
  version = "1.0.0",
  maintainer = "Bioconductor Package Maintainer <maintainer@bioconductor.org>",
  author = "Bioconductor Core Team",
  destDir = ".",
  license = "Artistic-2.0")
}



% ## set up the data.frame with the relevant relationships:
% xDbs <- c("GO.db","org.Rn.eg.db")
% yDbs <- c("org.Rn.eg.db","TxDb.Rnorvegicus.UCSC.rn4.ensGene")
% xKeys <- c("GOID","ENSEMBL")
% yKeys <- c("GO","GENEID")
% gd <- data.frame(cbind(xDbs, yDbs, xKeys, yKeys))

% ## makes an Organism package for human called Homo.sapiens
% makeOrganismPackage(pkgname = "Rattus.norvegicus",
%   graphData = gd,
%   organism = "Rattus norvegicus",
%   version = "1.0.0",
%   maintainer = "Bioconductor Package Maintainer <maintainer@bioconductor.org>",
%   author = "Bioconductor Core Team",
%   destDir = ".",
%   license = "Artistic-2.0")


