## general functions go here:

## Still need to define my object in AnnotationDbi


###############################################################################
## OrganismDb class Objects
## What is the big idea?
## These objects will contain a list of resources along with theprimary
## (package name and hence the name of the object that gets made when that
## package loads) and foreign keys for each resource.  The foreign key is to
## be listed as a keytype.  The basic notion is that these objects will allow
## us to write select() methods etc. that join contents of multiple databases.
## The relevant select ,cols, keys and keytypes methods should all call the
## base methods for each of the relevant packages, and then cat them together.
## Results should come back to the users and the object itself ought to be
## able to call select on each of the tables, and then merge it to the rest
## based on the foreign keytype etc.


## So lets start with an object where we have a slot for an OrgDb and also a
## slot for a GODb, then lets add in a slot for a TranscriptDb.  Since
## AnnotationDbi does not know about TranscriptDbs, I will have to define the
## base class here (with the aim of creating a new software package to hold
## this stuff later on.


## require(AnnotationDbi)
## require(GenomicFeatures)




## I need an initialize method just to allow me to do things like
## require(GO.db) etc.
OrganismDb <-
    setClass("OrganismDb",
             representation(keys="data.frame",
                            graph="graphNEL")
)




## A generalized constructor (in the style of loadDb).

## helper to convert text strings into objects
.makeReal <- function(x){
  eval(parse(text=x))
}

## helpers to get all supporting libs loaded
.initPkg <- function(pkg){
  if (missing(pkg)){
    stop(paste(pkg,"is strictly required, please make sure you have installed it"))
  }else{
    require(pkg, character.only = TRUE)
  }
}

## early sanity checks for graphData
.testGraphData <- function(graphData){
  if(dim(graphData)[2] !=4){stop("graphData must contain exactly 4 columns.")}
  ## enforce colnames of graphData to always be uniform.
  colnames(graphData) <- c("xDbs","yDbs","xKeys","yKeys")
  graphData
}

## helper for extracting pkgs and cols as a vector
.extractPkgsAndCols <- function(gd){
  gd <- as.matrix(gd)
  res <- c(gd[,3],gd[,4])
  pkgs <-  c(gd[,1],gd[,2])
  names(res) <- pkgs 
  res
}

## test keys for graphData (do this after making sure that pkgs are present)
.testKeys <- function(fkeys){
  pkgs <- unlist(lapply(names(fkeys), .makeReal))
  res <- logical(length(pkgs))
  for(i in seq_len(length(pkgs))){
    res[i] <- fkeys[i] %in% cols(pkgs[[i]]) 
  }  
  if(!all(res)){
    stop("some of the foreign keys supplied are not present in their associated databases.")
  }
}


## Constructor 
OrganismDb <- function(dbType, graphData){
  ## Check the graphData object
  graphData <- .testGraphData(graphData)
    
  ## make graphData into a graphNEL
  gd <- as.matrix(graphData)    
  graph <- ftM2graphNEL(gd[,1:2], edgemode="undirected")
  
  ## We should try to call require on all the supporting packages.
  pkgs <- unique(names(.extractPkgsAndCols(gd)))
  lapply(pkgs, .initPkg)

  ## Check that the fkeys are really cols for the graphData
  fkeys <- .extractPkgsAndCols(gd)
  .testKeys(fkeys)
  
  ## Then make the object.
  new("OrganismDb", keys=graphData, graph=graph)
}



## planned usage.  I think the constructor should take a data.frame.  The user
## just doesn't ever need to see that.  But internally I do want to store a
## graphNEL, because I want to be able to call shortest path algorithms later
## on...

## HERE we define an example of what that data.frame might look like (one row per edge):

##  xDbs <- c("GO.db","org.Hs.eg.db")
##  yDbs <- c("org.Hs.eg.db","TxDb.Hsapiens.UCSC.hg19.knownGene")
##  xKeys <- c("GOID","ENTREZID")
##  yKeys <- c("GO","GENEID")
##  gd <- data.frame(cbind(xDbs, yDbs, xKeys, yKeys))

## Constructor should look like:
##  hs <- OrganismDbi:::OrganismDb(dbType= "Homo.sapiens", graphData=gd)




## How about this: (OR: Just how general IS this???)

##  xDbs <- c("GO.db","org.Rn.eg.db")
##  yDbs <- c("org.Hs.eg.db","TxDb.Rnorvegicus.UCSC.rn4.ensGene")
##  xKeys <- c("GOID","ENSEMBL")
##  yKeys <- c("GO","GENEID")
##  gd <- data.frame(cbind(xDbs, yDbs, xKeys, yKeys))

## Constructor should look like:
##  rn <- OrganismDbi:::OrganismDb(dbType= "Rattus.norvegicus", graphData=gd)



 

## Other TODOs:
## 1) add checks to make sure that the data.frame entered is reasonable.






###########################################################
## Convenience function that will load the package.
## IOW, there will be a call to this in zzz.R
.loadOrganismDbiPkg <- function(pkgname,
                                graphData){
  obj <- OrganismDbi:::OrganismDb(dbType= pkgname,
                                  graphData=graphData)
  ns <- asNamespace(pkgname)
  assign(pkgname, obj, envir=ns)
  namespaceExport(ns, pkgname)
}







### TODO:
## 1) Fix the show method (here or in select.R, call .getDbObjs and show those.
## 2) 






###############################################################################
## Rules for these kinds of packages:
###############################################################################


## 1) You cannot have more than one example of each field that can be retrieved from each type of package that is included.  So basically, all values for cols must be unique across ALL supporting packages.  You cannot combine resources together if the cols are not unique.  So if one package has a cols value that is "FOO", you cannot add any other packages that have a value of "FOO" for their cols.

## 2) You cannot have more than one example of each object type.  So you cannot have two org packages (for example).

## 3) You cannot have cycles in the graph.  Or maybe you can, but it is a bad idea because it can generate unpredictable results when the algorithm for walking along the tree nodes is used to interpolate nodes.  IOW, whenever the algorithm has to traverse a cycle the route it takes will be consistent, but may not be the route that the user intended).
