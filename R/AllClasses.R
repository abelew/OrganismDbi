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
             representation(keys="matrix",
                            graph="graphNEL")
)




## A generalized constructor (in the style of loadDb).


## helpers to get all supporting libs loaded
.initPkg <- function(pkg){
  if (missing(pkg)){
    stop("'", pkg, "' is required, please install it")
  }else{
    require(pkg, character.only = TRUE)
  }
}

## helper for extracting pkgs and cols as a vector
.extractPkgsAndCols <- function(gd){
  ##  gd <- as.matrix(gd)
  setNames(as.vector(gd[,3:4]), as.vector(gd[,1:2]))
}


## Constructor 
OrganismDb <- function(dbType, graphData, ...){
  ## make graphData into a graphNEL
  ## FIXME: validate graphData -- required columns?
  gd <- graphData    
  graph <- ftM2graphNEL(gd[,1:2], edgemode="undirected")
  
  ## We should try to call require on all the supporting packages.
  pkgs <- unique(names(.extractPkgsAndCols(gd)))
  for (pkg in pkgs)
      .initPkg(pkg)

  ## Then make the object.
  new("OrganismDb", ..., keys=graphData, graph=graph)
}




###########################################################
## Convenience function that will load the package.
## IOW, there will be a call to this in zzz.R
.loadOrganismDbiPkg <- function(pkgname,
                                graphData){
  obj <- OrganismDb(pkgname, graphData)
  ns <- asNamespace(pkgname)
  assign(pkgname, obj, envir=ns)
  namespaceExport(ns, pkgname)
}




## Some getter methods to access the slots
setGeneric("keyFrame", function(x) standardGeneric("keyFrame"))
setMethod("keyFrame", "OrganismDb",
    function(x){x@keys}
)

setGeneric("dbGraph", function(x) standardGeneric("dbGraph"))
setMethod("dbGraph", "OrganismDb",
    function(x){x@graph}
)

## Then some helpers to process some of these results a bit
.getDbObjNames <- function(x){
  gd <- keyFrame(x)
  unique(c(gd[,1],gd[,2]))
}

.getDbObjs <- function(x){
  dbs <- .getDbObjNames(x)
  setNames(lapply(dbs, .makeReal), dbs)
}


## Show method (I am not really sure what to put here)
setMethod("show", "OrganismDb",
    function(object)
    {
        cat("class:", class(object), "\n")
        cat("Annotation resources:\n")
        objs <- .getDbObjNames(object)
        show(objs)
        cat("Annotation relationships:\n")
        kf <- keyFrame(object)
        show(kf)
    }
)








###############################################################################
## Rules for these kinds of packages:
###############################################################################
## 1/2) There must be a select interface implemented.

## 1) You cannot have more than one example of each field that can be retrieved from each type of package that is included.  So basically, all values for cols must be unique across ALL supporting packages.  You cannot combine resources together if the cols are not unique.  So if one package has a cols value that is "FOO", you cannot add any other packages that have a value of "FOO" for their cols.

## 2) You cannot have more than one example of each object type.  So you cannot have two org packages (for example).

## 3) You cannot have cycles in the graph.  Or maybe you can, but it is a bad idea because it can generate unpredictable results when the algorithm for walking along the tree nodes is used to interpolate nodes.  IOW, whenever the algorithm has to traverse a cycle the route it takes will be consistent, but may not be the route that the user intended).
