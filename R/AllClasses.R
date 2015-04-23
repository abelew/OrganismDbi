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
## slot for a GODb, then lets add in a slot for a TxDb.  Since
## AnnotationDbi does not know about TxDbs, I will have to define the
## base class here (with the aim of creating a new software package to hold
## this stuff later on.


## require(AnnotationDbi)
## require(GenomicFeatures)




## I need an initialize method just to allow me to do things like
## require(GO.db) etc.

## Original class
OrganismDb <-
    setClass("OrganismDb",
             representation(keys="matrix",
                            graph="graphNEL")
)

## ## Reference style class that does not work
## ## Also this recapitulates an odd bug that keeps popping up...
## OrganismDb <-
##     setRefClass("OrganismDb",
##              contains="TxDb",
##              fields=list(keys="matrix", graph="graphNEL")
## )

## ## Reference style class that does work (but no inheritance)
## OrganismDb <-
##     setRefClass("OrganismDb",
##           ##   contains="TxDb",
##              fields=list(keys="matrix", graph="graphNEL")
## )


 
.cantFindResourceMsg <- function(pkg){
    (paste0("The '", pkg, "' object is required but not ",
            "found. Please either install it (if it's from ",
            "a package) or ensure that it's loaded to the ",
            "search path."))
}

## helpers to get all supporting libs loaded
## .initPkg needs to:
## 1) see if the package is on the search path?
## 2) if not on search path, try to see if it's installed (and load if needed).
## 3) emit an appropriate warning in either case. 
## NOTE:
## This function deliberately does not use my .biocAnnPackages()
## function because I have no way of knowing if someone else
## has made a custom one or is using one from another repos. etc.
.initPkg <- function(pkg, OrganismDbPkgName, ns=NULL){
##    message("pkg is:", pkg)
##    message("OrganismDbPkgName is:", OrganismDbPkgName)
    
    if(!exists(pkg)){ ## IOW there is no object
        if(suppressWarnings(!library(pkg,
                                     character.only = TRUE,
                                     logical.return=TRUE))){
            if(!is.null(OrganismDbPkgName)){
##    message("The '", pkg, "' pkg is now trying to load from 'inst/extdata/'.")
                pkgPath <- system.file("extdata", paste0(pkg,".sqlite"),
                           package=OrganismDbPkgName)
##    message("SEARCHING this path:", pkgPath)
                msg <- .cantFindResourceMsg(pkg)
                tryCatch(loadDb(pkgPath),
                         error = function(e){stop(wmsg(msg))} )
            }
        }
    }
}

## helper for extracting pkgs and cols as a vector
.extractPkgsAndCols <- function(gd){
  ##  gd <- as.matrix(gd)
  setNames(as.vector(gd[,3:4]), as.vector(gd[,1:2]))
}


## Constructor 
OrganismDb <- function(dbType=NULL, graphData, ns=NULL, ...){
  ## make graphData into a graphNEL
  ## FIXME: validate graphData -- required columns?
  gd <- graphData    
  graph <- ftM2graphNEL(gd[,1:2], edgemode="undirected")
  
  ## We should try to call require on all the supporting packages.
  pkgs <- unique(names(.extractPkgsAndCols(gd)))
  for (pkg in pkgs)
      .initPkg(pkg, dbType, ns=ns)

  ## Then make the object.
  new("OrganismDb", ..., keys=graphData, graph=graph)
}


## helper that is just used for those resources that are not separate
## packages, but which need to have loadDb called (and be sealed into
## the namespace for the OrganismDb object)
.addLocalPkgsToNamespace <- function(pkgname, graphData, ns){
    pkgs <- unique(names(.extractPkgsAndCols(graphData)))
    xsts <- sapply(pkgs, exists)
    pkgs <- pkgs[!xsts]
##    pkgVals <- list()
    ## then get the ones that we don't already have and seal to the namespace.
    for(pkg in pkgs){
        msg <- .cantFindResourceMsg(pkg)
        pkgPath <- system.file("extdata", paste0(pkg,".sqlite"),
                               package=pkgname)
        tryCatch({assign(eval(pkg), loadDb(pkgPath))},
                  error = function(e){stop(wmsg(msg))} )
        ## save what we are getting
##        pkgVals <- c(pkgVals, get(eval(pkg)))
        ## and assign it to the namespace
        assign(pkg, get(eval(pkg)), envir=ns)
        namespaceExport(ns, pkg)
    }
    ## then return the namespace...
##    return(list(ns=ns, vars=pkgVals))
}

###########################################################
## Convenience function that will load the package.
## IOW, there will be a call to this in zzz.R
.loadOrganismDbiPkg <- function(pkgname,
                                graphData){
  ns <- asNamespace(pkgname)
  ## trouble: you can only do namespaceExport once... So I need a
  ## separate helper (not the constructor) to add any local objects
  ## into the namespace before we seal it (below).
  .addLocalPkgsToNamespace(pkgname, graphData, ns)
  ## THEN make the object...
  obj <- OrganismDb(pkgname, graphData, ns)
  assign(pkgname, obj, envir=ns)
  ## now export it
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
        cat("For more details, please see the show methods for the component objects listed above.\n")
    }
)

setMethod("metadata", "OrganismDb",
    function(x)
    {
        objs <- .getDbObjs(x)
        res <- data.frame()
        for(i in seq_len(length(objs))){
            meta <- metadata(objs[[i]])
            res <- rbind(meta,res) ## inneficient (but i is normally <= 4)
        }
        ## Then get rid of any rows that are repeated
        d <- duplicated(res[,1])
        dupNames <- unique(res[d,1])
        res <- res[!(res[,1] %in% dupNames),]
        message("Only unique values are returned by this metadata method.  For individual metadata values that may share a key, such as the Db type, be sure to call metadata on the individual objects. \n")
        show(res)
    }
)         






###############################################################################
## Rules for these kinds of packages:
###############################################################################
## 1/2) There must be a select interface implemented.

## 1) You cannot have more than one example of each field that can be retrieved from each type of package that is included.  So basically, all values for cols must be unique across ALL supporting packages.  You cannot combine resources together if the cols are not unique.  So if one package has a cols value that is "FOO", you cannot add any other packages that have a value of "FOO" for their cols.

## 2) You cannot have more than one example of each object type.  So you cannot have two org packages (for example).

## 3) You cannot have cycles in the graph.  Or maybe you can, but it is a bad idea because it can generate unpredictable results when the algorithm for walking along the tree nodes is used to interpolate nodes.  IOW, whenever the algorithm has to traverse a cycle the route it takes will be consistent, but may not be the route that the user intended).
