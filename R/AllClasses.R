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



require(AnnotationDbi)
require(GenomicFeatures)


## I need an initialize method just to allow me to do things like
## require(GO.db) etc.
AnnotationOrganismDb <-
    setClass("AnnotationOrganismDb",
        representation(GODb="GODb", OrgDb="OrgDb", TranscriptDb="TranscriptDb")
)




.loadOrganismDbiPkg <- function(pkgname,
                                GOPkg,
                                OrgPkg,
                                TranscriptPkg){
  obj <- OrganismDbi:::AnnotationOrganismDb(dbType= pkgname,
                                           GOPkg=GOPkg,
                                           OrgPkg=OrgPkg,
                                           TranscriptPkg=TranscriptPkg)
  ns <- asNamespace(pkgname)
  assign(pkgname, obj, envir=ns)
  namespaceExport(ns, pkgname)
}




## A generalized constructor (in the style of loadDb), except notice that
## there can be ONLY one type of DB per package.  This constraint is not
## shared by loadDb, where one package can just be like an instance of a
## particular type.  Here we only want one Homo.sapiens package.  As currently
## imagined, users can toggle the contents of the slots if they want to change
## it.

## helpers
.initPkg <- function(pkg){
  if (missing(pkg)){
    stop(paste("The",pkg,"argument is strictly required, please supply it"))
  }else{
    require(pkg, character.only = TRUE)
    res <- eval(parse(text=pkg))
  }
  return(res)
}

## generalize this to loop over a vector of DBs ...
.initAnnotationDbs <- function(GOPkg, OrgPkg, TranscriptPkg){
  .GODb <- .initPkg(GOPkg)
  .OrgDb <- .initPkg(OrgPkg)
  .TranscriptDb <- .initPkg(TranscriptPkg)
  list(GODb=.GODb, OrgDb=.OrgDb, TranscriptDb=.TranscriptDb)
}

## Constructor helper (guts) OR 
## function version of constructor 
AnnotationOrganismDb <- function(dbType, GOPkg, OrgPkg, TranscriptPkg){
    Dbs <- .initAnnotationDbs(GOPkg, OrgPkg, TranscriptPkg)
    #obj <- getRefClass(dbType, where=getNamespace(dbType))
    obj <- getClassDef(dbType, where=getNamespace(dbType))
    #obj$new(GODb=Dbs[['GODb']],OrgDb=Dbs[['OrgDb']],
    #        TranscriptDb=Dbs[['TranscriptDb']])
    new(obj,GODb=Dbs[['GODb']],OrgDb=Dbs[['OrgDb']],
        TranscriptDb=Dbs[['TranscriptDb']])
}

## Turns out that I just don't need a method here at all.  I will just pass
## the string I need for the unique organism name and (providing that there is
## a type for that defined above (or in someone elses code), then it should
## work fine with a generic constructor function.

## Constructor Generic
## setGeneric("AnnotationOrganismDb",
##            function(dbType, GOPkg, OrgPkg, TranscriptPkg)
##            standardGeneric("AnnotationOrganismDb"))


## Constructor for method (so that later I can call
## AnnotationOrganismDb("Mus.musculus" ...) etc.
## setMethod(AnnotationOrganismDb,
##           c("character", "character", "character", "character"),
##           .AnnotationOrganismDb)




## planned usage.
## Works
##  hs <- AnnotationOrganismDb(dbType= "Homo.sapiens", GOPkg="GO.db", OrgPkg="org.Hs.eg.db", TranscriptPkg="TxDb.Hsapiens.UCSC.hg19.knownGene")



