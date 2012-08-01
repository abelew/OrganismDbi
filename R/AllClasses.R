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


## helpers
## TODO: modify this so that it inits all the members of graph

.initPkg <- function(pkg){
  if (missing(pkg)){
    stop(paste(pkg,"is strictly required, please make sure you have installed it"))
  }else{
    require(pkg, character.only = TRUE)
  }
}



## Constructor helper (guts) OR 
## function version of constructor 
OrganismDb <- function(dbType, graphData){
    ## make graphData into a graphNEL
    gd <- as.matrix(graphData)
    graph <- ftM2graphNEL(gd[,1:2], edgemode="undirected")
    
    ## We should try to call require on all the supporting packages.
    pkgs <- unique(c(gd[,1],gd[,2]))
    lapply(pkgs, .initPkg)
    ## Then make the object.

    ## TODO: getClassDef() presumes that a package exists that defines
    ## this already (so you can't just make one on the fly this way)
    ## If possible I think I would like to make it so that I can just generate
    ## an object here without requiring that a package exist ahead of time...
    obj <- getClassDef(dbType, where=getNamespace(dbType))
    
    ## BUT I cannot (for example) do this:
    ## obj <- setClass(dbType, contains="OrganismDb")
    ## assign(dbType, obj)
    ## Because the environment is sealed once the package is loaded, so by the
    ## time I am defining this, it's too late...  The only way around this is
    ## probably to use reference classes, but those had other problems.  For
    ## now, this is proably sufficient.
    new(obj, keys=graphData, graph=graph)
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


