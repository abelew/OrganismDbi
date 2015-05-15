## New helper to lookup which org object or package should be used
## based on the taxonomy ID.  It takes a tax ID and returns an appropriate OrgDb object.
.packageTaxIds <- function(){
    c('180454'='org.Ag.eg.db',
      '3702'='org.At.tair.db',
      '9913'='org.Bt.eg.db',
      '9615'='org.Cf.eg.db',
      '9031'='org.Gg.eg.db',
      '9598'='org.Pt.eg.db',
      '511145'='org.EcK12.eg.db',
      '386585'='org.EcSakai.eg.db',
      '7227'='org.Dm.eg.db',
      '9606'='org.Hs.eg.db',
      '10090'='org.Mm.eg.db',
      '9823'='org.Ss.eg.db',
      '10116'='org.Rn.eg.db',
      '9544'='org.Mmu.eg.db',
      '6239'='org.Ce.eg.db',
      '8355'='org.Xl.eg.db',
      '559292'='org.Sc.sgd.db',
      '7955'='org.Dr.eg.db',
      '36329'='org.Pf.plasmo.db')
}
    
.taxIdToOrgDb <- function(taxid){
    ## These are the packaged TaxIds
    packageTaxIds <- .packageTaxIds() 
    if(taxid %in% names(packageTaxIds)){
        pkg <- packageTaxIds[names(packageTaxIds) %in% taxid]
        library(pkg, character.only = TRUE)
        res <- get(pkg)
    }else{
        ## If we don't have a package, then lets get the taxIds and AHIds
        ## for the hub objects
        require(AnnotationHub)
        ah <- AnnotationHub()
        ah <- subset(ah, ah$rdataclass=='OrgDb') 
        mc <- mcols(ah)[,'taxonomyid', drop=FALSE]
        ## Then just get the object
        AHID <- rownames(mc[mc$taxonomyid==taxid,,drop=FALSE])
        res <- ah[[AHID]]
    }
    res
}
## examples:
## .taxIdToOrgDb(9606)
## .taxIdToOrgDb(9986)


## Need another helper to get us from taxID to the OrgDbName...
.taxIdToOrgDbName <- function(taxid){
    packageTaxIds <- .packageTaxIds() 
    if(taxid %in% names(packageTaxIds)){
        pkg <- packageTaxIds[names(packageTaxIds) %in% taxid]
        library(pkg, character.only = TRUE)
        obj <- get(pkg)
        path <- dbfile(obj)
        pathSplit <- unlist(strsplit(path, split=.Platform$file.sep))
        res <- sub("sqlite","db", pathSplit[length(pathSplit)])
    }else{
        ## If we don't have a package, then lets get the taxIds and AHIds
        ## for the hub objects
        require(AnnotationHub)
        ah <- AnnotationHub()
        ah <- subset(ah, ah$rdataclass=='OrgDb') 
        mc <- mcols(ah)[,c('taxonomyid','title'), drop=FALSE]
        ## Then just get the object
        data <- mc[mc$taxonomyid==taxid,,drop=FALSE]
        res <- sub("sqlite","db", data$title) 
    }
    res
}
## examples
## .taxIdToOrgDbName(9606)
## .taxIdToOrgDbName(9986)


#############################################################################

## simplify DB retrievals from metadata table
.getMetaDataValue <- function(db, name){
  con <- AnnotationDbi::dbconn(db)
  res <- dbGetQuery(con,
    paste0("SELECT value FROM metadata WHERE name='", name,"'"))[[1]]
  if(!is.character(res))
      stop("missing metadata table value for: ", name)
  res
}

## function to allow us to convert a list into the inernally preferred form...
.mungeGraphData <- function(graphData){
  ##pkgs <- sapply(graphData, names) ## This is no good.   :(
  pkgs <- matrix(unlist(lapply(graphData, names)), ncol=2, byrow=TRUE)
  keys <- matrix(unlist(graphData), ncol=2, byrow=TRUE)
  graphData <- cbind(pkgs, keys)
  colnames(graphData) <- c("xDbs","yDbs","xKeys","yKeys")
  graphData
}

## early sanity checks for graphData
.testGraphData <- function(graphData){
  if(ncol(graphData) !=4L)
      stop("'graphData' must contain exactly 4 columns.")
}

## test keys for graphData 
.testKeys <- function(fkeys){
  pkgs <- unlist(lapply(names(fkeys), .makeReal))
  res <- logical(length(pkgs))
  for(i in seq_len(length(pkgs))){
    res[i] <- fkeys[i] %in% columns(pkgs[[i]]) 
  }  
  if(!all(res))
    stop("some foreign keys are not present in their associated databases")
}

## helper to list bioc Annot packages (for filling in things like
## suggests fields)
.biocAnnPackages <- function(){
    availAnns <- as.data.frame(available.packages(
                               contrib.url(biocinstallRepos()[["BioCann"]],
                               "source")))
    as.character(availAnns[["Package"]])
}


## Helper to set up to just load packages that need loading.
.extractDbFiles <- function(gd, deps){
    pkgs <- unique(names(.extractPkgsAndCols(gd)))
    ## Before we can proceed, we may need to call library on the deps...
    lapply(deps, library, character.only = TRUE)
    files <- unlist(lapply(pkgs, function(x){dbfile(get(x))}))
    setNames(files, pkgs)
}


## We want makeOrganismPackage to be self contained (have all it needs)
## IOW we want to store .sqlite files in a local inst/extdata when
## they are not known packages.
makeOrganismPackage <- function(pkgname,
                                graphData,
                                organism,
                                version,
                                maintainer,
                                author,
                                destDir,
                                license="Artistic-2.0"){
   ## there should only be one template
   template_path <- system.file("OrgPkg-template",package="OrganismDbi")
   ## We need to get a list of dependencies:
   ## 1st convert graphData into a data.frame
   gd <- .mungeGraphData(graphData)
   allDeps <- unique(as.vector(gd[,1:2]))
   ## Filter dependencies to make sure they are really package names
   biocPkgNames <- .biocAnnPackages()
   deps <- allDeps[allDeps %in% biocPkgNames]
   depsStr <- paste(deps,collapse=", ")   
   ## We need to define some symbols in order to have the
   ## template filled out correctly. 
   symvals <- list(
    PKGTITLE=paste("Annotation package for the",pkgname,"object"),
    PKGDESCRIPTION=paste("Contains the",pkgname,"object",
      "to access data from several related annotation packages."),
    PKGVERSION=version,
    AUTHOR=author,
    MAINTAINER=maintainer,
    LIC=license,        
    ORGANISM=organism,
    ORGANISMBIOCVIEW=gsub(" ","_",organism),
    PKGNAME=pkgname,
    DEPENDENCIES=depsStr
   )
   ## Check the graphData object and rename if needed
   .testGraphData(gd)
## Try to call require on all the supporting packages.
## pkgs <- unique(names(.extractPkgsAndCols(gd)))
## for (pkg in pkgs)
##     .initPkg(pkg)
   
   ## ######################################################################### 
   ## Extract the dbFile information from each object and store that
   ## into resources
   resources <- .extractDbFiles(gd, deps)
   ## ######################################################################### 
   
   ## Also check that the fkeys are really columns for the graphData
   fkeys <- .extractPkgsAndCols(gd)
   .testKeys(fkeys)
   
   ## Should never have duplicates
   if (any(duplicated(names(symvals))))
       stop("'symvals' contains duplicated symbols")
   ## All symvals should by single strings (non-NA)
   is_OK <- sapply(symvals, isSingleString)
   if (!all(is_OK)) {
       bad_syms <- paste(names(is_OK)[!is_OK], collapse="', '")
       stop("values for symbols '", bad_syms, "' are not single strings")
   }
   createPackage(pkgname=pkgname,
                 destinationDir=destDir,
                 originDir=template_path,
                 symbolValues=symvals)

   ## Now just do work to save the data.frame (pass that in instead of the
   ## other stuff) in /data as a serialized R file.
   ## There will already be a /data dir in the template
   ## So just save to it:
   graphData <- gd
   graphInfo <- list(graphData=graphData, resources=resources)
   ## create data dir (because R CMD build removes empty dirs)
   ## And then save the data there.
   dir.create(file.path(destDir,pkgname,"data"))   
   save(graphInfo, file=file.path(destDir,pkgname,"data","graphInfo.rda"))

   ## Get and other things that need to be saved and stash them into
   ## /inst/extdata
   otherDeps <- allDeps[!allDeps %in% biocPkgNames]
   .saveFromStr <- function(x, file){saveDb(x=get(x), file=file)}
   if(length(otherDeps)>0){
       ## Then we have to save stuff locally
       extPath <- file.path(destDir,pkgname,"inst","extdata")
       dir.create(extPath, recursive=TRUE)
       mapply(.saveFromStr, x=otherDeps, file=file.path(extPath,
                                           paste0(otherDeps,".sqlite")))
   }
}






################################################################################
## Now for some create functions that are more specialized:
## To create OrgansimDb objects from UCSC or from biomaRt
## the initial versions of these will just create the object (and not
## a package)

## Also: there are some pre-agreed upon expectations for this
## function.  It will make a specific kind of OrganismDb object.  One
## that contains certain specific elements (GO, OrgDb and TxDb - to
## start).  And the expectation is that we will always have those
## things when this function finishes.  So the contract is less
## general than before.


########
## Helper to set up to just load packages that need loading.  BUT this
## version of this function is less agressive and doesn't try to load
## a file if it already exists.  This is a different behavior than we
## want for packaging where things should be more strict.
.gentlyExtractDbFiles <- function(gd, deps){
    pkgs <- unique(names(.extractPkgsAndCols(gd)))
    ## Before we can proceed, we may need to call library on the deps...
    .library <- function(dep){
        if(!exists(dep)){
            library(dep, character.only = TRUE)
        }
    }
    lapply(deps, .library)
    files <- unlist(lapply(pkgs, function(x){dbfile(get(x))}))
    setNames(files, pkgs)
}


## from TxDb
## There are some issues with this as it's currently implemented:
## Because its not using the txdb that is passed in but is instead
## making one up and then assigning it to the global namespace.  This
## is what we want for makeOrgansimDbFromXXX when XXX is UCSC or
## BiomaRt, but its *not* what we want for a simple exsting txdb.
## It's straightforward to do something different when the txdb
## exists().  But: how to I look up it's original name in the
## .GlobalEnv ???  I need to know it's name for the graphData list
## below.
## The other big problem is that if my TxDb objects don't have a
## dbfile() then they can't be saved and re-loaded later.  But this
## problem is not "new" for this object.
## Still TODO: 1) either put the assigned functions in a special env that is accesible to these objects OR else make a subclass that can hold those objects.
## 2) make a save method that complains if any of the objects has not had saveDb called on it before calling a constructor.  And also add a warning to the constructor when somethings dbfile() is an empty string.

### Private environment for storing TxDb objects as needed (failed strategy)
## To work I would have to make the env public.  - It was never a
## great idea fo r this use case.
## TxDbObjs <- new.env(hash=TRUE, parent=emptyenv())
## I may still want to go with the other option of stashing this data
## into a subclass...  For one thing, that option can't have name
## clashes...
## Also: the name this local TxDb gets assigned to cannot be the same as is used by a package.  Otherwise a shortened 'custom' TxDb can be overwritten by a name clash with a package name...  This could end up being true even if I store the TxDb locally inside of a named sub-class.
## Also also: the name should not be made 'special' in the case where makeOrganismDbFromTxDb is called as a helper function from within makeOrganismDbFromUCSC or makeOrganismDbFromBiomart.


makeOrganismDbFromTxDb <- function(txdb){
    ## Then assign that object value to the appropriate name:
    txdbName <- GenomicFeatures:::.makePackageName(txdb)
    ## assign to global scope (b/c you need it there if you 'generated' it)
    ## Is there a better way?
    ## assign(txdbName, txdb, envir=TxDbObjs) #.GlobalEnv)  
    assign(txdbName, txdb, .GlobalEnv)  
    ## Then get the tax ID:
    taxId <- taxonomyId(txdb)
    
    ## Then get the name and valued for the OrgDb object
    orgdbName <- OrganismDbi:::.taxIdToOrgDbName(taxId)
    orgdb <- OrganismDbi:::.taxIdToOrgDb(taxId)
    assign(orgdbName, orgdb)
    ## get the primary key for the OrgDb object:
    geneKeyType <- AnnotationDbi:::.chooseCentralOrgPkgSymbol(orgdb)
    
    graphData <- list(join1 = setNames(object=c('GOID', 'GO'),
                                       nm=c('GO.db', orgdbName)),
                      join2 = setNames(object=c(geneKeyType, 'GENEID'),
                                       nm=c(orgdbName, txdbName)))
    
    ## get the organism
    organism <- organism(txdb)

    #############################################################
    ## Process and then test the graph Data
    gd <- OrganismDbi:::.mungeGraphData(graphData)
    OrganismDbi:::.testGraphData(gd)    
    allDeps <- unique(as.vector(gd[,1:2]))
    biocPkgNames <- OrganismDbi:::.biocAnnPackages()
    deps <- allDeps[allDeps %in% biocPkgNames]
    resources <- OrganismDbi:::.gentlyExtractDbFiles(gd, deps)    
    ## Check that the fkeys are really columns for the graphData
    fkeys <- OrganismDbi:::.extractPkgsAndCols(gd)
    OrganismDbi:::.testKeys(fkeys)
    ## Then make the object:
    graphInfo <- list(graphData=gd, resources=resources)
    OrganismDbi:::OrganismDb(graphInfo=graphInfo)
}


## from UCSC
makeOrganismDbFromUCSC <- function(genome="hg19",
                                   tablename="knownGene",
                                   transcript_ids=NULL,
                                   circ_seqs=DEFAULT_CIRC_SEQS,
                                   url="http://genome.ucsc.edu/cgi-bin/",
                     goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath",
                                   miRBaseBuild=NA){

    ## So call the function to make that TxDb
    txdb <- makeTxDbFromUCSC(genome=genome,
                             tablename=tablename,
                             transcript_ids=transcript_ids,
                             circ_seqs=circ_seqs,
                             url=url,
                             goldenPath_url=goldenPath_url,
                             miRBaseBuild=miRBaseBuild)
    makeOrganismDbFromTxDb(txdb)
    ## ## Then assign that object value to the appropriate name:
    ## txdbName <- GenomicFeatures:::.makePackageName(txdb)
    ## assign(txdbName, txdb) ## txdbName still represents the correct thing
    ## ## Then get the tax ID:
    ## taxId <- taxonomyId(txdb)
    
    ## ## Then get the name and valued for the OrgDb object
    ## orgdbName <- OrganismDbi:::.taxIdToOrgDbName(taxId)
    ## orgdb <- OrganismDbi:::.taxIdToOrgDb(taxId)
    ## assign(orgdbName, orgdb)
    ## ## get the primary key for the OrgDb object:
    ## geneKeyType <- AnnotationDbi:::.chooseCentralOrgPkgSymbol(orgdb)
    
    ## graphData <- list(join1 = setNames(object=c('GOID', 'GO'),
    ##                                    nm=c('GO.db', orgdbName)),
    ##                   join2 = setNames(object=c(geneKeyType, 'GENEID'),
    ##                                    nm=c(orgdbName, txdbName)))
    
    ## ## get the organism
    ## organism <- organism(txdb)

    ## #######################################################################3
    ## ## Process and then test the graph Data
    ## gd <- OrganismDbi:::.mungeGraphData(graphData)
    ## OrganismDbi:::.testGraphData(gd)    
    ## allDeps <- unique(as.vector(gd[,1:2]))
    ## biocPkgNames <- OrganismDbi:::.biocAnnPackages()
    ## deps <- allDeps[allDeps %in% biocPkgNames]
    ## resources <- OrganismDbi:::.gentlyExtractDbFiles(gd, deps)
    ## ## Check that the fkeys are really columns for the graphData
    ## fkeys <- OrganismDbi:::.extractPkgsAndCols(gd)
    ## OrganismDbi:::.testKeys(fkeys)
    ## ## Then make the object:
    ## graphInfo <- list(graphData=gd, resources=resources)
    ## OrganismDbi:::OrganismDb(graphInfo=graphInfo)
}

## ## Usage/testing:
## library(OrganismDbi)
## ODb <- OrganismDbi:::makeOrganismDbFromUCSC(genome="hg19",
##                                    tablename="knownGene",
##                                    transcript_ids=NULL,
##                                    circ_seqs=DEFAULT_CIRC_SEQS,
##                                    url="http://genome.ucsc.edu/cgi-bin/",
##                      goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath",
##                                    miRBaseBuild=NA)


## TODO: Document this.
## TODO: remove all the uncessary OrganismDbi::: 's


## TODO: make one of these for BiomaRt.

makeOrganismDbFromBiomart <- function(biomart="ensembl",
                                      dataset="hsapiens_gene_ensembl",
                                      transcript_ids=NULL,
                                      circ_seqs=DEFAULT_CIRC_SEQS,
                                      filters="",
                                      id_prefix="ensembl_",
                                      host="www.biomart.org",
                                      port=80,
                                      miRBaseBuild=NA){

    ## So call the function to make that TxDb
    txdb <- makeTxDbFromBiomart(biomart=biomart,
                                dataset=dataset,
                                transcript_ids=transcript_ids,
                                circ_seqs=circ_seqs,
                                filters=filters,
                                id_prefix=id_prefix,
                                host=host,
                                port=port,
                                miRBaseBuild=miRBaseBuild)
    makeOrganismDbFromTxDb(txdb)
    ## ## Then assign that object value to the appropriate name:
    ## txdbName <- GenomicFeatures:::.makePackageName(txdb)
    ## ## assign to global scope (b/c you need it there if you 'generated' it)
    ## ## Is there a better way?
    ## assign(txdbName, txdb,envir = .GlobalEnv)  
    ## ## Then get the tax ID:
    ## taxId <- taxonomyId(txdb)
    
    ## ## Then get the name and valued for the OrgDb object
    ## orgdbName <- OrganismDbi:::.taxIdToOrgDbName(taxId)
    ## orgdb <- OrganismDbi:::.taxIdToOrgDb(taxId)
    ## assign(orgdbName, orgdb)
    ## ## get the primary key for the OrgDb object:
    ## geneKeyType <- AnnotationDbi:::.chooseCentralOrgPkgSymbol(orgdb)
    
    ## graphData <- list(join1 = setNames(object=c('GOID', 'GO'),
    ##                                    nm=c('GO.db', orgdbName)),
    ##                   join2 = setNames(object=c(geneKeyType, 'GENEID'),
    ##                                    nm=c(orgdbName, txdbName)))
    
    ## ## get the organism
    ## organism <- organism(txdb)

    ## #######################################################################
    ## ## Process and then test the graph Data
    ## gd <- OrganismDbi:::.mungeGraphData(graphData)
    ## OrganismDbi:::.testGraphData(gd)    
    ## allDeps <- unique(as.vector(gd[,1:2]))
    ## biocPkgNames <- OrganismDbi:::.biocAnnPackages()
    ## deps <- allDeps[allDeps %in% biocPkgNames]
    ## resources <- OrganismDbi:::.gentlyExtractDbFiles(gd, deps)    
    ## ## Check that the fkeys are really columns for the graphData
    ## fkeys <- OrganismDbi:::.extractPkgsAndCols(gd)
    ## OrganismDbi:::.testKeys(fkeys)    
    ## ## Then make the object:
    ## graphInfo <- list(graphData=gd, resources=resources)
    ## OrganismDbi:::OrganismDb(graphInfo=graphInfo)
}



## ## Usage/testing:
## library(OrganismDbi)
## transcript_ids <- c(
##     "ENST00000013894",
##     "ENST00000268655",
##     "ENST00000313243",
##     "ENST00000435657",
##     "ENST00000384428",
##     "ENST00000478783"
## )
## ODb <- OrganismDbi:::makeOrganismDbFromBiomart(biomart="ensembl",
##                                             dataset="hsapiens_gene_ensembl",
##                                             transcript_ids=transcript_ids,
##                                             circ_seqs=DEFAULT_CIRC_SEQS,
##                                             filters="",
##                                             id_prefix="ensembl_",
##                                             host="www.biomart.org",
##                                             port=80,
##                                             miRBaseBuild=NA)

## TxDb.Hsapiens.BioMart.ensembl.GRCh38.p2


## PROBLEM: OrganismDbi:::.extractDbFiles(gd, deps) requires (strictly) that all objects be available as files somewhere (no exceptions allowed)
## This means that when I get to this stage, with biomaRt, it fails because there is not a TxDb on disc...


