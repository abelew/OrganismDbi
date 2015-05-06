## New helper to lookup which org object or package should be used
## based on the taxonomy ID.  It takes a tax ID and returns an appropriate OrgDb object.
.taxIdToOrgDb <- function(taxid){
    ## These are the packaged TaxIds
    packageTaxIds <- c('180454'='org.Ag.eg.db',
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



