## simplify DB retrievals from metadata table
.getMetaDataValue <- function(db, name){
  con <- AnnotationDbi:::dbConn(db)
  res <- dbGetQuery(con,
    paste0("SELECT value FROM metadata WHERE name='", name,"'"))[[1]]
  if(!is.character(res))
      stop("missing metadata table value for: ", name)
  res
}

## early sanity checks for graphData
.testGraphData <- function(graphData){
  if(ncol(graphData) !=4L)
      stop("'graphData' must contain exactly 4 columns.")
  ## enforce colnames of graphData to always be uniform.
  colnames(graphData) <- c("xDbs","yDbs","xKeys","yKeys")
  graphData
}

## test keys for graphData 
.testKeys <- function(fkeys){
  pkgs <- unlist(lapply(names(fkeys), .makeReal))
  res <- logical(length(pkgs))
  for(i in seq_len(length(pkgs))){
    res[i] <- fkeys[i] %in% cols(pkgs[[i]]) 
  }  
  if(!all(res))
    stop("some foreign keys are not present in their associated databases")
}


makeOrganismPackage <- function(pkgname,
                                graphData,
                                organism,
                                version,
                                maintainer,
                                author,
                                destDir=".",
                                license="Artistic-2.0"){
   ## there should only be one template
   template_path <- system.file("OrgPkg-template",package="OrganismDbi")
   ## We need to get a list of dependencies:
   gd <- as.matrix(graphData)
   deps <- paste(unique(as.vector(gd[,1:2])),collapse=", ")
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
    DEPENDENCIES=deps
   )
   ## Check the graphData object and rename if needed
   graphData <- .testGraphData(graphData)
   ## Try to call require on all the supporting packages.
   pkgs <- unique(names(.extractPkgsAndCols(gd)))
   for (pkg in pkgs)
       .initPkg(pkg)
   ## Also check that the fkeys are really cols for the graphData
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
   save(graphData, file=file.path(destDir,pkgname,"data","graphData.Rda"))
}


