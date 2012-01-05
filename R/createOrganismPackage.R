## simplify DB retrievals from metadata table
.getMetaDataValue <- function(db, name){
  con <- AnnotationDbi:::dbConn(db)
  res <- dbGetQuery(con,
    paste("SELECT value FROM metadata WHERE name='",
      name,"'", sep=""))[[1]]
  if(!is.character(res))error("Your metadata table is missing a value for:",
    name,".")
  res
}

makeOrganismPackage <- function(pkgname,
                                OrgPkg,
                                TxDbPkg,
                                version,
                                maintainer,
                                author,
                                destDir=".",
                                license="Artistic-2.0"){
   ## there should only be one template
   template_path <- system.file("OrgPkg-template",package="OrganismDbi")
   ## We need to define some symbols in order to have the
   ## template filled out correctly.
   require(OrgPkg, character.only=TRUE)
   db <- eval(parse(text=OrgPkg))
   symvals <- list(
    ORGPKG = OrgPkg,
    TXDBPKG = TxDbPkg,
    PKGTITLE=paste("Annotation package for the",pkgname,
      "object"),
    PKGDESCRIPTION=paste("Contains the",pkgname,"object",
      "which allows access to data from the GO.db,",OrgPkg,"and",TxDbPkg,
      "packages."),
    PKGVERSION=version,
    AUTHOR=author,
    MAINTAINER=maintainer,
    LIC=license,        
    ORGANISM=.getMetaDataValue(db,'ORGANISM'),
    ORGANISMBIOCVIEW=gsub(" ","_",.getMetaDataValue(db,'ORGANISM')),
    PKGNAME=pkgname 
   )
   ## Should never have duplicates
   if (any(duplicated(names(symvals)))) {
       str(symvals)
       stop("'symvals' contains duplicated symbols")
   }
   ## All symvals should by single strings (non-NA)
   is_OK <- sapply(symvals, isSingleString)
   if (!all(is_OK)) {
       bad_syms <- paste(names(is_OK)[!is_OK], collapse=", ")
       stop("values for symbols ", bad_syms, " are not single strings")
   }
   createPackage(pkgname=pkgname,
                 destinationDir=destDir,
                 originDir=template_path,
                 symbolValues=symvals)

}
