## Optionally:
## If I change to reference classes, then I could just create this object on
## the fly (nice bonus).
## The constructor already takes a data.frame (it just can't spawn itself) 

## For now: this step below can only happen here (because of using standard
## classes)
assign("@PKGNAME@", setClass("@PKGNAME@", contains="OrganismDb"))

## .onLoad gets the data.frame from the /data directory
.onLoad <- function(libname, pkgname) {
  load(system.file("data","graphData.Rda",package=pkgname,
                         lib.loc=libname))
  OrganismDbi:::.loadOrganismDbiPkg(pkgname=pkgname,
                          graphData=graphData)
}


