## .onLoad gets the data.frame from the /data directory
.onLoad <- function(libname, pkgname) {
  load(system.file("data","graphData.Rda",package=pkgname,
                         lib.loc=libname))
  OrganismDbi:::.loadOrganismDbiPkg(pkgname=pkgname,
                          graphData=graphData)
}


