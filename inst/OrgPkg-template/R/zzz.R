## I want this to be "minimalist"
assign("@PKGNAME@", setClass("@PKGNAME@", contains="AnnotationOrganismDb"))

.onLoad <- function(libname, pkgname) {
  OrganismDbi:::.loadOrganismDbiPkg(pkgname=pkgname,
                          GOPkg="GO.db",
                          OrgPkg="@ORGPKG@",
                          TranscriptPkg="@TXDBPKG@")
}
