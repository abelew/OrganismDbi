## script to test my package code generator
require(OrganismDbi)
## makeOrganismPackage(pkgname = "Mus.musculus",
##                                 OrgPkg = "org.Mm.eg.db",
##                                 TxDbPkg = "TxDb.Mmusculus.UCSC.mm9.knownGene",
##                                 version = "1.0.0",
##                                 maintainer = "me,<me@someplace.org>",
##                                 author = "me",
##                                 destDir = ".",
##                                 license = "Artistic-2.0")

makeOrganismPackage(pkgname = "Homo.sapiens",
                                OrgPkg = "org.Hs.eg.db",
                                TxDbPkg = "TxDb.Hsapiens.UCSC.hg19.knownGene",
                                version = "1.0.0",
                                maintainer = "Bioconductor Package Maintainer <maintainer@bioconductor.org>",
                                author = "Bioconductor Core Team",
                                destDir = ".",
                                license = "Artistic-2.0")
