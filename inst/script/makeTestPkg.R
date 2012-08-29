## script to test my package code generator
require(OrganismDbi)
version = "1.0.0"

## for human
gd <- data.frame(xDbs=c("GO.db","org.Hs.eg.db"),
                 yDbs=c("org.Hs.eg.db","TxDb.Hsapiens.UCSC.hg19.knownGene"),
                 xKeys=c("GOID","ENTREZID"),
                 yKeys= c("GO","GENEID"))

makeOrganismPackage(pkgname = "Homo.sapiens",
                    graphData = gd,
                    organism = "Homo sapiens",
                    version = version,
                    maintainer =
              "Bioconductor Package Maintainer <maintainer@bioconductor.org>",
                    author = "Bioconductor Core Team",
                    destDir = ".",
                    license = "Artistic-2.0")



## for mouse
gd <- data.frame(xDbs = c("GO.db","org.Mm.eg.db"),
                 yDbs = c("org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.ensGene"),
                 xKeys = c("GOID","ENSEMBL"),
                 yKeys = c("GO","GENEID"))

makeOrganismPackage(pkgname = "Mus.musculus",
                    graphData = gd,
                    organism = "Mus musculus",
                    version = version,
                    maintainer =
              "Bioconductor Package Maintainer <maintainer@bioconductor.org>",
                    author = "Bioconductor Core Team",
                    destDir = ".",
                    license = "Artistic-2.0")


## for rat
gd <- data.frame(xDbs = c("GO.db","org.Rn.eg.db"),
                 yDbs = c("org.Rn.eg.db","TxDb.Rnorvegicus.UCSC.rn4.ensGene"),
                 xKeys = c("GOID","ENSEMBL"),
                 yKeys = c("GO","GENEID"))

makeOrganismPackage(pkgname = "Rattus.norvegicus",
                    graphData = gd,
                    organism = "Rattus norvegicus",
                    version = "1.0.0",
                    maintainer =
              "Bioconductor Package Maintainer <maintainer@bioconductor.org>",
                    author = "Bioconductor Core Team",
                    destDir = ".",
                    license = "Artistic-2.0")
