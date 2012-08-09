xDbs <- c("GO.db","org.Hs.eg.db")
yDbs <- c("org.Hs.eg.db","TxDb.Hsapiens.UCSC.hg19.knownGene")
xKeys <- c("GOID","ENTREZID")
yKeys <- c("GO","GENEID")
gd <- data.frame(cbind(xDbs, yDbs, xKeys, yKeys))
require("RUnit")
require("Homo.sapiens")


test_extractPkgsAndCols <- function(){
  res <- OrganismDbi:::.extractPkgsAndCols(gd)
  checkEquals(names(res), c("GO.db","org.Hs.eg.db","org.Hs.eg.db",
                          "TxDb.Hsapiens.UCSC.hg19.knownGene"))
  checkTrue(all(res %in%  c("GOID","ENTREZID","GO","GENEID")))
}

test_testKeys <- function(){
  res <- OrganismDbi:::.extractPkgsAndCols(gd)
  ## the following should not blow up
  OrganismDbi:::.testKeys(res)

  ## the following should blow up
  res2 <- c("GO.db","org.Hs.eg.db","org.Hs.eg.db",
            "TxDb.Hsapiens.UCSC.hg19.knownGene")
  names(res2) <- c("WRONGVALUE","ENTREZID","GO","GENEID")
  checkException(OrganismDbi:::.testKeys(res2))

  ## note that the case of a bad DB would have been caught elsewhere
  ## (and earlier)
}

