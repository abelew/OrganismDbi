gd <- list(join1 = c(GO.db="GOID", org.Hs.eg.db="GO"),
           join2 = c(org.Hs.eg.db="ENTREZID",
                     TxDb.Hsapiens.UCSC.hg19.knownGene="GENEID"))

require("Homo.sapiens")

test_extractPkgsAndCols <- function(){
  gdm <- OrganismDbi:::.mungeGraphData(gd)
  res <- OrganismDbi:::.extractPkgsAndCols(gdm)
  checkEquals(names(res), c("GO.db","org.Hs.eg.db","org.Hs.eg.db",
                          "TxDb.Hsapiens.UCSC.hg19.knownGene"))
  checkTrue(all(res %in%  c("GOID","ENTREZID","GO","GENEID")))
}

test_testKeys <- function(){
  gdm <- OrganismDbi:::.mungeGraphData(gd)
  res <- OrganismDbi:::.extractPkgsAndCols(gdm)
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

test_TxDb <- function(){
    txdbMouse <- loadDb(system.file('extdata','myTxDb.sqlite',
                                    package='OrganismDbi'))
    library(TxDb.Mmusculus.UCSC.mm9.knownGene)
    txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
    odb <- makeOrganismDbFromTxDb(txdb=txdbMouse)
    res <- resources(odb)
    checkTrue(res[['TxDb.Mmusculus.UCSC.mm9.knownGene']] == dbfile(txdbMouse))
    ## NOW change it:
    TxDb(odb) <- txdb
    res <- resources(odb)
    checkTrue(res[['TxDb.Mmusculus.UCSC.mm9.knownGene']] == dbfile(txdb))
}


## Above stuff works BUT: this does NOT work:
##  transcript_ids <- c(
##      "uc009uzf.1",
##      "uc009uzg.1",
##      "uc009uzh.1",
##      "uc009uzi.1",
##      "uc009uzj.1"
##  )
##  txdbMouse <- makeTxDbFromUCSC(genome="mm9", tablename="knownGene",
##                            transcript_ids=transcript_ids)
## odb <- makeOrganismDbFromTxDb(txdb=txdbMouse)
## TxDb(odb)
