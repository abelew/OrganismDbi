## unit tests for the "meta-select"
## These tests are testing the software and not the indiv. packages
## Will base testing on humans for now
require("Homo.sapiens")
x <- Homo.sapiens
require("RUnit")


test_keytypes <- function(){
  res <- OrganismDbi:::.keytypes(x)
  checkTrue("GOID" %in% res)
  checkTrue("TXID" %in% res)
  checkTrue("ENTREZID" %in% res)
}

test_makekeytypeMapping <- function(){
  res <- OrganismDbi:::.makekeytypeMapping(x)
  checkTrue("GO.db" == names(res)[res=='GOID'] )
  checkTrue("TxDb.Hsapiens.UCSC.hg19.knownGene" == names(res)[res=='TXID'] )
  checkTrue("org.Hs.eg.db" == names(res)[res=='ENTREZID'] )
}

test_lookupDbFromKeytype <- function(){
  res <- OrganismDbi:::.lookupDbFromKeytype(x, "GOID")
  checkTrue(class(res)=="GODb")  
  res <- OrganismDbi:::.lookupDbFromKeytype(x, "TXID")
  checkTrue(class(res)=="TranscriptDb")
  res <- OrganismDbi:::.lookupDbFromKeytype(x, "ENTREZID")
  checkTrue(class(res)=="OrgDb")
}

test_lookupDbFromKeytype <- function(){
  res <- OrganismDbi:::.keys(x, "GOID")
  checkTrue(is.character(head(res)))
  checkTrue(length(res) > 33000)

  res <- OrganismDbi:::.keys(x, "TXID")
  checkTrue(is.integer(head(res)))
  checkTrue(length(res) > 70000)

  res <- OrganismDbi:::.keys(x, "ENTREZID")
  checkTrue(is.character(head(res)))
  checkTrue(length(res) > 42000)
}


test_makecolMapping <- function(){
  res <- OrganismDbi:::.makecolMapping(x)
  checkTrue("GO.db" == names(res)[res=='TERM'] )
  checkTrue("TxDb.Hsapiens.UCSC.hg19.knownGene" == names(res)[res=='TXID'] )
  checkTrue("org.Hs.eg.db" == names(res)[res=='ENTREZID'] )
}


test_resortDbs <- function(){
  cls <- c("GOID","PATH")
  pkgs <- unique(OrganismDbi:::.lookupDbNamesFromCols(x, cls))
  keytype <- "ENTREZID"
  res <- OrganismDbi:::.resortDbs(x, pkgs, keytype)
  checkTrue(length(res) == 2)
  checkTrue(names(res)[1] == "org.Hs.eg.db") ## shuold always start here
  checkTrue(names(res)[2] == "GO.db") ## should always go here

  cls <- c("GOID","PATH", "TXNAME")
  pkgs <- unique(OrganismDbi:::.lookupDbNamesFromCols(x, cls))
  keytype <- "TXNAME"
  res <- OrganismDbi:::.resortDbs(x, pkgs, keytype)
  checkTrue(length(res) == 3)
  checkTrue(names(res)[1] == "TxDb.Hsapiens.UCSC.hg19.knownGene")
  checkTrue(names(res)[2] == "org.Hs.eg.db") 
  checkTrue(names(res)[3] == "GO.db")

}



test_getForeignEdgeKeys <- function(){
  keytype <- c("ENTREZID")
  cls <- c("OBSOLETE") ## this is a defunct ID
  ## And this should therefore throw an error
  checkException(OrganismDbi:::.getForeignEdgeKeys(x,cls,"GENEID"))

  ## This method needs to be able to interpolate between nodes.
  ## So this should work out OK:
  keytype <- c("TXNAME")
  cls <- c("GOID")
  res <- OrganismDbi:::.getForeignEdgeKeys(x, cls, keytype)
  checkTrue("GENEID" == unlist(res)[1]) 
  checkTrue("ENTREZID" == unlist(res)[2])
  checkTrue("GO" == unlist(res)[3]) 
  checkTrue("GOID" == unlist(res)[4])

  keytype <- c("ENTREZID")
  cls <- c("GOID")
  res <- OrganismDbi:::.getForeignEdgeKeys(x, cls, keytype)
  checkTrue("GO" == unlist(res)[1])
  checkTrue("GOID" == unlist(res)[2])
  
  keytype <- c("TXNAME")
  cls <- c("GO")
  res <- OrganismDbi:::.getForeignEdgeKeys(x, cls, keytype)
  checkTrue("GENEID" == unlist(res)[1])
  checkTrue("ENTREZID" == unlist(res)[2])

  ## Here is a case for when the start node is in the middle of the graph 
  keytype <- c("ENTREZID") 
  cls = c("GOID","SYMBOL","TXNAME") 
  res <- OrganismDbi:::.getForeignEdgeKeys(x, cls, keytype) 
  checkTrue("GO" == unlist(res)[1]) 
  checkTrue("GOID" == unlist(res)[2]) 
  checkTrue("ENTREZID" == unlist(res)[3]) 
  checkTrue("GENEID" == unlist(res)[4]) 
} 

test_mkeys <- function(){
  tbl1 = "TxDb.Hsapiens.UCSC.hg19.knownGene"
  tbl2 = "org.Hs.eg.db"
  res <- OrganismDbi:::.mkeys(x, tbl1, tbl2, key="tbl1")
  checkTrue("GENEID"==res)
  res <- OrganismDbi:::.mkeys(x, tbl1, tbl2, key="tbl2")
  checkTrue("ENTREZID"==res)

  tbl1 = "GO.db"
  tbl2 = "org.Hs.eg.db"
  res <- OrganismDbi:::.mkeys(x, tbl1, tbl2, key="tbl1")
  checkTrue("GOID"==res)
  res <- OrganismDbi:::.mkeys(x, tbl1, tbl2, key="tbl2")
  checkTrue("GO"==res)

  tbl1 = "GO.db"
  tbl2 = "org.Hs.eg.db"
  res <- OrganismDbi:::.mkeys(x, tbl1, tbl2, key="both")
  checkEquals(res, c("GOID","GO"))
}


test_getSelects <- function(){
  cls <- c("TERM", "ALIAS")
  kt <- "GENEID"
  cls <- OrganismDbi:::.getForeignEdgeKeys(x, cls, kt)
  dbs <-  OrganismDbi:::.lookupDbsFromCols(x, cls, kt)  
  keys <- head(keys(x, kt), n=2)
  res <- OrganismDbi:::.getSelects(x, dbs, keys, cls, kt)
  checkTrue(length(res)==3)
  checkTrue(class(res)=="list")
  checkTrue("GENEID" %in% colnames(res[[1]]))
  checkTrue("GO" %in% colnames(res[[2]]))
  checkTrue("TERM" %in% colnames(res[[3]]))
  
  cls <- c("SYMBOL")
  kt <- "OMIM"
  cls <- OrganismDbi:::.getForeignEdgeKeys(x, cls, kt)
  dbs <-  OrganismDbi:::.lookupDbsFromCols(x, cls, kt)   
  keys <- head(keys(x, kt), n=2)
  res <- OrganismDbi:::.getSelects(x, dbs, keys, cls, kt)  
  checkTrue(length(res)==1)
  checkTrue(class(res)=="list")
  checkTrue("OMIM" %in% colnames(res[[1]]))
  checkTrue("SYMBOL" %in% colnames(res[[1]]))
}

test_mergeSelectResults <- function(){
  cls <- c("TERM", "ALIAS")
  kt <- "GENEID"
  cls <- OrganismDbi:::.getForeignEdgeKeys(x, cls, kt)
  dbs <-  OrganismDbi:::.lookupDbsFromCols(x, cls, kt)  
  keys <- head(keys(x, kt), n=3)
  sels <- OrganismDbi:::.getSelects(x, dbs, keys, cls, kt)

  res <- OrganismDbi:::.mergeSelectResults(x, sels)
  checkTrue(dim(res)[2]==6)
  checkTrue(class(res)=="data.frame")
  checkTrue("GO" %in% colnames(res)) 
  checkTrue("GENEID" %in% colnames(res))
  checkTrue("ALIAS" %in% colnames(res))  
}
## TODO: Investigate why the above query has so many NA rows for the gene related data???


## MANY more tests
test_select <- function(){
  cls <- c("GO","ALIAS","CHR","CHRLOC")
  keys <- head(keys(x, "ENTREZID"))
  keytype <- "ENTREZID"
  res <- OrganismDbi:::.select(x, keys, cls, keytype)
  checkTrue(dim(res)[2]==8)
  checkTrue(class(res)=="data.frame")
  checkTrue("GO" %in% colnames(res)) 
  checkTrue("EVIDENCE" %in% colnames(res)) 
  checkTrue("ENTREZID" %in% colnames(res))
  checkTrue("CHRLOCCHR" %in% colnames(res))  
  checkTrue("ALIAS" %in% colnames(res))

  cls <- c("IPI", "ALIAS", "CDSSTART") 
  res <- OrganismDbi:::.select(x, keys, cls, keytype) 
  checkTrue(dim(res)[2]==4) 
  checkTrue("IPI" %in% colnames(res)) 
  checkTrue("ENTREZID" %in% colnames(res)) 
  checkTrue("ALIAS" %in% colnames(res)) 
  checkTrue("CDSSTART" %in% colnames(res)) 

  ## bugged because I have to change the way I sort these things so that they always start with the keytype dbs...
  cls <- c("GOID","ENTREZID")
  res <- OrganismDbi:::.select(x, keys, cls, keytype)
  checkTrue(dim(res)[2]==4)
  checkTrue("GO" %in% colnames(res)) 
  checkTrue("ENTREZID" %in% colnames(res))
 
  cls <- c("ALIAS","CHR","EXONNAME")
  res <- OrganismDbi:::.select(x, keys, cls, keytype) 
  checkTrue(dim(res)[2]==4) 
  checkTrue("ALIAS" %in% colnames(res)) 
  checkTrue("ENTREZID" %in% colnames(res)) 
  checkTrue("EXONNAME" %in% colnames(res)) 
  
  cls <- c("ACCNUM","CDSSTART") 
  res <- OrganismDbi:::.select(x, keys, cls, keytype)
  checkTrue(dim(res)[2]==3)
  checkTrue("ENTREZID" %in% colnames(res))
  checkTrue("ACCNUM" %in% colnames(res))
  checkTrue("CDSSTART" %in% colnames(res))

  cls <- c("ACCNUM", "ALIAS")
  res <- OrganismDbi:::.select(x, keys, cls, keytype)
  checkTrue(dim(res)[2]==3)
  checkTrue("ENTREZID" %in% colnames(res))
  checkTrue("ACCNUM" %in% colnames(res))
  checkTrue("ALIAS" %in% colnames(res))

  cls <- c("CDSSTART","CDSEND")
  res <- OrganismDbi:::.select(x, keys, cls, keytype)
  checkTrue(dim(res)[2]==3)
  checkTrue("ENTREZID" %in% colnames(res))
  checkTrue("CDSSTART" %in% colnames(res))
  checkTrue("CDSEND" %in% colnames(res))

  cls <- c("CDSSTART")
  res <- OrganismDbi:::.select(x, keys, cls, keytype)
  checkTrue(dim(res)[2]==2)
  checkTrue("ENTREZID" %in% colnames(res))
  checkTrue("CDSSTART" %in% colnames(res))

  cls <- c("ENTREZID")
  res <- OrganismDbi:::.select(x, keys, cls, keytype)
  checkTrue(dim(res)[2]==1)
  checkTrue("ENTREZID" %in% colnames(res))


  
  ## Why do the following cases not work?
  ## Because .getSelects also (and maybe even mergeSelectResults) need to walk
  ## along the graph intelligently (and not just based on distances).
  ## In this case, two nodes were equal distance from the central node, so
  ## just sorting them based on distance will not put them into the correct
  ## order for processing.  I need to preserve the order from my initial walk
  ## (which happened when I added appropriate cols to gather the foreign
  ## keys).  And I need to follow that ordering in these other methods.
  ## OR I could generalize the.nodeWalker code even more to take 3 walks
  ## OR maybe I still don't need to be calling .mkeys like this?
  ## OR maybe I just need to be a tiny bit smarter about what I pass into
  ## these other two methods (so that the keys are grouped into edges)
  
  keys <- head(keys(x, "ENTREZID"))
  keytype <- "ENTREZID"
  cls = c("GOID" ,  "SYMBOL", "TXNAME")
  res <- select(Homo.sapiens, keys, cls, keytype)

##   ## same exact bug comes up when I try to use a homology package
##   cls = c("ALIAS", "ORYZA_SATIVA")
##   res <- select(Homo.sapiens, keys, cls, keytype)
  
}


## Also need to test a species with other keys to join DBs
require("Rattus.norvegicus") 
r <- Rattus.norvegicus 


test_rattus <- function(){ 
  cls <- c("GO","ALIAS","CHR") 
  k <- head(keys(r, "ENTREZID")) 
  keytype <- "ENTREZID" 
  res <- OrganismDbi:::.select(r, k, cls, keytype) 
  checkTrue(dim(res)[2]==6) 
  checkTrue("ENTREZID" %in% colnames(res)) 
  checkTrue("GO" %in% colnames(res)) 
  checkTrue("ONTOLOGY" %in% colnames(res)) 
  checkTrue("CHR" %in% colnames(res)) 
  checkTrue("ALIAS" %in% colnames(res)) 

  cls <- c("GO","ALIAS","CHR","TXNAME") 
  res <- OrganismDbi:::.select(r, k, cls, keytype) 
  checkTrue(dim(res)[2]==7) 
  checkTrue("ENTREZID" %in% colnames(res)) 
  checkTrue("GO" %in% colnames(res)) 
  checkTrue("ALIAS" %in% colnames(res)) 
  checkTrue("CHR" %in% colnames(res)) 
  checkTrue("TXNAME" %in% colnames(res)) 

  cls <- c("CHR","TXNAME") 
  res <- OrganismDbi:::.select(r, k, cls, keytype) 
  checkTrue(dim(res)[2]==3) 
  checkTrue("ENTREZID" %in% colnames(res)) 
  checkTrue("CHR" %in% colnames(res)) 
  checkTrue("TXNAME" %in% colnames(res)) 

  ## now test different keytype
  k = head(keys(r, keytype="ENSEMBL"))
  keytype="ENSEMBL"
  cls <- c("GO","ALIAS","CHR","TXNAME") 
  res <- OrganismDbi:::.select(r, k, cls, keytype) 
  checkTrue(dim(res)[2]==7) 
  checkTrue("ENSEMBL" %in% colnames(res)) 
  checkTrue("GO" %in% colnames(res)) 
  checkTrue("ONTOLOGY" %in% colnames(res)) 
  checkTrue("CHR" %in% colnames(res)) 
  checkTrue("ALIAS" %in% colnames(res)) 
  checkTrue("TXNAME" %in% colnames(res)) 

  ## now test key that starts us from TxDb
  k = head(keys(r, keytype="TXNAME"))
  keytype="TXNAME"
  cls <- c("GO","ALIAS","CHR") 
  res <- OrganismDbi:::.select(r, k, cls, keytype) 
  checkTrue(dim(res)[2]==6) 
  checkTrue("TXNAME" %in% colnames(res)) 
  checkTrue("GO" %in% colnames(res)) 
  checkTrue("EVIDENCE" %in% colnames(res)) 
  checkTrue("CHR" %in% colnames(res)) 
  checkTrue("ALIAS" %in% colnames(res)) 

  ## now test key that starts us from Go
  ## TODO: bug??? row of NAs for entire 1st line...
  k = head(keys(r, keytype="GOID"))
  keytype="GOID"
  cls <- c("GOID","ALIAS","CHR") 
  res <- OrganismDbi:::.select(r, k, cls, keytype) 
  checkTrue(dim(res)[2]==5) 
  checkTrue("GOID" %in% colnames(res)) 
  checkTrue("ONTOLOGY" %in% colnames(res)) 
  checkTrue("CHR" %in% colnames(res)) 
  checkTrue("ALIAS" %in% colnames(res))


##   ## Also I see the bug here:
##   ## what happens when we use a key from the middle?
##   ## Our algorithm will walk 
##   k = keys=head(keys(r,keytype="ENTREZID"))
##   cls = c("GOID","SYMBOL","TXNAME")
##   keytype = "ENTREZID"
##   res <- OrganismDbi:::.select(r, k, cls, keytype)
  
  
} 











