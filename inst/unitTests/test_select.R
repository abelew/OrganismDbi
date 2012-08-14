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


## test_resortDbs <- function(){
##   cls <- c("GOID","PATH")
##   pkgs <- unique(OrganismDbi:::.lookupDbNamesFromCols(x, cls))
##   keytype <- "ENTREZID"
##   res <- OrganismDbi:::.resortDbs(x, pkgs, keytype)
##   checkTrue(length(res) == 2)
##   checkTrue(names(res)[1] == "org.Hs.eg.db") ## shuold always start here
##   checkTrue(names(res)[2] == "GO.db") ## should always go here

##   cls <- c("GOID","PATH", "TXNAME")
##   pkgs <- unique(OrganismDbi:::.lookupDbNamesFromCols(x, cls))
##   keytype <- "TXNAME"
##   res <- OrganismDbi:::.resortDbs(x, pkgs, keytype)
##   checkTrue(length(res) == 3)
##   checkTrue(names(res)[1] == "TxDb.Hsapiens.UCSC.hg19.knownGene")
##   checkTrue(names(res)[2] == "org.Hs.eg.db") 
##   checkTrue(names(res)[3] == "GO.db")

## } 



test_getForeignEdgeKeys <- function(){
  keytype <- c("ENTREZID")
  cls <- c("OBSOLETE") ## this is a defunct ID
  ## And this should therefore throw an error
  checkException(OrganismDbi:::.getForeignEdgeKeys(x,cls,"GENEID"))

  ## What if there is only DBS required?
  keytype <- c("TXNAME")
  cls <- c("TXID")
  res <- OrganismDbi:::.getForeignEdgeKeys(x, cls, keytype)
  checkTrue(length(res) == 0) ## Should give me an empty list (IOW no edges)
  
  ## This method needs to be able to interpolate between nodes.
  ## So this should work out OK:
  keytype <- c("TXNAME")
  cls <- c("GOID")
  res <- OrganismDbi:::.getForeignEdgeKeys(x, cls, keytype)
  checkTrue("GENEID" == res[1]) 
  checkTrue("ENTREZID" == res[2])
  checkTrue("GO" == res[3]) 
  checkTrue("GOID" == res[4])

  keytype <- c("ENTREZID")
  cls <- c("GOID")
  res <- OrganismDbi:::.getForeignEdgeKeys(x, cls, keytype)
  checkTrue("GO" == res[1])
  checkTrue("GOID" == res[2])
  
  keytype <- c("TXNAME")
  cls <- c("GO")
  res <- OrganismDbi:::.getForeignEdgeKeys(x, cls, keytype)
  checkTrue("GENEID" == res[1])
  checkTrue("ENTREZID" == res[2])

  ## Here is a case for when the start node is in the middle of the graph 
  keytype <- c("ENTREZID") 
  cls = c("GOID","SYMBOL","TXNAME") 
  res <- OrganismDbi:::.getForeignEdgeKeys(x, cls, keytype) 
  checkTrue("GO" == res[1]) 
  checkTrue("GOID" == res[2]) 
  checkTrue("ENTREZID" == res[3]) 
  checkTrue("GENEID" == res[4]) 
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
  res2 <- c("GOID","GO")
  names(res2) <- c("GO.db","org.Hs.eg.db")
  checkEquals(res, res2)
} 


test_getSelects <- function(){
  cls <- c("TERM", "ALIAS")
  kt <- "GENEID"
  fkys <- OrganismDbi:::.getForeignEdgeKeys(x, cls, kt)
  dbs <-  OrganismDbi:::.lookupDbsFromFkeys(x, fkys, kt)  ## reverse order???
  keys <- head(keys(x, kt), n=2)
  cols <- unique(c(kt, cls, fkys))
  res <- OrganismDbi:::.getSelects(x, dbs, keys, cols, kt)
  checkTrue(length(res)==3)
  checkTrue(class(res)=="list")
  checkTrue("GENEID" %in% colnames(res[[1]]))
  checkTrue("GO" %in% colnames(res[[2]]))
  checkTrue("TERM" %in% colnames(res[[3]]))
  
  cls <- c("SYMBOL")
  kt <- "OMIM"
  fkys <- OrganismDbi:::.getForeignEdgeKeys(x, cls, kt)
  dbs <-  OrganismDbi:::.lookupDbsFromFkeys(x, fkys, kt)   
  keys <- head(keys(x, kt), n=2)
  cols <- unique(c(kt, cls, fkys))
  res <- OrganismDbi:::.getSelects(x, dbs, keys, cols, kt)  
  checkTrue(length(res)==1)
  checkTrue(class(res)=="list")
  checkTrue("OMIM" %in% colnames(res[[1]]))
  checkTrue("SYMBOL" %in% colnames(res[[1]]))

  ## Then there is this case:
  cls = c("GOID" ,  "SYMBOL", "TXNAME")
  kt <- "ENTREZID"
  fkys <- OrganismDbi:::.getForeignEdgeKeys(x, cls, kt)

  ## This shows that in fact .lookupDbNamesFromCols is not fully respecting
  ## the order of fkys...
  pkgs <- OrganismDbi:::.lookupDbNamesFromCols(x, fkys)
  ## And then .lookupDbsFromFkeys is calling unique, which is something we
  ## want it to do, but ONLY when the same thing is repeated twice IN A ROW.
  ## So I will need to call a custom method from .lookupDbsFromFkeys()
  
  
  dbs <-  OrganismDbi:::.lookupDbsFromFkeys(x, fkys, kt)   
  keys <- head(keys(x, "ENTREZID"))
  cols <- unique(c(kt, cls, fkys))
  res <- OrganismDbi:::.getSelects(x, dbs, keys, cols, kt)  
  ## fails not because of path issue (resolved) but because of problem with a
  ## need to reset whenever we start gathering from the root again...
  
}

test_mergeSelectResults <- function(){
  cls <- c("TERM", "ALIAS")
  kt <- "GENEID"
  fkys <- OrganismDbi:::.getForeignEdgeKeys(x, cls, kt)
  dbs <-  OrganismDbi:::.lookupDbsFromFkeys(x, fkys, kt)  
  keys <- head(keys(x, kt), n=3)
  cols <- unique(c(kt, cls, fkys))
  sels <- OrganismDbi:::.getSelects(x, dbs, keys, cols, kt)

  res <- OrganismDbi:::.mergeSelectResults(x, sels)
  checkTrue(dim(res)[2]==6)
  checkTrue(class(res)=="data.frame")
  checkTrue("GO" %in% colnames(res)) 
  checkTrue("GENEID" %in% colnames(res))
  checkTrue("ALIAS" %in% colnames(res))  
}


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

  ## So the path problem sorting problem above has been dealt with.  But now
  ## what is remaining is just that sometimes when we are getting data via
  ## getSelects, we are not getting the same thing in a row.  IOW, sometimes
  ## we are starting a new walk (to leaves and from root node).  When this
  ## happens, .getSelects() needs to somehow know about it so that it can
  ## behave differently.



  ## one idea from here is:  
  ## So actually .lookupDbsFromFkeys() needs to NOT call unique before
  ## getSelects gets to it, since it needs to repeat the root node (for
  ## example) whenever it sees that node again so that later on I can use that
  ## node to know I have started another path...  The names on this modified
  ## vector will then become the set of pkgs to get data from.

  ## For .getSelects, I will just need to track each node as they are added, so
  ## that I can be sure not to add the same node twice.  This is OK, because
  ## my cols statement will be general enough to ensure that I get all
  ## relevant data (and keys) for each node the 1st time I hit it (because the
  ## cols have been pre-computed during the walk)

  ## And then for .mergeSelects, I will have one last problem where I can have
  ## missing gaps (on account of having dropped the nodes in .getSelects
  ## above).
  
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











