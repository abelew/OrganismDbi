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
  checkTrue("GODb" == names(res)[res=='GOID'] )
  checkTrue("TranscriptDb" == names(res)[res=='TXID'] )
  checkTrue("OrgDb" == names(res)[res=='ENTREZID'] )
}

test_lookupDbsFromKeytype <- function(){
  res <- OrganismDbi:::.lookupDbsFromKeytype(x, "GOID")
  checkTrue(class(res)=="GODb")  
  res <- OrganismDbi:::.lookupDbsFromKeytype(x, "TXID")
  checkTrue(class(res)=="TranscriptDb")
  res <- OrganismDbi:::.lookupDbsFromKeytype(x, "ENTREZID")
  checkTrue(class(res)=="OrgDb")
}

test_lookupDbsFromKeytype <- function(){
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
  checkTrue("GODb" == names(res)[res=='TERM'] )
  checkTrue("TranscriptDb" == names(res)[res=='TXID'] )
  checkTrue("OrgDb" == names(res)[res=='ENTREZID'] )
}

test_lookupDbsFromCol <- function(){
  res <- OrganismDbi:::.lookupDbsFromCol(x,"TERM")
  checkTrue(class(res) == "GODb")
  res <- OrganismDbi:::.lookupDbsFromCol(x,"TXID")
  checkTrue(class(res) == "TranscriptDb")
  res <- OrganismDbi:::.lookupDbsFromCol(x,"ENTREZID")
  checkTrue(class(res) == "OrgDb")  
}


## x in this method is actually not an OrganismDb object, but a list of the
## AnnotationDb objects...
test_resortDbs <- function(){
  xlist <- list(x@GODb, x@TranscriptDb, x@OrgDb)
  res <- OrganismDbi:::.resortDbs(xlist)
  checkTrue(class(res[[2]]) == "OrgDb") 
}


test_lookupDbsFromCols <- function(){
  res <- OrganismDbi:::.lookupDbsFromCols(x,
                                          c("TXID","TERM","ENTREZID"),
                                          "GENEID")
  checkTrue(length(res)==3)
  checkTrue(class(res)=="list")
  checkTrue(class(res[[2]])=="OrgDb")
}


test_addAppropriateCols <- function(){
  cols <- c("TERM","ALIAS")
  db <- "OrgDb"
  dbs <-  OrganismDbi:::.lookupDbsFromCols(x,cols,"GENEID")
  res <- OrganismDbi:::.addAppropriateCols(db, dbs, cols)
  checkTrue("ENTREZID" %in% res)
  checkTrue("GO" %in% res)
  
  cols <- c("ACCNUM","ALIAS")
  dbs <-  OrganismDbi:::.lookupDbsFromCols(x,cols,"GENEID")
  res <- OrganismDbi:::.addAppropriateCols(db, dbs, cols)
  checkTrue("ENTREZID" %in% res)

  cols <- c("GO","ALIAS")
  dbs <-  OrganismDbi:::.lookupDbsFromCols(x,cols,"GENEID")
  res <- OrganismDbi:::.addAppropriateCols(db, dbs, cols)
  checkTrue("ENTREZID" %in% res)

  cols <- c("OBSOLETE")
  db <- "GODb"
  dbs <-  OrganismDbi:::.lookupDbsFromCols(x,cols,"GENEID")
  res <- OrganismDbi:::.addAppropriateCols(db, dbs, cols)
  checkTrue("TERM" %in% res)

  cols <- c("CDSID","CDSCHROM")
  db <- "TranscriptDb"
  dbs <-  OrganismDbi:::.lookupDbsFromCols(x,cols,"GENEID")
  res <- OrganismDbi:::.addAppropriateCols(db, dbs, cols)
  checkTrue("GENEID" %in% res)
}


test_getSelects <- function(){
  cols <- c("TERM","ALIAS")
  keytype <- "GENEID"
  dbs <-  OrganismDbi:::.lookupDbsFromCols(x,cols,keytype)
  keys <- head(keys(x, keytype),n=2)
  mkeys <- OrganismDbi:::.mkeys()
  res <- OrganismDbi:::.getSelects(dbs, keys, cols, keytype, mkeys)
  checkTrue(length(res)==3)
  checkTrue(class(res)=="list")
  checkTrue("GENEID" %in% colnames(res[[1]]))
  checkTrue("GO" %in% colnames(res[[2]]))
  checkTrue("Term" %in% colnames(res[[3]]))
  
  cols <- c("SYMBOL")
  keytype <- "OMIM"
  dbs <-  OrganismDbi:::.lookupDbsFromCols(x,cols,keytype)
  keys <- head(keys(x, keytype),n=2)
  res <- OrganismDbi:::.getSelects(dbs, keys, cols, keytype, mkeys)  
  checkTrue(length(res)==1)
  checkTrue(class(res)=="list")
  checkTrue("OMIM" %in% colnames(res[[1]]))
  checkTrue("SYMBOL" %in% colnames(res[[1]]))
}


test_mergeSelectResults <- function(){
  cols <- c("TERM","ALIAS")
  keytype <- "GENEID"
  dbs <-  OrganismDbi:::.lookupDbsFromCols(x,cols,keytype)
  keys <- head(keys(x, keytype),n=3)
  mkeys <- OrganismDbi:::.mkeys()
  sels <- OrganismDbi:::.getSelects(dbs, keys, cols, keytype, mkeys)
  
  res <- OrganismDbi:::.mergeSelectResults(sels, mkeys)
  checkTrue(dim(res)[2]==10)
  checkTrue(class(res)=="data.frame")
  checkTrue("GO" %in% colnames(res)) 
  checkTrue("GENEID" %in% colnames(res))
  checkTrue("ALIAS" %in% colnames(res))  
}


## MANY more tests
test_select <- function(){
  cols <- cols(x)[c(7,10,11,12)]
  keys <- head(keys(x, "ENTREZID"))
  keytype <- "ENTREZID"
  res <- OrganismDbi:::.select(x, keys, cols, keytype)
  checkTrue(dim(res)[2]==11)
  checkTrue(class(res)=="data.frame")
  checkTrue("GO" %in% colnames(res)) 
  checkTrue("ENTREZID" %in% colnames(res))
  checkTrue("ACCNUM" %in% colnames(res))  
  checkTrue("ALIAS" %in% colnames(res))

  cols <- cols(x)[c(7,10,11,37)]
  res <- OrganismDbi:::.select(x, keys, cols, keytype)
  checkTrue(dim(res)[2]==11)
  checkTrue("GO" %in% colnames(res)) 
  checkTrue("ENTREZID" %in% colnames(res))
  checkTrue("ACCNUM" %in% colnames(res))  
  checkTrue("CDSSTART" %in% colnames(res))
  
  cols <- cols(x)[c(7,8)]
  res <- OrganismDbi:::.select(x, keys, cols, keytype)
  checkTrue(dim(res)[2]==14)
  checkTrue("GO" %in% colnames(res)) 
  checkTrue("ENTREZID" %in% colnames(res))
  checkTrue("Term" %in% colnames(res))  
  
  cols <- cols(x)[c(10,11,37)]
  res <- OrganismDbi:::.select(x, keys, cols, keytype)
  checkTrue(dim(res)[2]==3)
  checkTrue("ENTREZID" %in% colnames(res))
  checkTrue("ACCNUM" %in% colnames(res))
  checkTrue("CDSSTART" %in% colnames(res))

  cols <- cols(x)[c(10,11,12)]
  res <- OrganismDbi:::.select(x, keys, cols, keytype)
  checkTrue(dim(res)[2]==3)
  checkTrue("ENTREZID" %in% colnames(res))
  checkTrue("ACCNUM" %in% colnames(res))
  checkTrue("ALIAS" %in% colnames(res))

  cols <- cols(x)[c(37,38)]
  res <- OrganismDbi:::.select(x, keys, cols, keytype)
  checkTrue(dim(res)[2]==3)
  checkTrue("ENTREZID" %in% colnames(res))
  checkTrue("CDSSTART" %in% colnames(res))
  checkTrue("CDSEND" %in% colnames(res))

  cols <- cols(x)[c(37)]
  res <- OrganismDbi:::.select(x, keys, cols, keytype)
  checkTrue(dim(res)[2]==2)
  checkTrue("ENTREZID" %in% colnames(res))
  checkTrue("CDSSTART" %in% colnames(res))

  cols <- cols(x)[c(10)]
  res <- OrganismDbi:::.select(x, keys, cols, keytype)
  checkTrue(dim(res)[2]==1)
  checkTrue("ENTREZID" %in% colnames(res))
  
}


## TODO: I may want to change what I am passing to reqCols when calling
## AnnotationDbi:::.resort().  This may be a better place than elsewhere to
## clean up the columns.  But it may not be a good way to get rid of
## duplicated columns...

## TODO: add more tests from other kinds of IDs
## Right now there seem to be some performance issues for these.
## test_select__otherIDTypes <- function(){

## cols <- cols(x)[c(7,12)]
## keys <- head(keys(x, "ALIAS"))
## keytype <- "ALIAS"
## res <- OrganismDbi:::.select(x, keys, cols, keytype)

## debug(OrganismDbi:::.getSelects)
## debug(OrganismDbi:::.mergeSelectResults)
## debug(AnnotationDbi:::.resort)


##   cols <- cols(x)[c(7,10,11,12)]
##   keys <- head(keys(x, "ALIAS"))
##   keytype <- "ALIAS"
##   res <- OrganismDbi:::.select(x, keys, cols, keytype)

##   checkTrue(dim(res)[2]==11)
##   checkTrue(class(res)=="data.frame")
##   checkTrue("GO" %in% colnames(res)) 
##   checkTrue("ENTREZID" %in% colnames(res))
##   checkTrue("ACCNUM" %in% colnames(res))  
##   checkTrue("ALIAS" %in% colnames(res))
 
## }









## library(OrganismDbi); OrganismDbi:::.test()
