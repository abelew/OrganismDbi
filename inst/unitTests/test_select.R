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
  checkTrue("GODb" == names(res)[res=='TERM'] )
  checkTrue("TranscriptDb" == names(res)[res=='TXID'] )
  checkTrue("OrgDb" == names(res)[res=='ENTREZID'] )
}



## x in this method is actually not an OrganismDb object, but a list of the
## AnnotationDb objects...
test_resortDbs <- function(){
  xlist <- list(x@GODb, x@TranscriptDb, x@OrgDb)
  res <- OrganismDbi:::.resortDbs(xlist)
  checkTrue(class(res[[2]]) == "OrgDb") 
}



test_addAppropriateCols <- function(){
##   cls <- c("TERM","ALIAS")
##   db <- "OrgDb"
##   dbs <-  OrganismDbi:::.lookupDbsFromCols(x,cls,"GENEID")
##   res <- OrganismDbi:::.addAppropriateCols(db, dbs, cls)
##   checkTrue("ENTREZID" %in% res)
##   checkTrue("GO" %in% res)
  
##   cls <- c("ACCNUM","ALIAS")
##   dbs <-  OrganismDbi:::.lookupDbsFromCols(x,cls,"GENEID")
##   res <- OrganismDbi:::.addAppropriateCols(db, dbs, cls)
##   checkTrue("ENTREZID" %in% res)

##   cls <- c("GO","ALIAS")
##   dbs <-  OrganismDbi:::.lookupDbsFromCols(x,cls,"GENEID")
##   res <- OrganismDbi:::.addAppropriateCols(db, dbs, cls)
##   checkTrue("ENTREZID" %in% res)

##   cls <- c("OBSOLETE") ## this is a defunct ID
##   db <- "GODb"
##   ## And this should therefore throw an error
##   checkException(OrganismDbi:::.lookupDbsFromCols(x,cls,"GENEID"))

##   cls <- c("CDSID","CDSCHROM")
##   db <- "TranscriptDb"
##   dbs <-  OrganismDbi:::.lookupDbsFromCols(x,cls,"GENEID")
##   res <- OrganismDbi:::.addAppropriateCols(db, dbs, cls)
##   checkTrue("GENEID" %in% res)

  keytype <- c("TXNAME")
  cls <- c("GOID")
  res <- OrganismDbi:::.addAppropriateCols(x, cls, keytype)
  checkTrue("GENEID" %in% res)
}


test_getSelects <- function(){
  cls <- c("TERM","ALIAS")
  keytype <- "GENEID"
  ## dbs now has names to indicate which package for each
  dbs <-  OrganismDbi:::.lookupDbsFromCols(x,cls,keytype)
  keys <- head(keys(x, keytype),n=2)
#  mkeys <- OrganismDbi:::.mkeys()
  res <- OrganismDbi:::.getSelects(dbnames, keys, cls, keytype)
  checkTrue(length(res)==3)
  checkTrue(class(res)=="list")
  checkTrue("GENEID" %in% colnames(res[[1]]))
  checkTrue("GO" %in% colnames(res[[2]]))
  checkTrue("TERM" %in% colnames(res[[3]]))
  
  cls <- c("SYMBOL")
  keytype <- "OMIM"
  dbs <-  OrganismDbi:::.lookupDbsFromCols(x,cls,keytype)
  keys <- head(keys(x, keytype),n=2)
  res <- OrganismDbi:::.getSelects(dbs, keys, cls, keytype, mkeys)  
  checkTrue(length(res)==1)
  checkTrue(class(res)=="list")
  checkTrue("OMIM" %in% colnames(res[[1]]))
  checkTrue("SYMBOL" %in% colnames(res[[1]]))
}


test_mergeSelectResults <- function(){
  cls <- c("TERM","ALIAS")
  keytype <- "GENEID"
  dbs <-  OrganismDbi:::.lookupDbsFromCols(x,cls,keytype)
  keys <- head(keys(x, keytype),n=3)
  mkeys <- OrganismDbi:::.mkeys()
  sels <- OrganismDbi:::.getSelects(dbs, keys, cls, keytype, mkeys)
  
  res <- OrganismDbi:::.mergeSelectResults(sels, mkeys)
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

  cls <- c("GOID","ENTREZID")
  res <- OrganismDbi:::.select(x, keys, cls, keytype)
  checkTrue(dim(res)[2]==5)
  checkTrue("GO" %in% colnames(res)) 
  checkTrue("ENTREZID" %in% colnames(res))
  checkTrue("TERM" %in% colnames(res))  
 
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
  
}


## TODO: I may want to change what I am passing to reqCols when calling
## AnnotationDbi:::.resort().  This may be a better place than elsewhere to
## clean up the columns.  But it may not be a good way to get rid of
## duplicated columns...

## TODO: add more tests from other kinds of IDs
## Right now there seem to be some performance issues for these.
## test_select__otherIDTypes <- function(){

## cls <- cols(x)[c(7,12)]
## keys <- head(keys(x, "ALIAS"))
## keytype <- "ALIAS"
## res <- OrganismDbi:::.select(x, keys, cls, keytype)

## debug(OrganismDbi:::.getSelects)
## debug(OrganismDbi:::.mergeSelectResults)
## debug(AnnotationDbi:::.resort)


##   cls <- cols(x)[c(7,10,11,12)]
##   keys <- head(keys(x, "ALIAS"))
##   keytype <- "ALIAS"
##   res <- OrganismDbi:::.select(x, keys, cls, keytype)

##   checkTrue(dim(res)[2]==11)
##   checkTrue(class(res)=="data.frame")
##   checkTrue("GO" %in% colnames(res)) 
##   checkTrue("ENTREZID" %in% colnames(res))
##   checkTrue("ACCNUM" %in% colnames(res))  
##   checkTrue("ALIAS" %in% colnames(res))
 
## }









## library(OrganismDbi); OrganismDbi:::.test()

