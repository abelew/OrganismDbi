## TODO:
## 1) outline and export (etc.) the methods.
## 2) hook them up to the proper subcomponents (x$txdb)
## 3) need a helper to make the extraction of the txdb SAFE and warn if this is not possible.

.safelyGetTxDb <- function(x){
    if("TXID" %in% columns(x)){
        return(.getTxDb(x))
    }else{
        stop("The MultiDb object does not have an embedded TxDb.")
    }
}

## with the exception of seqnames (which doesn't make sense, this
## whole family will work for MultiDb objects (even though I have
## not defined them all explicitely (because they all call seqinfo()
## and have ANY methods)

## So I only need to alias seqinfo and seqnameStyle in the manual page
## (and then put in a see also reference)


#############################################################################

## To get all these functions, I ONLY need to make this setter and getter work

setMethod("seqinfo","MultiDb", 
          function(x){
              txdb <- .safelyGetTxDb(x)
              seqinfo(txdb)		
})


## This can work once I have a local object that I can modify?

## setReplaceMethod("seqinfo", "MultiDb",
##           function(x, new2old=NULL, force=FALSE, value){
##               txdb <- .safelyGetTxDb(x)
##               seqinfo(txdb, new2old=NULL, force=FALSE) <- value	
## })



## #############################################################################

## setMethod("seqlevels","MultiDb", 
##           function(x){ seqlevels(.safelyGetTxDb(x))		
## })

## setReplaceMethod("seqlevels", "MultiDb",
##                  function(x, force=FALSE, value){
##                      seqlevels(.safelyGetTxDb(x), force=FALSE) <- value	
## })

## #############################################################################

## setMethod("seqlengths","MultiDb", 
##           function(x){ seqlengths(.safelyGetTxDb(x))		
## })

## setReplaceMethod("seqlengths", "MultiDb",
##                  function(x, value){
##                      seqlengths(.safelyGetTxDb(x)) <- value		
## })

## #############################################################################

## setMethod("isCircular","MultiDb", 
##           function(x){ isCircular(.safelyGetTxDb(x))		
## })

## setReplaceMethod("isCircular", "MultiDb",
##                  function(x, value){
##                      isCircular(.safelyGetTxDb(x)) <- value		
## })

## #############################################################################

## setMethod("genome","MultiDb", 
##           function(x){ genome(.safelyGetTxDb(x))		
## })

## setReplaceMethod("genome", "MultiDb",
##                  function(x, value){
##                      genome(.safelyGetTxDb(x)) <- value		
## })

#############################################################################

## setMethod("seqnameStyle","MultiDb", 
##           function(x){
##               txdb <- .safelyGetTxDb(x)
##               seqnameStyle(txdb)	
## })

## setReplaceMethod("seqnameStyle", "MultiDb",
##                  function(x, value){
##                      txdb <- .safelyGetTxDb(x)
##                      seqnameStyle(txdb) <- value
##                      x@
## })


