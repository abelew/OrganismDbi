## TODO:
## 1) outline and export (etc.) the methods.
## 2) hook them up to the proper subcomponents (x$txdb)
## 3) need a helper to make the extraction of the txdb SAFE and warn if this is not possible.

.safelyGetTxDb <- function(x){
    if("TXID" %in% cols(x)){
        return(.getTxDb(x))
    }else{
        stop("The OrganismDb object does not have an embedded TranscriptDb.")
    }
}

## with the exception of seqnames (which doesn't make sense, this
## whole family will work for OrganismDb objects (even though I have
## not defined them all explicitely (because they all call seqinfo()
## and have ANY methods)

## So I only need to alias seqinfo and seqnameStyle in the manual page
## (and then put in a see also reference)


#############################################################################

setMethod("seqinfo","OrganismDb", 
          function(x){
              txdb <- .safelyGetTxDb(x)
              seqinfo(txdb)		
})

setReplaceMethod("seqinfo", "OrganismDb",
          function(x, new2old=NULL, force=FALSE, value){
              txdb <- .safelyGetTxDb(x)
              seqinfo(txdb, new2old=NULL, force=FALSE) <- value	
})


## #############################################################################

## setMethod("seqlevels","OrganismDb", 
##           function(x){ seqlevels(.safelyGetTxDb(x))		
## })

## setReplaceMethod("seqlevels", "OrganismDb",
##                  function(x, force=FALSE, value){
##                      seqlevels(.safelyGetTxDb(x), force=FALSE) <- value	
## })

## #############################################################################

## setMethod("seqlengths","OrganismDb", 
##           function(x){ seqlengths(.safelyGetTxDb(x))		
## })

## setReplaceMethod("seqlengths", "OrganismDb",
##                  function(x, value){
##                      seqlengths(.safelyGetTxDb(x)) <- value		
## })

## #############################################################################

## setMethod("isCircular","OrganismDb", 
##           function(x){ isCircular(.safelyGetTxDb(x))		
## })

## setReplaceMethod("isCircular", "OrganismDb",
##                  function(x, value){
##                      isCircular(.safelyGetTxDb(x)) <- value		
## })

## #############################################################################

## setMethod("genome","OrganismDb", 
##           function(x){ genome(.safelyGetTxDb(x))		
## })

## setReplaceMethod("genome", "OrganismDb",
##                  function(x, value){
##                      genome(.safelyGetTxDb(x)) <- value		
## })

#############################################################################

setMethod("seqnameStyle","OrganismDb", 
          function(x){
              txdb <- .safelyGetTxDb(x)
              seqnameStyle(txdb)	
})

setReplaceMethod("seqnameStyle", "OrganismDb",
                 function(x, value){
                     txdb <- .safelyGetTxDb(x)
                     seqnameStyle(txdb) <- value	
})


