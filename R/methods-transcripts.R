## This is where I will put methods to overload things like
## transcripts() and exons()...



## new argument: cols here can be any legit value for cols. (not just
## tx_id and tx_name etc.)

## vals will just pass through to the internal transcripts call.

## columns arg is just for b/c support and will just pass through to
## the internal transcripts call



## How will we merge the results from select() and transcripts()?  We
## will join on tx_id (for transcripts)
.transcripts <- function(x, cols, vals=NULL, columns=c("tx_id", "tx_name")){
    ## 1st get the TranscriptDb object.
    ## trick: there will *always* be a TXID
    txdb <- .lookupDbFromKeytype(x, "TXID")
    
    ## call transcripts method (on the TxDb)
    columns <- unique(c(columns, "tx_id"))
    txs <- transcripts(txdb, vals=vals, columns=columns)
    
    ## call select on the rest and use tx_id as keys 
    meta <- select(x, keys=mcols(txs)$tx_id, cols, "TXID")
    
    ## assemble it all together.
    mcols(txs) <- merge(as(mcols(txs), "data.frame"), meta,
                        by.x="tx_id", by.y="TXID")
    txs
}

setMethod("transcripts", "OrganismDb",
          function(x, cols, vals=NULL, columns=c("tx_id", "tx_name")){
              .transcripts(x, cols, vals=NULL, columns=c("tx_id", "tx_name"))})


## test usage:
## library(Homo.sapiens); h = Homo.sapiens; cols = c("TXNAME","SYMBOL")
## transcripts(h, cols)

