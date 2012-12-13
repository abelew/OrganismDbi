## This is where I will put methods to overload things like
## transcripts() and exons()...



## new argument: cols here can be any legit value for cols. (not just
## tx_id and tx_name etc.)

## vals will just pass through to the internal transcripts call.

## columns arg is just for b/c support and will just pass through to
## the internal transcripts call


.getTxDb <- function(x){
    ## trick: there will *always* be a TXID
    .lookupDbFromKeytype(x, "TXID")
}


## This helper processes data.frame data into a DataFrame with compressed chars
.compressMetadata <- function(rngs, meta, avoidID, joinID){
    ## make a special factor
    f <- factor(meta[[avoidID]],levels=mcols(rngs)[[joinID]])
    ## call splitAsList (using factor) on all cols except avoidId
    cols <- meta[,!colnames(meta) %in% avoidID]
    res <- lapply(cols, splitAsList, f) ## fast
    ## call unique on all cols
    res <- lapply(res, unique) 
    resf <- as(res, "DataFrame") ##do.call(DataFrame, res)
    ## attach to mcols values. from before.
    if(dim(mcols(rngs))[1] == dim(resf)[1]){
        return(c(mcols(rngs),resf))
    }else{
        stop("Ranges and annotations retrieved are not of matching lengths.")
    }
}



## How will we merge the results from select() and transcripts()?  We
## will join on tx_id (for transcripts)
.transcripts <- function(x, cols, vals=NULL, columns=c("tx_id", "tx_name")){
    ## 1st get the TranscriptDb object.
    txdb <- .getTxDb(x)
    ## call transcripts method (on the TxDb)
    columns <- unique(c(columns, "tx_id")) ## tx_id always exists
    txs <- transcripts(txdb, vals=vals, columns=columns)  
    ## call select on the rest and use tx_id as keys 
    meta <- select(x, keys=mcols(txs)$tx_id, cols, "TXID")    
    ## assemble it all together.
    mcols(txs) <- .compressMetadata(txs,meta,avoidID="TXID",joinID="tx_id") 
    txs
}

setMethod("transcripts", "OrganismDb",
          function(x, cols, vals=NULL, columns=c("tx_id", "tx_name")){
              .transcripts(x, cols, vals=NULL, columns=c("tx_id", "tx_name"))})


## test usage:
## library(Homo.sapiens); h = Homo.sapiens; cols = c("TXNAME","SYMBOL")
## transcripts(h, cols)


## How will we merge the results from select() and transcripts()?  We
## will join on tx_id (for transcripts)
.exons <- function(x, cols, vals=NULL, columns="exon_id"){
    ## 1st get the TranscriptDb object.
    txdb <- .getTxDb(x)
    
    ## call transcripts method (on the TxDb)
    columns <- unique(c(columns, "exon_id")) ## exon_id always exists
    exs <- exons(txdb, vals=vals, columns=columns)
    
    ## call select on the rest and use tx_id as keys 
    meta <- select(x, keys=mcols(exs)$exon_id, cols, "EXONID")
    
    ## assemble it all together.
    mcols(exs) <- .compressMetadata(exs,meta,avoidID="EXONID",joinID="exon_id")
    exs
}

setMethod("exons", "OrganismDb",
          function(x, cols, vals=NULL, columns="exon_id"){
              .exons(x, cols, vals=NULL, columns="exon_id")})


## test usage:
## library(Homo.sapiens); h = Homo.sapiens; cols = c("CHR","SYMBOL")
## exons(h, cols)


## How will we merge the results from select() and transcripts()?  We
## will join on tx_id (for transcripts)
.cds <- function(x, cols, vals=NULL, columns="cds_id"){
    ## 1st get the TranscriptDb object.
    txdb <- .getTxDb(x)
    
    ## call transcripts method (on the TxDb)
    columns <- unique(c(columns, "cds_id")) ## cds_id always exists
    cds <- cds(txdb, vals=vals, columns=columns)
    
    ## call select on the rest and use tx_id as keys 
    meta <- select(x, keys=mcols(cds)$cds_id, cols, "CDSID")
    
    ## assemble it all together.
    mcols(cds) <- .compressMetadata(cds,meta,avoidID="CDSID",joinID="cds_id")
    cds
}

setMethod("cds", "OrganismDb",
          function(x, cols, vals=NULL, columns="cds_id"){
              .cds(x, cols, vals=NULL, columns="cds_id")})


## test usage:
## library(Homo.sapiens); h = Homo.sapiens; cols = c("GENENAME","SYMBOL")
## cds(h, cols)




########################################################################
## General problem: I will usually have more stuff to cram into mcols
## than I have rows of ranges...  How should we handle this?


## You can see this in action by doing this (for example)
## library(Homo.sapiens); h = Homo.sapiens; cols = c("TXNAME","SYMBOL")
## exonsh, cols)


## If I am lucky, there will be a method already to squish a
## DataFrame() and so I will be able to just call that in a helper
## when I am merging...
