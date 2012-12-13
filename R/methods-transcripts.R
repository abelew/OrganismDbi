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

## TODO: .compressMetadata() might be useful to move into IRanges, as
## a complement to expand() methods?  - Discuss this with Val (who
## apparently may have similar issues in vcf...

## .compressMetadata() processes data.frame data into a DataFrame with
## compressed chars

## It does so by taking a special factor (f) and then applying it to
## ALL of the columns in a data.frame (meta) except for the one that
## was the basis for the special factor (avoidID)
.compressMetadata <- function(f, meta, avoidID){
    cols <- meta[,!colnames(meta) %in% avoidID]
    ## call splitAsList (using factor) on all cols except avoidId
    res <- lapply(cols, splitAsList, f) ## fast
    ## call unique on all cols
    res <- lapply(res, unique)  ## slower
    as(res, "DataFrame") 
}

## This helper does book keeping that is relevant to my situation here.
.combineMetadata <- function(rngs, meta, avoidID, joinID){
    ## make a special factor
    f <- factor(meta[[avoidID]],levels=mcols(rngs)[[joinID]])
    ## compress the metadata by splitting according to f
    res <- .compressMetadata(f, meta, avoidID)
    ## attach to mcols values. from before.
    if(dim(mcols(rngs))[1] == dim(res)[1]){
        return(c(mcols(rngs),res))
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
    mcols(txs) <- .combineMetadata(txs,meta,avoidID="TXID",joinID="tx_id") 
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
    mcols(exs) <- .combineMetadata(exs,meta,avoidID="EXONID",joinID="exon_id")
    exs
}

setMethod("exons", "OrganismDb",
          function(x, cols, vals=NULL, columns="exon_id"){
              .exons(x, cols, vals=NULL, columns="exon_id")})


## test usage:
## library(Homo.sapiens); h = Homo.sapiens; cols = c("CHR","REFSEQ")
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
    mcols(cds) <- .combineMetadata(cds,meta,avoidID="CDSID",joinID="cds_id")
    cds
}

setMethod("cds", "OrganismDb",
          function(x, cols, vals=NULL, columns="cds_id"){
              .cds(x, cols, vals=NULL, columns="cds_id")})


## test usage:
## library(Homo.sapiens); h = Homo.sapiens; cols = c("GENENAME","SYMBOL")
## cds(h, cols)




########################################################################
########################################################################
##                       The "By" methods
########################################################################
########################################################################

## "By" methods will just cram the same metadata into the INTERNAL
## metadata slot so that it appears with the show method.
## No attempt will be made to manage the insanity of knowing which
## metadata types belong in which spot...


.transcriptsBy <- function(x, by, cols){
    ## 1st get the TranscriptDb object.
    txdb <- .getTxDb(x)
    ## call transcriptsBy with use.names set to FALSE
    txby <- transcriptsBy(txdb, by=by, use.names=FALSE)

    ## get the tx_ids from the transcripts
    ## AND I need to one from the internal slot.
    gr <- txby@unlistData
    k  <- mcols(gr)$tx_id
    
    ## call select on the rest and use tx_id as keys 
    meta <- select(x, keys=k, cols, "TXID")    
    ## assemble it all together.
    mcols(gr) <- .combineMetadata(gr, meta, avoidID="TXID", joinID="tx_id") 

    ## now cram it back in there.
    txby@unlistData <- gr
    txby
}

setMethod("transcriptsBy", "OrganismDb",
          function(x, by, cols){
              .transcriptsBy(x, by, cols)})



## library(Homo.sapiens);h=Homo.sapiens;by="gene";cols = c("GENENAME","SYMBOL")
## transcriptsBy(h, by="gene", cols)







.exonsBy <- function(x, by, cols){
    ## 1st get the TranscriptDb object.
    txdb <- .getTxDb(x)
    ## call transcriptsBy with use.names set to FALSE
    exby <- exonsBy(txdb, by=by, use.names=FALSE)

    ## get the tx_ids from the transcripts
    ## AND I need to one from the internal slot.
    gr <- exby@unlistData
    k  <- mcols(gr)$exon_id
    
    ## call select on the rest and use tx_id as keys 
    meta <- select(x, keys=k, cols, "EXONID")    
    ## assemble it all together.
    mcols(gr) <- .combineMetadata(gr, meta, avoidID="EXONID", joinID="exon_id") 

    ## now cram it back in there.
    exby@unlistData <- gr
    exby
}

setMethod("exonsBy", "OrganismDb",
          function(x, by, cols){
              .exonsBy(x, by, cols)})



## library(Homo.sapiens);h=Homo.sapiens;by="gene";cols = c("GENENAME","SYMBOL")
## exonsBy(h, by="tx", cols)





.cdsBy <- function(x, by, cols){
    ## 1st get the TranscriptDb object.
    txdb <- .getTxDb(x)
    ## call transcriptsBy with use.names set to FALSE
    cdsby <- cdsBy(txdb, by=by, use.names=FALSE)

    ## get the tx_ids from the transcripts
    ## AND I need to one from the internal slot.
    gr <- cdsby@unlistData
    k  <- mcols(gr)$cds_id
    
    ## call select on the rest and use tx_id as keys 
    meta <- select(x, keys=k, cols, "CDSID")    
    ## assemble it all together.
    mcols(gr) <- .combineMetadata(gr, meta, avoidID="CDSID", joinID="cds_id") 

    ## now cram it back in there.
    cdsby@unlistData <- gr
    cdsby
}

setMethod("cdsBy", "OrganismDb",
          function(x, by, cols){
              .cdsBy(x, by, cols)})



## library(Homo.sapiens);h=Homo.sapiens;by="gene";cols = c("GENENAME","SYMBOL")
## cdsBy(h, by="tx", cols)





## TODO: (known issues)
## 1) cols don't come back in same order that the went in

## 2) some values (tx_id and tx_name come to mind) are not relabeled
## in a pretty way and may not have been requested (to solve this we
## have to adress issue #3)

## 3) I now have a columns AND a cols argument for the transcripts()
## family of methods.  This is totally redundant.  Proposed fix:
## rename arguments base method to be cols (maybe this is also an
## opportunity to rename cols everywhere), but rename it so that it's
## consistent, and then here, just only have one argument...

## 4) exonsBy and cdsBy may have some extra issues that I am missing...

