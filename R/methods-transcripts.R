## This is where I will put methods to overload things like
## transcripts() and exons()...

## new argument: columns here can be any legit value for columns. (not just
## tx_id and tx_name etc.)

## vals will just pass through to the internal transcripts call.

## columns arg is just for b/c support and will just pass through to
## the internal transcripts call


## For consistency, the helper columns just wraps around cols method...
setMethod("columns", "OrganismDb", cols)


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
    columns <- meta[,!colnames(meta) %in% avoidID, drop=FALSE]
    ## call splitAsList (using factor) on all columns except avoidId
    res <- lapply(columns, splitAsList, f) ## fast 
    ## call unique on all columns
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
        res <- c(mcols(rngs),res)
        ## throw out joining IDs
        res <- res[!(colnames(res) %in% c("tx_id","exon_id","cds_id"))]
        return(res)
    }else{
        stop("Ranges and annotations retrieved are not of matching lengths.")
    }
}



## How will we merge the results from select() and transcripts()?  We
## will join on tx_id (for transcripts)
.transcripts <- function(x, vals, columns){
    ## 1st get the TranscriptDb object.
    txdb <- .getTxDb(x)
    ## call transcripts method (on the TxDb)
    txs <- transcripts(txdb, vals, columns="tx_id")  
    ## call select on the rest and use tx_id as keys 
    meta <- select(x, keys=mcols(txs)$tx_id, columns, "TXID")    
    ## assemble it all together.
    mcols(txs) <- .combineMetadata(txs,meta,avoidID="TXID",joinID="tx_id") 
    txs
}

setMethod("transcripts", "OrganismDb",
          function(x, vals=NULL, columns=c("TXID", "TXNAME")){
              .transcripts(x, vals, columns)})


## test usage:
## library(Homo.sapiens); h = Homo.sapiens; columns = c("TXNAME","SYMBOL")
## transcripts(h, columns)


## How will we merge the results from select() and transcripts()?  We
## will join on tx_id (for transcripts)
.exons <- function(x, vals, columns){
    ## 1st get the TranscriptDb object.
    txdb <- .getTxDb(x)
    
    ## call transcripts method (on the TxDb)
    exs <- exons(txdb, vals, columns="exon_id")
    
    ## call select on the rest and use tx_id as keys 
    meta <- select(x, keys=mcols(exs)$exon_id, columns, "EXONID")
    
    ## assemble it all together.
    mcols(exs) <- .combineMetadata(exs,meta,avoidID="EXONID",joinID="exon_id")
    exs
}

setMethod("exons", "OrganismDb",
          function(x, vals=NULL, columns="EXONID"){
              .exons(x, vals, columns)})


## test usage:
## library(Homo.sapiens); h = Homo.sapiens; columns = c("CHR","REFSEQ")
## exons(h, columns)


## How will we merge the results from select() and transcripts()?  We
## will join on tx_id (for transcripts)
.cds <- function(x, vals, columns){
    ## 1st get the TranscriptDb object.
    txdb <- .getTxDb(x)
    
    ## call transcripts method (on the TxDb)
    cds <- cds(txdb, vals, columns="cds_id")
    
    ## call select on the rest and use tx_id as keys 
    meta <- select(x, keys=mcols(cds)$cds_id, columns, "CDSID")
    
    ## assemble it all together.
    mcols(cds) <- .combineMetadata(cds,meta,avoidID="CDSID",joinID="cds_id")
    cds
}

setMethod("cds", "OrganismDb",
          function(x, vals=NULL, columns="CDSID"){
              .cds(x, vals, columns)})


## test usage:
## library(Homo.sapiens); h = Homo.sapiens; columns = c("GENENAME","SYMBOL")
## cds(h, columns)




########################################################################
########################################################################
##                       The "By" methods
########################################################################
########################################################################

## "By" methods will just cram the same metadata into the INTERNAL
## metadata slot so that it appears with the show method.
## No attempt will be made to manage the insanity of knowing which
## metadata types belong in which spot...


.transcriptsBy <- function(x, by, columns){
    ## 1st get the TranscriptDb object.
    txdb <- .getTxDb(x)
    ## call transcriptsBy with use.names set to FALSE
    txby <- transcriptsBy(txdb, by=by, use.names=FALSE)

    ## get the tx_ids from the transcripts
    ## AND I need to one from the internal slot.
    gr <- txby@unlistData
    k  <- mcols(gr)$tx_id
    
    ## call select on the rest and use tx_id as keys 
    meta <- select(x, keys=k, columns, "TXID")    
    ## assemble it all together.
    mcols(gr) <- .combineMetadata(gr, meta, avoidID="TXID", joinID="tx_id") 

    ## now cram it back in there.
    txby@unlistData <- gr
    txby
}

setMethod("transcriptsBy", "OrganismDb",
          function(x, by, columns){
              .transcriptsBy(x, by, columns)})



## library(Homo.sapiens);h=Homo.sapiens;by="gene";columns = c("GENENAME","SYMBOL")
## transcriptsBy(h, by="gene", columns)







.exonsBy <- function(x, by, columns){
    ## 1st get the TranscriptDb object.
    txdb <- .getTxDb(x)
    ## call transcriptsBy with use.names set to FALSE
    exby <- exonsBy(txdb, by=by, use.names=FALSE)

    ## get the tx_ids from the transcripts
    ## AND I need to one from the internal slot.
    gr <- exby@unlistData
    k  <- mcols(gr)$exon_id
    
    ## call select on the rest and use tx_id as keys 
    meta <- select(x, keys=k, columns, "EXONID")    
    ## assemble it all together.
    mcols(gr) <- .combineMetadata(gr, meta, avoidID="EXONID", joinID="exon_id") 

    ## now cram it back in there.
    exby@unlistData <- gr
    exby
}

setMethod("exonsBy", "OrganismDb",
          function(x, by, columns){
              .exonsBy(x, by, columns)})



## library(Homo.sapiens);h=Homo.sapiens;by="gene";columns = c("GENENAME","SYMBOL")
## exonsBy(h, by="tx", columns)





.cdsBy <- function(x, by, columns){
    ## 1st get the TranscriptDb object.
    txdb <- .getTxDb(x)
    ## call transcriptsBy with use.names set to FALSE
    cdsby <- cdsBy(txdb, by=by, use.names=FALSE)

    ## get the tx_ids from the transcripts
    ## AND I need to one from the internal slot.
    gr <- cdsby@unlistData
    k  <- mcols(gr)$cds_id
    
    ## call select on the rest and use tx_id as keys 
    meta <- select(x, keys=k, columns, "CDSID")    
    ## assemble it all together.
    mcols(gr) <- .combineMetadata(gr, meta, avoidID="CDSID", joinID="cds_id") 

    ## now cram it back in there.
    cdsby@unlistData <- gr
    cdsby
}

setMethod("cdsBy", "OrganismDb",
          function(x, by, columns){
              .cdsBy(x, by, columns)})



## library(Homo.sapiens);h=Homo.sapiens;by="gene";columns = c("GENENAME","SYMBOL")
## cdsBy(h, by="tx", columns)





## TODO: (known issues)
## 1) columns don't come back in same order that the went in

## 2) some values (tx_id and tx_name come to mind) are not relabeled
## in a pretty way and may not have been requested (to solve this we
## have to adress issue #3)

## 3) I now have a columns AND a columns argument for the transcripts()
## family of methods.  This is totally redundant.  Proposed fix:
## rename arguments base method to be columns (maybe this is also an
## opportunity to rename columns everywhere), but rename it so that it's
## consistent, and then here, just only have one argument...

## 4) exonsBy and cdsBy may have some extra issues that I am missing...
