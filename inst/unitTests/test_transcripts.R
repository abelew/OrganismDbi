require("Homo.sapiens")

## TODO: to speed up unit tests I need some small pieces of annotation
## for testing a small subset of granges and the matching metadata
## Ideally I want a Homo.sapiens package that uses the small subset DB
## from GenomicFeatures.  Perhaps a "mini-me" package for testing? -
## but making this is a bit of a project.

x <- Homo.sapiens
txdb <- OrganismDbi:::.getTxDb(x)


## some internal testing (make sure helpers work as expected)

test_compressMetadata <- function(){
    cols <- c("SYMBOL","GENENAME", "CHR", "PMID")
    txs <- transcripts(txdb, vals=NULL, columns="tx_id")[1:100]  ## shortened
    meta <- select(x, keys=mcols(txs)$tx_id, cols, "TXID") 
    f <- factor(meta[["TXID"]],levels=mcols(txs)[["tx_id"]])
    res <- OrganismDbi:::.compressMetadata(f, meta, "TXID")
    checkTrue(class(res)== "DataFrame")
    checkTrue(dim(res)[2] ==4)
    checkTrue(dim(res)[1] ==100)
    checkTrue(all(colnames(res) %in% cols))
}


test_combineMetadata <- function(){
    cols <- c("SYMBOL","GENENAME", "CHR", "PMID")
    txs <- transcripts(txdb, vals=NULL, columns="tx_id")[1:100]  ## shortened
    meta <- select(x, keys=mcols(txs)$tx_id, cols, "TXID") 
    res <- OrganismDbi:::.combineMetadata(txs,meta,avoidID="TXID",
                                          joinID="tx_id")
    checkTrue(class(res)== "DataFrame")
    checkTrue(dim(res)[2] ==4)
    checkTrue(dim(res)[1] ==100)
    checkTrue(all(colnames(res) %in% c(cols))) 
}


## These tests are slow so I will need a smaller thing to test with...
test_transcripts <- function(){
    library(Homo.sapiens); h = Homo.sapiens; cols = c("TXNAME","SYMBOL")
    res <- transcripts(h, columns=cols)
    checkTrue(class(res) == "GRanges")
    checkTrue(length(res) == 80922)
    checkTrue(all(colnames(mcols(res)) %in%
                  c("TXNAME","SYMBOL")))
}

test_exons <- function(){
    library(Homo.sapiens); h = Homo.sapiens; cols = c("CHR","REFSEQ")
    res <- exons(h, columns=cols)
    checkTrue(class(res) == "GRanges")
    checkTrue(length(res) == 286852)
    checkTrue(all(colnames(mcols(res)) %in%
                  c("CHR","REFSEQ")))
}

test_cds <- function(){
    library(Homo.sapiens); h = Homo.sapiens; cols = c("GENENAME","SYMBOL")
    res <- cds(h, columns=cols)
    checkTrue(class(res) == "GRanges")
    checkTrue(length(res) == 235842)
    checkTrue(all(colnames(mcols(res)) %in%
                  c("GENENAME","SYMBOL")))
}





test_transcriptsBy <- function(){
    library(Homo.sapiens);h=Homo.sapiens;by="gene";cols = c("GENENAME","SYMBOL")
    res <- transcriptsBy(h, by="gene", cols)    
    checkTrue(class(res) == "GRangesList")
    checkTrue(length(res) == 22932)
    checkTrue(all(colnames(mcols(res)) %in%
                  c("GENENAME","SYMBOL")))

    ## extra check for case where we only have one field.
    cols = c("SYMBOL")
    res2 <- transcriptsBy(h, by="gene", cols)
    checkTrue(class(res) == "GRangesList")
    checkTrue(length(res) == 22932)
    checkTrue(all(colnames(mcols(res)) %in%
                  c("SYMBOL")))

    ## 
}

test_exonsBy <- function(){
    library(Homo.sapiens);h=Homo.sapiens;by="gene";cols = c("GENENAME","SYMBOL")
    res <- exonsBy(h, by="gene", cols)
    ## TODO: look more closely at this one.  The metadata looks off...
    checkTrue(class(res) == "GRangesList")
    checkTrue(length(res) == 22932) 
    checkTrue(all(colnames(mcols(res)) %in%
                  c("GENENAME","SYMBOL")))
}

test_cdsBy <- function(){
    library(Homo.sapiens);h=Homo.sapiens;by="gene";cols = c("GENENAME","SYMBOL")
    res <- cdsBy(h, by="gene", cols)
    checkTrue(class(res) == "GRangesList")
    checkTrue(length(res) == 19511)
    checkTrue(all(colnames(mcols(res)) %in%
                  c("GENENAME","SYMBOL")))
}






#################
## Another bug.  The following should work? (but doesn't)
## it's probably something to do with expansion of TXIDs by select?
## It's a corner case...  Cause because TXID is on the "avoid" list
## and is also requested here...  The problems caused by this will
## change as I remove the default values for tx_id etc. from the
## internal call to transcriptsBy()...  So that means that 'avoding'
## TXID will also be unecessary.

## THEREFORE: TODO = resolve the other issues 1st.
## x = Homo.sapiens; txby= transcriptsBy(x, by="gene", "TXID")
