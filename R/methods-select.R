# This will just hold code for the initial implementation of select and friends

## helper to convert text strings (Db pkgs names) into real objects
.makeReal <- function(x){
    eval(parse(text=x))
}

## Standard methods:
.keytypes <- function(x){
    dbs <- .getDbObjs(x)
    unique(unlist(lapply(dbs, keytypes)))
}

setMethod("keytypes", "OrganismDb", .keytypes)

## Usage:
## keytypes(Homo.sapiens)

.cols <- function(x){
    dbs <- .getDbObjs(x)
    unique(unlist(lapply(dbs, columns)))
}

setMethod("columns", "OrganismDb", function(x){.cols(x)})

## Usage:
## columns(Homo.sapiens)

## Strategy for keys: I need a lookup function that can 1) generate the keys
## for each slot and then lookup which slot I should be tapping based on a
## keytype.  2) This needs to be general purpose (will be needed again in
## select) and 3) it also may need to be able to return multiple hits in the
## event that there are eventually multiple IDs named the same way (depends on
## whether or not we allow repeat ID names).  I think we WILL want to allow
## this, which means I will have to do some kind of name-spacing scheme.


.makekeytypeMapping <- function(x){
    objs <- .getDbObjs(x)
    unlist2(lapply(objs, keytypes))
}

.lookupDbNameFromKeytype <- function(x, keytype){
    res <- .makekeytypeMapping(x)
    ## no duplicates so I can just return the name
    names(res)[res %in% keytype]  
}

.lookupDbFromKeytype <- function(x, keytype){
    db <- .lookupDbNameFromKeytype(x, keytype)
    eval(parse(text=db))
}

.keys <- function(x, keytype, ...){
    AnnotationDbi:::.testForValidKeytype(x, keytype)
    db <- .lookupDbFromKeytype(x, keytype)
    ## And then we can just call keys...
    as.character(keys(db, keytype, ...))
}

setMethod("keys", "OrganismDb", .keys)


## Usage: 
## head(keys(Homo.sapiens, keytype="PMID"))
## the use case for GOID will present a special challenge...
## head(keys(Homo.sapiens, keytype="GOID"))

## This method just gets me the pkg names as names and vals are fkeys
.getDbNameFKeys <- function(x){
    gd <- keyFrame(x)
    ## now give all the keys as a vector, but named by their databases.
    .extractPkgsAndCols(gd)
}


## .mkeys will return appropriate value "on the fly" based on the
## contents of keyFrame().  It will take at least three arguments: the two
## tables plus an indicator for which of the two keys 1st or 2nd table key is
## needed.

## tbl1,tbl2 wil be actual package names like 'org.Hs.eg.db' or 'GO.db'
## tbl1 = "TxDb.Hsapiens.UCSC.hg19.knownGene"
## tbl2 = "org.Hs.eg.db"
## key = "tbl1"
.parseCol <- function(piece, str) grepl(str, piece)

.mkeys <- function(x, tbl1, tbl2, key=c("tbl1","tbl2", "both")){
    if(length(tbl1) != 1L || length(tbl2) != 1L)
        stop("specify only one pair of tables at a time")
    key <- match.arg(key)
    kf <- keyFrame(x)
    ## process for a double match of tbl1 and tbl2 (in any order)
    ## note: (we should ALWAYS have one when this function is called)
    
    res <- apply(kf[,1:2], MARGIN=2, FUN=.parseCol, tbl1)
    res2 <- apply(kf[,1:2], MARGIN=2, FUN=.parseCol, tbl2)
    fin <- res | res2 
    resRowIdx <- fin[,1] & fin[,2]
    matchRow <- kf[resRowIdx,]
    if(length(matchRow) == 0L)
        stop("no relationship found for ",tbl1," and ",tbl2)
    
    ## now the tricky part is that in returning the keys I have to get the
    ## correct keys back to the user...  And this is based on whether tbl1 was
    ## one thing or another.
    if(length(matchRow[["xDbs"]]) >1L)
        stop("failed to limit choices to 1")
    if(key=="tbl1"){
        if(grepl(tbl1,matchRow[["xDbs"]])){
            ans <- as.character(matchRow[["xKeys"]])
        }else{ ## then its reversed of the order in the row...
            ans <- as.character(matchRow[["yKeys"]])
        }
    }else if(key=="tbl2"){
        if(grepl(tbl2,matchRow[["yDbs"]])){ 
            ans <- as.character(matchRow[["yKeys"]])
        }else{ ## and the reverse case
            ans <- as.character(matchRow[["xKeys"]])
        }
    }else if(key=="both"){
        ans <- c(as.character(matchRow[["xKeys"]]),
                 as.character(matchRow[["yKeys"]]))
        names(ans) <- c(as.character(matchRow[["xDbs"]]),
                        as.character(matchRow[["yDbs"]]))
        ## When we say "both" we still want keys returned in same order as
        ## original packages.  IOW, if tbl1 goes with key 1, then we should list
        ## key 1 1st in the result...
        ans <- ans[match(c(tbl1,tbl2),names(ans))]
    }
    ans
}


## helper for getting all cols by all nodes
.colsByNodes <- function(x){
    gr <- dbGraph(x)
    allCols <- lapply(nodes(gr), function(elt) columns(.makeReal(elt)))
    names(allCols) <- nodes(gr)
    allCols
}
## library(Homo.sapiens)
## library(RBGL)
## library(graph)a
## x = Homo.sapiens
## allCols <- .colsByNodes(x)

## helper to get the subgraph
.getRelevantSubgraph <- function(x, cols, keys, keytype){
    gr <- dbGraph(x)
    allCols <- .colsByNodes(x)
    inSubgraph <- sapply(allCols,
           function(cols, keys) any(keys %in% cols),  union(keytype, cols))
    subGraph(names(inSubgraph)[inSubgraph], gr)
}
## kt <- "ENTREZID"
## cls = c("GOID" ,  "SYMBOL", "TXNAME")
## keys <- head(keys(x, "ENTREZID"))
## subgr <- .getRelevantSubgraph(x, cols=cls, keys, keytype=kt)
## We will also need the root
## root = OrganismDbi:::.lookupDbNameFromKeytype(x, kt)


## I think this is meant to be an lapply
.getForeignKeys <- function(x, subgr){

    fKeys <- lapply(strsplit(edgeNames(subgr), "~"),
                    function(tables, x, key)
                    .mkeys(x, tables[[1]], tables[[2]], "both"),
                    x)

    unlist(fKeys, use.names=FALSE)
}
## fKeys <- .getForeignKeys(x, subgr)


## now combine all the keys together
## selectCols = unique(c(kt, fKeys, cls))

## sort the needed cols by their nodes
.getColsByNodes <- function(subgr, selectCols, allCols){
    lapply(allCols[nodes(subgr)],
           function(col, selectCols) col[col %in% selectCols], selectCols)
}
## needCols <- .getColsByNodes(subgr, selectCols, allCols)


## get list of nodes to visit
.bfs <- function(object, node)
    ## names are bfs order; values are 'from' nodes
{
    bfs <- bfs(object, node)
    from <- sapply(edges(object)[bfs], function(table, x) {
        x[which.max(x %in% table)]
    }, bfs)
    from[1] <- NA
    from
}
## So our visitNodes then becomes:
## visitNodes = .bfs(subgr, root)



## new version of .getSelects()
## ## select 
.getSelects <- function(x, keytype, keys, needCols, visitNodes, fks=NULL){
    ## set up an empty list with names that match what we want to fill...
    selected <- setNames(
                         vector("list", length(visitNodes)),
                         names(visitNodes))
    ## in 1st case we only need the name
    node1 <- names(visitNodes)[[1]]
    selected[[node1]] <- 
        select(.makeReal(node1),
               keys=as.character(keys),
               columns=needCols[[node1]],
               keytype=keytype,
               fks=fks) ## will usually be NULL
    ## but here we need to use the name and the value of visitNodes
    otherNodes <- visitNodes[-1] 
    for (i in seq_len(length(otherNodes))) {
        nodeName <- names(otherNodes)[i]
        fromNode <- otherNodes[i] 
        fromKey <- .mkeys(x, fromNode, nodeName, "tbl1")
        fromKeys <- unique(selected[[fromNode]][[fromKey]])
        fromKeys <- fromKeys[!is.na(fromKeys)]
        toKey <- .mkeys(x, fromNode, nodeName, "tbl2")
        selected[[nodeName]] <- 
            select(.makeReal(nodeName),
                   keys=as.character(fromKeys),
                   columns=needCols[[nodeName]],
                   keytype=toKey,
                   fks=fks)
        ## TODO: see fileWithProblemsForValidityKeys.R for more details...
    }
    selected
}
## selected <- .getSelect(kt,keys,needCols, visitNodes)


## new version of .mergeSelectResults
## merge
.mergeSelectResults <- function(x, selected, visitNodes, oriCols){
    final <- selected[[1]]
    otherNodes <- visitNodes[-1] 
    for (i in seq_len(length(otherNodes))) {
        nodeName <- names(otherNodes)[i]
        fromNode <- otherNodes[i] 
        fromKey <- .mkeys(x, fromNode, nodeName, "tbl1")
        toKey <- .mkeys(x, fromNode, nodeName, "tbl2")
        final <- merge(final, selected[[nodeName]],
                       by.x=fromKey, by.y=toKey, all=TRUE)
        ## recover the col that is lost from the merge
        ## (header is sometimes needed)
        lostKeys <- data.frame(toKey=final[[fromKey]])
        colnames(lostKeys) <- toKey
        final <- cbind(final, lostKeys) ## bind b.c lostKeys is post-merge clone 
    }
    final
}
## res <- .mergeSelectResults(selected, visitNodes)

## helper to get fks when it's needed.  This returns NULL when not
## appropriate and a compound list of keys whenever it is.
.getFksWhenAppropriate <- function(x, keytype){
    forgnKeys<-as.data.frame(x@keys,stringsAsFactors=FALSE)[,c('xKeys','yKeys')]
    keysIdx<-grepl(keytype, forgnKeys)
    if(any(keysIdx)){
        getKeys <- as.character(forgnKeys[keysIdx,]) ## will always be one row
        res <- unique(c(keys(x, keytype=getKeys[1]),
                        keys(x, keytype=getKeys[2])))
    }else{ res <- NULL}
 res   
}

.select <- function(x, keys, cols, keytype, ...){
    ## Argument checking:
    if(missing(keys)){stop("You must provide a keys argument")}
    if(missing(cols)){stop("You must provide columns argument")}
    if(missing(keytype)){stop("You must provide a keytype argument")}
    ## Some more argument checking
##    fks <- .getFksWhenAppropriate(x, keytype)
    fks <- NULL ## This returns us (temporarily) to the older bahavior...
    AnnotationDbi:::.testSelectArgs(x, keys=keys, cols=cols, keytype=keytype,
                                    fks)
    ## if asked for what they have, just return that.
    if(all(cols %in% keytype)  && length(cols)==1L){
        res <- data.frame(keys=keys)
        colnames(res) <- cols
        return(res)
    }
    
    ## Preserve original cols (we will be adding some to get our results
    ## along the way
    oriCols <- cols  
    
    ## New methods make more use of graph objects.
    allCols <- .colsByNodes(x)
    subgr <- .getRelevantSubgraph(x, cols=cols, keys, keytype=keytype)
    root <- .lookupDbNameFromKeytype(x, keytype)
    fKeys <- .getForeignKeys(x, subgr)
    selectCols <- unique(c(keytype, fKeys, cols))
    needCols <- .getColsByNodes(subgr, selectCols, allCols)
    visitNodes <- .bfs(subgr, root)
    selected <- .getSelects(x, keytype,keys,needCols, visitNodes, fks=fks)
    res <- .mergeSelectResults(x, selected, visitNodes, oriCols)
    
    ## Next we need to filter out all columns that we didn't ask for.  
    ## Actually that is not quite right, what we want to do is make a blacklist
    ## of columns that were added (in fkeys) and that were NOT requested
    ## (oriCols and keytype).
    
    extraKeys <- .getDbNameFKeys(x)
    blackList <- extraKeys[!(extraKeys %in% unique(c(oriCols, keytype)))]
    ## if they asked for one of the GO items, then GO is not blacklisted
    ##   if(any(columns(GO.db) %in% oriCols)){
    ##     blackList <- blackList[!(blackList %in% "GO")]
    ##   }
    res <- res[,!(colnames(res) %in% blackList), drop=FALSE] 
    
    ## Then call code to clean up, reorder the rows (and add NA rows as needed).
    if(nrow(res) > 0L){
        res <- AnnotationDbi:::.resort(tab=res, keys=keys,
                                       jointype=keytype,
                                       reqCols=colnames(res))
    }
#    unique(res) ## NO! We don't want to do this.
    res
}


## ##  Remove this select warning function after 2.13 has released
## .selectWarnOrganismDb <- function(x, keys, columns, keytype, ...){
##     extraArgs <- list(...)
##     if("cols" %in% names(extraArgs)){ 
##         ## warn the user about the old argument
##         AnnotationDbi:::.colsArgumentWarning()
##         ## then call it using cols in place of columns
##         ## YES keytype=columns.  Really!, but ONLY when keytype is null...
##         ## if(missing(keytype)){
##         ##     .select(x, keys, extraArgs[["cols"]], keytype = columns, ...)
##         ## }else{
##         ##     .select(x, keys, extraArgs[["cols"]], keytype = keytype, ...)
##         ## }
##     }else{
##         .select(x, keys, columns, keytype, ...)
##     }
## }

setMethod("select", "OrganismDb",
          function(x, keys, columns, keytype, ...){
            ## .selectWarnOrganismDb(x, keys, columns, keytype, ...)
            .select(x, keys, columns, keytype, ...)
          }
)


##TODO: .mergeSelectResults is leaving incorrect labels on things:  Clean this up!


## methods for easy DB access:
.dbconn <- function(x){
    dbs <- .getDbObjs(x)
    res <- unique(unlist(lapply(dbs, dbconn)))
    names(res) <- names(dbs)
    res
}
setMethod("dbconn", "OrganismDb", function(x){.dbconn(x)})

.dbfile <- function(x){
    dbs <- .getDbObjs(x)
    res <- unique(unlist(lapply(dbs, dbfile)))
    names(res) <- names(dbs)
    res
}
setMethod("dbfile", "OrganismDb", function(x){.dbfile(x)})


## mapIds
## Standard methods:
.mapIds <- function(x, keys, column, keytype, ...,
                    multiVals=c("filter","asNA","first","list",
                      "CharacterList")){
    AnnotationDbi:::.testForValidKeytype(x, keytype)
    db <- .lookupDbFromKeytype(x, keytype)
    mapIds(db, keys, column, keytype, ...,
                    multiVals=multiVals)
}

setMethod("mapIds", "OrganismDb", .mapIds)

## library(Homo.sapiens);  mapIds(Homo.sapiens, keys=c('1','10'), column='ALIAS', keytype='ENTREZID',  multiVals="CharacterList")
## TODO: add some unit tests for this.


## taxonomyId for OrganismDb relies on the TxDb object.
## this could be changed to instead check the OrgDb object.
## but that would require adding a new helper "getOrgDbIfAvailable()"
.taxonomyId <- function(x){
    txdb <- getTxDbIfAvailable(x)
    taxonomyId(txdb)
}
setMethod("taxonomyId", "OrganismDb", function(x){.taxonomyId(x)})





#########################################################################
## New method (experimental) to just see if we can make it easier for
## people who have RANGES and then want to see the associated
## annotations.
## Eventually, we may want to let the user choose which annotation
## range accessor should be called to see if their ranges overlap
## (with an 'annotFUN' argument).
## BUT RIGHT NOW: this will just do the simplest possible thing:

## issues:
## 1) exons, transcripts returns redundant results... (I really want
## exonsBy(x, by='gene') and then collapse the result.  Unfortunately, this means that the metadata is not in the mcols.
## 2) use transcriptsBy() for 'genes' (more accurate) - (currently they are both offered)

## 3) utrs and introns (similar issues to #1 above)
## 4) for utrs and introns, what you get back is grouped by transcript.  So I need to be able to re-group the results by gene OR (failing that, just call findOverlaps on *that* and the do:
## 5) And then I also (separately) need to be able to get the columns for these genes by calling select and then compressing that result to a DataFrame that can be put into mcols.

.selectByRanges <- function(x, ranges, columns=c('ENTREZID','SYMBOL'),
                            overlaps=c('gene','tx','exon', 'cds',
                                       'intron','5utr','3utr'),
                            ignore.strand=FALSE){
    ## Make sure everyone is OK with overlaps as argument name
    ##    overlaps <- match.arg(overlaps, several.ok = TRUE)
    overlaps <- match.arg(overlaps)
    subj <- switch(overlaps,
                  gene=genes(x,columns=columns),
                  exon=exonsBy(x,columns=columns,by='gene',outerMcols=TRUE),
                  cds=cdsBy(x,columns=columns,by='gene',outerMcols=TRUE),
                  tx=transcriptsBy(x,columns=columns,by='gene',outerMcols=TRUE),
                   ## the next three all return GRL grouped by transcripts...
                   '5utr'=fiveUTRsByTranscript(x),
                   '3utr'=threeUTRsByTranscript(x),
                   intron=intronsByTranscript(x)
                   )
    ## Then get the mcols
    if(overlaps %in% c('gene','tx','exon', 'cds')){
    ## Next do the overlaps                    
        hits <- findOverlaps(query=ranges, subject=subj,
                             ignore.strand=ignore.strand)
        results <- ranges[queryHits(hits)]
        mcols(results) <- mcols(subj[subjectHits(hits)])
    }else{ ## then it's mapped to transcripts so:
        ## Here we have to get our metadata set up and compressed FIRST
        keys <- names(subj)  ## keys are NOT unique)
        ## Get basic metadata mapped to TXID
        meta <- select(x, keys=keys, columns=columns,
                       keytype='TXID')
        ## Then compress based on the TXID keytype
        fa <- factor(meta[['TXID']], levels=unique(as.character(keys)))
        metaC <- .compressMetadata(fa, meta, avoidID='TXID')
        ## Then attach this compressed data onto the subject
        mcols(subj) <- metaC
        ## THEN we can do our overlaps
        hits <- findOverlaps(query=ranges, subject=subj,
                             ignore.strand=ignore.strand)
        results <- ranges[queryHits(hits)]
        mcols(results) <- mcols(subj[subjectHits(hits)])
        ## because we mapped by TXIDs, we have to remove redundant results
        dfRes <- as(results,'data.frame')
        uniqueIdx = !duplicated(dfRes)
        results <- results[uniqueIdx]
    }
    results
}


setMethod("selectByRanges", "OrganismDb",
          function(x,ranges,columns,overlaps,ignore.strand){
              if(missing(overlaps)){ overlaps <- 'tx' }
              if(missing(columns)){ columns <- c('ENTREZID','SYMBOL') }
              if(missing(ignore.strand)){ ignore.strand <- FALSE }
       .selectByRanges(x,ranges,columns,overlaps,ignore.strand)})



## ## Some Testing
## library(Homo.sapiens);
## ranges <-  GRanges(seqnames=Rle(c('chr11'), c(2)),IRanges(start=c(107899550, 108025550), end=c(108291889, 108050000)), strand='+', seqinfo=seqinfo(Homo.sapiens))

## selectByRanges(Homo.sapiens, ranges, 'SYMBOL')
## selectByRanges(Homo.sapiens, ranges, 'SYMBOL', 'exon')
## selectByRanges(Homo.sapiens, ranges, 'ENTREZID')
## ## What if they ask for something more compex?
## selectByRanges(Homo.sapiens, ranges, 'ALIAS')
## ## What if they ask for a couple things?
## selectByRanges(Homo.sapiens, ranges, c('ENTREZID','ALIAS'))

## The following should all give the same basic answer (because ranges
## don't change)

## selectByRanges(Homo.sapiens, ranges, c('SYMBOL','PATH'), 'tx')
## selectByRanges(Homo.sapiens, ranges, c('SYMBOL','PATH'), 'exon')
## selectByRanges(Homo.sapiens, ranges, c('SYMBOL','PATH'), 'cds')

## selectByRanges(Homo.sapiens, ranges, c('SYMBOL','PATH'), '5utr')
## selectByRanges(Homo.sapiens, ranges, c('SYMBOL','PATH'), '3utr')
## selectByRanges(Homo.sapiens, ranges, c('SYMBOL','PATH'), 'intron')



## Current troubles:

## debug(OrganismDbi:::.selectByRanges)
## selectByRanges(Homo.sapiens, ranges, c('SYMBOL','PATH'), '5utr')


## 1) it doesn't currently seem to respect strand (strand is being lost somewhere?
## 2) for UTRs/introns there is a problem where the results are 'redundant'.  I need a simple way to simplify them (considering both the contents of mcols AND the ranges) before returning them to the user.

## 3) The documentation and unit tests need a big upgrade...

## 4) I need an early version of Vinces complementary function still (selectRangesBy)




## FOR LATER sticky I want to implement support for multiple values of
## 'overlaps'

## So I want to be able to do something kind of like this (to implement the geometry idea of multiple 'overlaps' arguments).  BUT: it doesn't respect the contents of mcols...
## foo = selectByRanges(Homo.sapiens, ranges, c('SYMBOL','PATH'), '5utr');
## bar =  selectByRanges(Homo.sapiens, ranges, c('SYMBOL','PATH'), 'tx');
## unique(c(foo, bar))

## ALTERNATIVELY: I *could* implement this by just using the transcript centered strategy that I used above (for UTRs/introns), but applying it to 'everything', THEN merging all the tx centered metadata into a tx ID'd list and then overlapping as the last step.


## Older notes about this from the sticky:
## 'overlaps' CAN be a vector (which will result in multiple ranges getting summed together). - this suggestion is going to have to be it's own entirely separate sticky b/c the standard mechanisms for combining the results do *not* currently have any mechanism for respecting the geometry.  EITHER THAT, or I am going to have to change the way that the whole function works (again), by handling everything more the way that I currently handle things for 5UTRs/3UTRs/introns.  (IOW get all the ranges, always grouped by transcript, then combine to form one transcript oriented list and *then* annotate them, and then (at the end): overlap.












#######################################################################
## I need a way to get the inner mcols back out to the outer mcols (and quickly)


## Herve has a helper for extracting 'inner' mcols out of a GRanges
## list very quickly.
makeOuterMcolFromInnerMcol <- function(x, colname)
{
    if (!is(x, "List"))
        stop("'x' must be a List object")
  ##  tmpOri <- unique(relist(mcols(unlist(x, use.names=FALSE))[[colname]], x))
    tmp <- unique(relist(as.character(mcols(unlist(x, use.names=FALSE))[[colname]]),x))
    if (any(elementLengths(tmp) != 1L))
       stop(colname, " inner metadata column cannot be made an outer metadata column")
    unlist(tmp)
}

## Let's try with the exon_rank metadata col on the object returned by
## exonsBy() (should return an error):

##   ex_by_tx <- exonsBy(txdb)
##   mcols(ex_by_tx)[["exon_rank"]] <- makeOuterMcolFromInnerMcol(ex_by_tx, "exon_rank")

## but this should work with an inner metadata column that is really
## an attribute of the top-level list elements.


## unfortunately, this function seems to have some bugs. (Which I think I have mostly fixed)

## BUT: The function needs to extract all the viable mcols, and to format them as a DataFrame so that they can be put into the outer mcols for the results object. 

## AND even if I get this function working perfectly, I need to check do one of TWO things for each column.  If the column is from the inner column level (exon_rank) then I can't return that data in the result (since it won't map back out to the gene level).  These inner columns are not lost.  You might think that there is no sensible way to map them out to the result in the function above but they can just go into a integerList object...  So things like 'EXONRANK' will have to be processed differently...  In the case there the data is actually repeated from the outer column level then I need to use a variant of this function to map it back out.  So: two different things need to happen based on whether the data is repeated or not...


## the helpers exonsBy and transcriptsBy etc. need to be 'fixed' so that (for viable mcols) they have their outer mcols populated.  This will help since for the annotation recover, the outer mcols are the only ones I will want to use anyways.  This is still true for things like EXONRANK since for EXONRANK I will have an integerList (for example). The bottom line is: everthing in that outer mcols needs to be annotated at the 'GROUP level' regardless of what is in the inner mcols...  Once I have these base methods acting better it should be easy for my methor to do the right thing...

## there appears to be another separate bug that happens with exonsBy(x, by='tx') where the extra columns are not fully populated.  This needs to be fixed but doesn't happen with by='gene'...

## stash some private variables to hold the information about which columns are viable and which ones are not (for this).  This will help me to dispatch on columns that nee to be treated separately

## I may need a different accessor to list 'outer' columns, or I may need to add an argument to columns (geneLevel=TRUE)
