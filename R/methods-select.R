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
.getSelects <- function(x, keytype, keys, needCols, visitNodes){
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
               keytype=keytype)
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
                   keytype=toKey)
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


.select <- function(x, keys, cols, keytype, ...){
    ## Argument checking:
    if(missing(keys)){stop("You must provide a keys argument")}
    if(missing(cols)){stop("You must provide columns argument")}
    if(missing(keytype)){stop("You must provide a keytype argument")}
    ## Some argument checking
    AnnotationDbi:::.testSelectArgs(x, keys=keys, cols=cols, keytype=keytype)
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
    selected <- .getSelects(x, keytype,keys,needCols, visitNodes)
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
