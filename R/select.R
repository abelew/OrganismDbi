## This will just hold code for the initial implementation of select and friends




## helper to convert text strings (Db pkgs names) into real objects
.makeReal <- function(x){
  eval(parse(text=x))
}



## Standard methods:
.keytypes <- function(x){
  dbs <- .getDbObjs(x)
  unique(unlist(lapply(dbs, keytypes)))
}

setMethod("keytypes", "OrganismDb",
    function(x) .keytypes(x)
)

## Usage:
## keytypes(Homo.sapiens)



.cols <- function(x){
  dbs <- .getDbObjs(x)
  unique(unlist(lapply(dbs, cols)))
}

setMethod("cols", "OrganismDb",
    function(x) .cols(x)
)

## Usage:
## cols(Homo.sapiens)




## Strategy for keys: I need a lookup function that can 1) generate the keys
## for each slot and then lookup which slot I should be tapping based on a
## keytype.  2) This needs to be general purpose (will be needed again in
## select) and 3) it also may need to be able to return multiple hits in the
## event that there are eventually multiple IDs named the same way (depends on
## whether or not we allow repeat ID names).  I think we WILL want to allow
## this, which means I will have to do some kind of name-spacing scheme.


.makekeytypeMapping <- function(x){
    objs <- .getDbObjs(x)
    res <- lapply(objs, keytypes)
    unlist2(res)
}

.lookupDbNameFromKeytype <- function(x, keytype){
  kts <- keytypes(x)
  if(!(keytype %in% kts)){
    stop("keytype must be a value returned by keytypes(x).")
  }  
    res <- .makekeytypeMapping(x)
    ## no duplicates so I can just return the name
    names(res)[res %in% keytype]  
}

.lookupDbFromKeytype <- function(x, keytype){
    db <- .lookupDbNameFromKeytype(x, keytype)
    eval(parse(text=db))
}

.keys <- function(x, keytype){
    kts <- keytypes(x)
    if(!(keytype %in% kts)){
      stop("supplied keytype is not allowed.  Please call the keytypes method to see which keytypes are allowed for this organism.")
    }
    ## We need to retrieve the keys from the relevant slot...
    ## So 1st we have to determine WHO has the keys we are after
    if(length(keytype) !=1){
      stop("The keys method can only accept one keytype at a time.")
    }
    db <- .lookupDbFromKeytype(x, keytype)
    ## And then we can just call keys...
    keys(db, keytype)
}

setMethod("keys", "OrganismDb",
    function(x, keytype){.keys(x, keytype)}
)


## Usage: (so far this works, but there is something dumb happening with the method
## head(keys(Homo.sapiens, keytype="PMID"))
## the use case for GOID will present a special challenge...
## head(keys(Homo.sapiens, keytype="GOID"))









## Select will work by just knowing how to merge together the results from
## the separate select() calls.  I will need know know who can be joined with
## who for this (and how).  That info should be stored in the
## OrganismDb when it is created.  

## vectorized keytype->DB matching BUT ALSO: we have to sort this so
## that we will do select() first for the DB that we actually have a key for
## .makecolMapping <- function(x){
##     objs <- OrganismDbi:::.getDbObjs(x)
##     res <- lapply(objs, cols)
##     unlist2(res)
## }

## .lookupDbNamesFromCols <- function(x, cls){
##   cols <- cols(x)
##   if(!all(cls %in% cols)){
##     stop("All values for cls must be a value returned by cols(x).")
##   }
##     res <- .makecolMapping(x)
##     ## no duplicates so I can just return the name
##     res <- res[res %in% cls]
##     ## BUT, the names need to come back in same order as cls
##     names(res)[match(cls,res)]
## }


## .getDistances <- function(x, keytype){
##   ## use the shortest path algorithm to get them into the order desired.
##   g <- dbGraph(x)
##   startNode <- .lookupDbNameFromKeytype(x, keytype)
##   sp <- dijkstra.sp(g, start=startNode)
##   sp$distances
## }


## newer version just uses names of DBs from fkeys
## Retrieves the list of DBs "in order of shortest path"
## .lookupDbsFromFkeys <- function(x, fkeys, keytype){
##   if(length(fkeys)==0){
##     pkgs <- .lookupDbNameFromKeytype(x, keytype)
##   }else{
##     pkgs <- .lookupDbNamesFromCols(x, fkeys)
##   }
## ##  pkgs <- unique(pkgs)
##   dbs <- lapply(pkgs, .makeReal)
##   names(dbs) <- pkgs
##   dbs
## }
## ## .lookupDbsFromFkeys(x, fkeys)


## ## library(Homo.sapiens); x= Homo.sapiens; cols=c("GOID","ENTREZID","TXNAME"); keytype="ENTREZID"; 


## This method just gets me the pkg names as names and vals are fkeys
.getDbNameFKeys <- function(x){
  gd <- keyFrame(x)
  ## now give all the keys as a vector, but named by their databases.
  .extractPkgsAndCols(gd)
}



## If the user chooses cols from GO and Txdb (for example), I
## still need to join those via a dbs that they need no data from (org pkg in
## this example).  So in that case, certain cols need to be incluced by a
## graph-aware .getForeignEdgeKeys() function.

## I have upgraded .getForeignEdgeKeys() so that it pays attention to the graph.
## IOW, the cols asked for should not just be all of them but only the ones
## that are needed based on the initial values of cols that was passed into
## select.

## ## helper function for matching and subsetting
## .matchSub <- function(x, m1, m2){
##   idx <- match(unlist(m1), unlist(m2))
##   idx <- idx[!is.na(idx)]
##   x[idx]
## }

## ## nodeWalker recurses to get back to the start node
## .nodeWalker <- function(g, dst, pkg, curDist, extraKeys=list()){ 
##   prevPkg <- pkg
##   ## get the edges 
##   e <- edges(g) 
##   enames <- names(e)
##   dstNames <- names(dst)
  
##   pkgConts <- .matchSub(e, pkg, enames) 
  
##   pkgDists <-  .matchSub(dst, pkgConts, dstNames)
##   pkg <- names(pkgDists)[pkgDists < curDist]    
##   curDist <- .matchSub(dst, pkg, dstNames)    
##   if(curDist !=0){
##     extraKeys[[length(extraKeys)+1]] <- rev(.mkeys(x, prevPkg, pkg, key="both"))
##     .nodeWalker(g, dst, pkg, curDist, extraKeys=extraKeys)
##   }else{
##     extraKeys[[length(extraKeys)+1]] <- rev(.mkeys(x, prevPkg, pkg, key="both"))
##     return(extraKeys)
##   }
## }

## ## This one will actually get the extra cols (foreign keys) for ALL the
## ## RELEVANT Dbs and then add that information to the keytype and cols.
## .getForeignEdgeKeys <- function(x, cols, keytype){ 
##   ## get the graph 
##   g <- dbGraph(x)
##   ## We want to run dijkstras 
##   dst <- .getDistances(x, keytype)
##   ## then we need to lookup the DBs we need (based on the cols alone). 
##   pkgs <- .lookupDbNamesFromCols(x, cols)
##   fKeys <- list()
##   ## master loop for all leaf pkgs (leaf nodes) 
##   for(i in seq_len(length(pkgs))){
##     pkg <- pkgs[i] 
##     curDist <- dst[names(dst) %in% pkg] 
##     if(curDist > 0){## then we are not there yet... 
##       fKeys <- c(fKeys,
##                  rev(.nodeWalker(g, dst, pkg, curDist))) 
##     }
##   }
##   ## then put the extrKeys together with the other things we need 
##   ## the listed form (not returned) is nice for looking at the edges formed 
##   unique(unlist(fKeys)) 
## } 




## ## also a helper to filter out cols that are duplicated after select()
## ## TODO: migrate this version of the method down to AnnotationDbi so that we
## ## never get column duplicates
## .dropDuplicatedgetSelectsCols <- function(tab){
##   cols <- colnames(tab)
##   cols <- cols[!duplicated(cols)]  
##   tab <- tab[,cols]
##   tab
## }



## .splitBy1stNode <- function(dbs, fkeys){
##   res <- names(dbs)
##   names(res) <- fkeys
##   root = res[1]
##   tf <- cumsum(as.numeric(res %in% root))
##   split(res, as.factor(tf))
## }


## I NEED to always start with the keytype...  Otherwise I don't get to take
## advantage of the efficiency from pulling only the keys I actually need into
## R...

    ## compute this path through the graph
##     require(graph)
##     gr <- dbGraph(x)
##     sgr <- subGraph(names(dbs), gr)
##     require(RBGL)
##     bfp <- bfs(sgr, names(dbs)[1]) ## should yield my path





## ##.dropDuplicatedMergeCols helper just takes advantage of the fact that my
## ##dupicates columns will always have the form a.x, a.y etc. and will drop the
## ##.y columns and then keep the .x ones.
## .dropDuplicatedMergeCols <- function(tab){
##   cols <- colnames(tab)  
##   cols <- cols[!grepl(".y", cols)]
##   tab <- tab[,cols]
##   colnames(tab) <- gsub(".x","",colnames(tab))
##   tab
## } 








## .mkeys will return appropriate value "on the fly" based on the
## contents of keyFrame().  It will take at least three arguments: the two
## tables plus an indicator for which of the two keys 1st or 2nd table key is
## needed.

## tbl1,tbl2 wil be actual package names like 'org.Hs.eg.db' or 'GO.db'
## tbl1 = "TxDb.Hsapiens.UCSC.hg19.knownGene"
## tbl2 = "org.Hs.eg.db"
## key = "tbl1"
.parseCol <- function(piece, str){
  grepl(str, piece)
} 

.mkeys <- function(x, tbl1, tbl2, key=c("tbl1","tbl2", "both")){
  if(length(tbl1) >1 || length(tbl2)>1) stop(".mkeys can only process one pair of tables at at time")
  key <- match.arg(key)
  kf <- keyFrame(x)
  ## process for a double match of tbl1 and tbl2 (in any order)
  ## note: (we should ALWAYS have one when this function is called)
  
  res <- apply(kf[,1:2], MARGIN=2, FUN=.parseCol, tbl1)
  res2 <- apply(kf[,1:2], MARGIN=2, FUN=.parseCol, tbl2)
  fin <- res | res2 
  resRowIdx <- fin[,1] & fin[,2]
  matchRow <- kf[resRowIdx,]
  if(dim(matchRow)[1]<1){stop("No relationship found for ",tbl1," and ",tbl2)}

  ## now the tricky part is that in returning the keys I have to get the
  ## correct keys back to the user...  And this is based on whether tbl1 was
  ## one thing or another.
  if(length(matchRow$xDbs) >1)stop("mkeys has failed to limit choices to 1.")
  if(key=="tbl1"){
    if(grepl(tbl1,matchRow$xDbs)){
      ans <- as.character(matchRow$xKeys)
    }else{ ## then its reversed of the order in the row...
      ans <- as.character(matchRow$yKeys)
    }
  }else if(key=="tbl2"){
    if(grepl(tbl2,matchRow$yDbs)){ 
      ans <- as.character(matchRow$yKeys)
    }else{ ## and the reverse case
      ans <- as.character(matchRow$xKeys)
    }
  }else if(key=="both"){
    ans <- c(as.character(matchRow$xKeys),as.character(matchRow$yKeys))
    names(ans) <- c(as.character(matchRow$xDbs),as.character(matchRow$yDbs))
    ## When we say "both" we still want keys returned in same order as
    ## original packages.  IOW, if tbl1 goes with key 1, then we should list
    ## key 1 1st in the result...
    ans <- ans[match(c(tbl1,tbl2),names(ans))]
  }
  ans
}





###############################################
## NEW plan helpers:
###############################################

## helper for getting all cols by all nodes
colsByNodes <- function(x){
  gr <- OrganismDbi:::dbGraph(x)
  allCols <- lapply(nodes(gr), function(elt) cols(OrganismDbi:::.makeReal(elt)))
  names(allCols) <- nodes(gr)
  allCols
}
## library(Homo.sapiens)
## library(RBGL)
## library(graph)a
## x = Homo.sapiens
## allCols <- colsByNodes(x)

## helper to get the subgraph
getRelevantSubgraph <- function(x, cols, keys, keytype){
  gr <- OrganismDbi:::dbGraph(x)
  allCols <- colsByNodes(x)
  inSubgraph = sapply(allCols,
    function(cols, keys) any(keys %in% cols),  union(keytype, cols))
  subgr = subGraph(names(inSubgraph)[inSubgraph], gr)
  subgr
}
## kt <- "ENTREZID"
## cls = c("GOID" ,  "SYMBOL", "TXNAME")
## keys <- head(keys(x, "ENTREZID"))
## subgr <- getRelevantSubgraph(x, cols=cls, keys, keytype=kt)
  

## now we will also need the root
## root = OrganismDbi:::.lookupDbNameFromKeytype(x, kt)


## I think this is meant to be an lapply
getForeignKeys <- function(x, subgr){
  fKeys = lapply(strsplit(edgeNames(subgr), "~"),
    function(tables, x, key)
    OrganismDbi:::.mkeys(x, tables[[1]], tables[[2]], "both"), x)
  unlist(fKeys, use.names=FALSE)
}
## fKeys <- getForeignKeys(x, subgr)


## now combine all the keys together
## selectCols = unique(c(kt, fKeys, cls))

## sort the needed cols by their nodes
getColsByNodes <- function(subgr, selectCols, allCols){
  lapply(allCols[nodes(subgr)],
    function(col, selectCols) col[col %in% selectCols], selectCols)
}
## needCols <- getColsByNodes(subgr, selectCols, allCols)


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
.getSelects <- function(x, keytype,keys,needCols, visitNodes){
  ## set up an empty list with names that match what we want to fill...
  selected = setNames(
    vector("list", length(visitNodes)),
    names(visitNodes))
  ## in 1st case we only need the name
  node1 = names(visitNodes)[[1]]
  selected[[node1]] =
    select(OrganismDbi:::.makeReal(node1),
           keys=keys,
           cols=needCols[[node1]],
           keytype=keytype)
  ## but here we need to use the name and the value of visitNodes
  otherNodes <- visitNodes[-1] 
  for (i in seq_len(length(otherNodes))) {
    nodeName <- names(otherNodes)[i]
    fromNode = otherNodes[i] 
    fromKey = OrganismDbi:::.mkeys(x, fromNode, nodeName, "tbl1")
    fromKeys = unique(selected[[fromNode]][[fromKey]])
    fromKeys = fromKeys[!is.na(fromKeys)]
    toKey = OrganismDbi:::.mkeys(x, fromNode, nodeName, "tbl2")
    selected[[nodeName]] =
      select(OrganismDbi:::.makeReal(nodeName),
             keys=fromKeys,
             cols=needCols[[nodeName]],
             keytype=toKey)
  }
  selected
}
## selected <- .getSelect(kt,keys,needCols, visitNodes)


## new version of .mergeSelectResults
## ## merge
.mergeSelectResults <- function(x, selected, visitNodes){
  final = selected[[1]]
  otherNodes <- visitNodes[-1] 
  for (i in seq_len(length(otherNodes))) {
    nodeName <- names(otherNodes)[i]
    fromNode = otherNodes[i] 
    fromKey = OrganismDbi:::.mkeys(x, fromNode, nodeName, "tbl1")
    toKey = OrganismDbi:::.mkeys(x, fromNode, nodeName, "tbl2")
    final = merge(final, selected[[nodeName]],
                  by.x=fromKey, by.y=toKey, all=TRUE)
  }
  final
}
## res <- .mergeSelectResults(selected, visitNodes)





.select <- function(x, keys, cols, keytype){
  ## if asked for what they have, just return that.
  if(all(cols %in% keytype)  && length(cols)==1){
    res <- data.frame(keys=keys)
    colnames(res) <- cols
    return(res)
  }
  
  ## Preserve original cols (we will be adding some to get our results along
  ## the way 
  oriCols <- cols  

  ## New methods make more use of graph objects.
  allCols <- colsByNodes(x)
  subgr <- getRelevantSubgraph(x, cols=cols, keys, keytype=keytype)
  root = OrganismDbi:::.lookupDbNameFromKeytype(x, keytype)
  fKeys <- getForeignKeys(x, subgr)
  selectCols = unique(c(keytype, fKeys, cols))
  needCols <- getColsByNodes(subgr, selectCols, allCols)
  visitNodes = .bfs(subgr, root)
  selected <- .getSelects(x, keytype,keys,needCols, visitNodes)
  res <- .mergeSelectResults(x, selected, visitNodes)
  
  ## Then we need to filter out all columns that we didn't ask for.  
  ## Actually that is not quite right, what we want to do is make a blacklist
  ## of columns that were added (in fkeys) and that were NOT requested
  ## (oriCols and keytype).

  extraKeys <- .getDbNameFKeys(x)
  blackList <- extraKeys[!(extraKeys %in% unique(c(oriCols, keytype)))]
  ## if they asked for one of the GO items, then GO is not blacklisted
  if(any(cols(GO.db) %in% oriCols)){
    blackList <- blackList[!(blackList %in% "GO")]
  }
  res <- res[,!(colnames(res) %in% blackList)]

  ## Then call code to clean up, reorder the rows (and add NA rows as needed).
  if(dim(res)[1]>0){
    res <- AnnotationDbi:::.resort(tab=res, keys=keys, jointype=keytype,
                                   reqCols=colnames(res))

  }
  unique(res)
}




setMethod("select", "OrganismDb",
    function(x, keys, cols, keytype) {
          .select(x, keys, cols, keytype)
        }
)





