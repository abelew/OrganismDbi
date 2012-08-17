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




## For keys I need to be a little smarter...  There will eventually be cases
## where the same keys are used in two locations.  What should we do then?  I
## think we HAVE to have a strategy for this.  One option would be to
## determine if the keys really "ARE" the same as each other (some tests could
## verify this).  If they are determined to be the same (one list is a subset
## of the other OR if there is substantial overlap - this second case would
## merit a warning), then we could just unique them together, but in this case
## we still have a much more complicated situation to deal with
## !!!
## I think that in the above case you just remove keys like that...


## Another strategy is just not to allow keytypes with the same name when you
## make the OrganismDb object...


## Strategy for keys: I need a lookup function that can 1) generate the keys
## for each slot and then lookup which slot I should be tapping based on a
## keytype.  2) This needs to be general purpose (will be needed again in
## select) and 3) it also may need to be able to return multiple hits in the
## event that there are eventually multiple IDs named the same way (depends on
## whether or not we allow repeat ID names).  I think we WILL want to allow
## this, which means I will have to return repeats, and then in the subsequent
## part of keys will have to call keys multiple times and unique them
## together.


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
.makecolMapping <- function(x){
    objs <- OrganismDbi:::.getDbObjs(x)
    res <- lapply(objs, cols)
    unlist2(res)
}

.lookupDbNamesFromCols <- function(x, cls){
  cols <- cols(x)
  if(!all(cls %in% cols)){
    stop("All values for cls must be a value returned by cols(x).")
  }
    res <- .makecolMapping(x)
    ## no duplicates so I can just return the name
    res <- res[res %in% cls]
    ## BUT, the names need to come back in same order as cls
    names(res)[match(cls,res)]
}


.getDistances <- function(x, keytype){
  ## use the shortest path algorithm to get them into the order desired.
  g <- dbGraph(x)
  startNode <- .lookupDbNameFromKeytype(x, keytype)
  sp <- dijkstra.sp(g, start=startNode)
  sp$distances
}


## newer version just uses names of DBs from fkeys
## Retrieves the list of DBs "in order of shortest path"
.lookupDbsFromFkeys <- function(x, fkeys, keytype){
  if(length(fkeys)==0){
    pkgs <- .lookupDbNameFromKeytype(x, keytype)
  }else{
    pkgs <- .lookupDbNamesFromCols(x, fkeys)
  }
##  pkgs <- unique(pkgs)
  dbs <- lapply(pkgs, .makeReal)
  names(dbs) <- pkgs
  dbs
}
## .lookupDbsFromFkeys(x, fkeys)


## library(Homo.sapiens); x= Homo.sapiens; cols=c("GOID","ENTREZID","TXNAME"); keytype="ENTREZID"; 


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

## helper function for matching and subsetting
.matchSub <- function(x, m1, m2){
  idx <- match(unlist(m1), unlist(m2))
  idx <- idx[!is.na(idx)]
  x[idx]
}

## nodeWalker recurses to get back to the start node
.nodeWalker <- function(g, dst, pkg, curDist, extraKeys=list()){ 
  prevPkg <- pkg
  ## get the edges 
  e <- edges(g) 
  enames <- names(e)
  dstNames <- names(dst)
  
  pkgConts <- .matchSub(e, pkg, enames) 
  
  pkgDists <-  .matchSub(dst, pkgConts, dstNames)
  pkg <- names(pkgDists)[pkgDists < curDist]    
  curDist <- .matchSub(dst, pkg, dstNames)    
  if(curDist !=0){
    extraKeys[[length(extraKeys)+1]] <- rev(.mkeys(x, prevPkg, pkg, key="both"))
    .nodeWalker(g, dst, pkg, curDist, extraKeys=extraKeys)
  }else{
    extraKeys[[length(extraKeys)+1]] <- rev(.mkeys(x, prevPkg, pkg, key="both"))
    return(extraKeys)
  }
}

## This one will actually get the extra cols (foreign keys) for ALL the
## RELEVANT Dbs and then add that information to the keytype and cols.
.getForeignEdgeKeys <- function(x, cols, keytype){ 
  ## get the graph 
  g <- dbGraph(x)
  ## We want to run dijkstras 
  dst <- .getDistances(x, keytype)
  ## then we need to lookup the DBs we need (based on the cols alone). 
  pkgs <- .lookupDbNamesFromCols(x, cols)
  fKeys <- list()
  ## master loop for all leaf pkgs (leaf nodes) 
  for(i in seq_len(length(pkgs))){
    pkg <- pkgs[i] 
    curDist <- dst[names(dst) %in% pkg] 
    if(curDist > 0){## then we are not there yet... 
      fKeys <- c(fKeys,
                 rev(.nodeWalker(g, dst, pkg, curDist))) 
    }
  }
  ## then put the extrKeys together with the other things we need 
  ## the listed form (not returned) is nice for looking at the edges formed 
  unique(unlist(fKeys)) 
} 




## also a helper to filter out cols that are duplicated after select()
## TODO: migrate this version of the method down to AnnotationDbi so that we
## never get column duplicates
.dropDuplicatedgetSelectsCols <- function(tab){
  cols <- colnames(tab)
  cols <- cols[!duplicated(cols)]  
  tab <- tab[,cols]
  tab
}



.splitBy1stNode <- function(dbs, fkeys){
  res <- names(dbs)
  names(res) <- fkeys
  root = res[1]
  tf <- cumsum(as.numeric(res %in% root))
  split(res, as.factor(tf))
}


## I NEED to always start with the keytype...  Otherwise I don't get to take
## advantage of the efficiency from pulling only the keys I actually need into
## R...

    ## compute this path through the graph
##     require(graph)
##     gr <- dbGraph(x)
##     sgr <- subGraph(names(dbs), gr)
##     require(RBGL)
##     bfp <- bfs(sgr, names(dbs)[1]) ## should yield my path



.getSelect <- function(x, dbs, cols, keytype, res){
  
  for(i in seq_len(length(dbs))){
    dbtype <- names(dbs)[[i]]
    ## in addition to looping over the dbs, the appropriate cols must be
    ## selected for EACH
    colsLocal <- cols[cols %in% cols(dbs[[i]])]
 
    ## start node is always the keytype
    if(i==1){
      ## prev records the db that we last used
      prev <- dbtype
      if(!(names(dbs)[[i]] %in% names(res))){
         sel <- select(dbs[[i]], keys, colsLocal, keytype=keytype)
         res <- c(res, sel)
#         names(res)[[i]] <- dbtype
       }
    }else{ ## more than one
      prev <- names(dbs)[[i-1]]
      kt <- .mkeys(x, prev, dbtype, key="tbl2")
#       kt <- fkeys[[i-1]][dbtype]
      prevKeyType <- .mkeys(x, prev, dbtype, key="tbl1")
#       prevKeyType <- fkeys[[i-1]][prev]
      ## THIS LINE RIGHT HERE needs to use actual res 
      ## but how to do that for the 1st pass??
      keys <- unique(res[[prev]][[prevKeyType]])
      if(!(names(dbs)[[i]] %in% names(res))){
        sel <- select(dbs[[i]], keys, colsLocal, keytype=kt)
        res <- c(res, sel)
#        names(res)[[i]] <- dbtype
      }
    }
   } 
  res
}

## dbs and fkeys should now be in SAME ORDER (dbs were derived from fkeys)
## Also, fkeys was in order of the chromosome walk.
.getSelects <- function(x, keys, cols, keytype){
  ## get the dbs
  fkeys <- OrganismDbi:::.getForeignEdgeKeys(x, cols, keytype)
  dbs <- .lookupDbsFromFkeys(x, fkeys, keytype)
  ## Then split that up according to the occurance of the 1st node.
  walks <- .splitBy1stNode(dbs, fkeys)
  ## from this point forward cols needs to be comprehensive
  cols <- unique(c(keytype, cols, fkeys))
  ## results will be in a list structure
  res <- list()

  ## Then for each list element call the code below.
  ## 1) track the nodes to avoid calling the same node twice.  Do this by
  ## naming the list elements as I fill them and checking against that for
  ## each node.
  
  ## 2) always assume for each node that we have to deduce the keys from what
  ## was selected before using the "walk" that we are currently on.
 
  ## 3) Add extra loop: For each "walk" call the following (for each of the
  ## dbs in that walk), but making sure to not call a dbs that we have already
  ## put into the result: "res"
  for(i in seq_len(length(walks))){
    walkDbNames <- walks[[i]]
    walkDbs <- dbs[match(walkDbNames,names(dbs))]
    ## Remove from walkDbs based on what is in res already
##    walkDbs <- walkDbs[!(names(walkDbs) %in% names(res))] ## NOT HERE.
    res <- c(res, .getSelect(x, walkDbs, cols, keytype, res))
    walkDbShort <- walkDbs[!(names(walkDbs) %in% names(res))]
    names(res) <- names(walkDbShort)
  }
##   names(res) <- unique(names(dbs)) ## looks scary but should be correct 
  res 
} 

##.dropDuplicatedMergeCols helper just takes advantage of the fact that my
##dupicates columns will always have the form a.x, a.y etc. and will drop the
##.y columns and then keep the .x ones.
.dropDuplicatedMergeCols <- function(tab){
  cols <- colnames(tab)  
  cols <- cols[!grepl(".y", cols)]
  tab <- tab[,cols]
  colnames(tab) <- gsub(".x","",colnames(tab))
  tab
} 


## This merges things based on the key relationships from mkeys
.mergeSelectResults <- function(x, sels){
  for(i in seq_len(length(sels))){
    if(i==1){
      ## mtype starts with table i, and just accumulates a history of tables
      ## that we have merged so far.
      mtype <- names(sels)[i]
      res <- sels[[1]]
    }else{## There is more than one, so we must merge...
        mtype <- c(mtype,names(sels)[i])
        ml <- length(mtype)
        xkey <- .mkeys(x, mtype[ml-1], mtype[ml], key="tbl1")
        ykey <- .mkeys(x, mtype[ml-1], mtype[ml], key="tbl2")
        res <- merge(res, sels[[i]],by.x=xkey, by.y=ykey, all=TRUE)
        res <- .dropDuplicatedMergeCols(res)
      } 
  }
  res
} 



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
  
  ## 1st add any missing foreign key cols (based on what our graph looks like,
  ## and also based on what was asked for)
#  fkeys <- .getForeignEdgeKeys(x, cols, keytype)
  
  ## Then I need to add back the keytype and cols into cols???
#  cols <- unique(c(keytype, cols, fkeys))
    
  ## next we get the data from each.
  sels <- .getSelects(x, keys, cols, keytype)
  ## Then we need to merge them together using the foreign keys
  res <- .mergeSelectResults(x, sels)
  
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





## new plan:

## library(RBGL)
## cls <- c("GOID" ,  "SYMBOL", "TXNAME")
## kt <- "ENTREZID"
## keys <- head(keys(x, "ENTREZID"))


## helper for getting all cols by all nodes
colsByNodes <- function(x){
  gr <- OrganismDbi:::dbGraph(x)
  allCols <- lapply(nodes(gr), function(elt) cols(OrganismDbi:::.makeReal(elt)))
  names(allCols) <- nodes(gr)
  allCols
}
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

## Derive the from to data.frame from the graph
## need all combinations of edges
getFTDF <- function(subgr){
  xx = strsplit(edgeNames(subgr, recipEdges="distinct"), "~")
  data.frame(
             from=sapply(xx, "[[", 1),
             to = sapply(xx, "[[", 2),
             stringsAsFactors=FALSE)
}
## ftDf <- getFTDF(subgr)

## get list of nodes to visit
## visitNodes = bfs(subgr, root)

## set up an empty list with names that match what we want to fill...
## selected = setNames(
##   vector("list", length(visitNodes)),
##   visitNodes)


## new version of .getSelects()
## ## select
## node = visitNodes[[1]]
## selected[[node]] =
##   select(OrganismDbi:::.makeReal(node),
##          keys, needCols[[node]], kt)

## for (node in visitNodes[-1]) {
##   fromNode = ftDf[ftDf$to == node, "from"]
##   fromKey = OrganismDbi:::.mkeys(x, fromNode, node, "tbl1")
##   fromKeys = unique(selected[[fromNode]][[fromKey]])
##   fromKeys = fromKeys[!is.na(fromKeys)]
##   toKey = OrganismDbi:::.mkeys(x, fromNode, node, "tbl2")
##   selected[[node]] =
##     select(OrganismDbi:::.makeReal(node),
##            fromKeys, needCols[[node]], toKey)
## }


## new version of .mergeSelects
## ## merge
## final = selected[[1]]
## for (node in visitNodes[-1]) {
##   fromNode = ftDf[ftDf$to == node, "from"]
##   fromKey = OrganismDbi:::.mkeys(x, fromNode, node, "tbl1")
##   toKey = OrganismDbi:::.mkeys(x, fromNode, node, "tbl2")
##   final = merge(final, selected[[node]], by.x=fromKey, by.y=toKey, all=TRUE)
## }

