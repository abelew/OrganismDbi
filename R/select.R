## This will just hold code for the initial implementation of select and friends

## Some helpers for retrieval of data.

## 1st some getters
setGeneric("keyFrame", function(x) standardGeneric("keyFrame"))
setMethod("keyFrame", "OrganismDb",
    function(x){x@keys}
)
setGeneric("dbGraph", function(x) standardGeneric("dbGraph"))
setMethod("dbGraph", "OrganismDb",
    function(x){x@graph}
)


## Then some helpers to process some of these results a bit
.getDbObjNames <- function(x){
  gd <- as.matrix(keyFrame(x))
  unique(c(gd[,1],gd[,2]))
}

.getDbObjs <- function(x){
  dbs <- .getDbObjNames(x)
  objs <- lapply(dbs, .makeReal)
  names(objs) <- dbs
  objs
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


















## For select then, I will need to also look up repeats for both cols and
## keytypes.  For keys, I want them uniqued together, but for cols I will need
## to fully qualify them when they are returned (TODO: changes to cols).  That
## means that users have to indicate which one they want with a fully
## qualified name when they call select.  Alternatively, they could just get
## "both" back...  Right NOW: I am still lucky (even for cols) in terms of
## name clashes.

## FOR NOW: we don't have any name clashes, so we will proceed as if there
## never will be any (when we get one, we can start to name-mangle etc.)


## So select will work by just knowing how to merge together the results from
## the separate select() calls.  I will need know know who can be joined with
## who for this (and how).  That info should be stored in the
## OrganismDb when it is created.  

## vectorized keytype->DB matching BUT ALSO: we have to sort this so
## that we will do select() first for the DB that we actually have a key for
.makecolMapping <- function(x){
    objs <- .getDbObjs(x)
    res <- lapply(objs, cols)
    unlist2(res)
}

.lookupDbNamesFromCols <- function(x, col){
  cols <- cols(x)
  if(!all(col %in% cols)){
    stop("col must be a value returned by cols(x).")
  }
    res <- .makecolMapping(x)
    ## no duplicates so I can just return the name
    names(res)[res %in% col]
}


.getDistances <- function(x, keytype){
  ## use the shortest path algorithm to get them into the order desired.
  g <- dbGraph(x)
  startNode <- .lookupDbNameFromKeytype(x, keytype)
  sp <- dijkstra.sp(g, start=startNode)
  sp$distances
}

## now uses graphs to order things better
.resortDbs <- function(x, pkgs, keytype){
  ## use the shortest path algorithm to get them into the order desired.
  res <- .getDistances(x, keytype)
  ## So now I just have to sort the names in order of distance.
  dbs <- names(res)[order(res)]
  ## Then convert into actual objects
  objs <- lapply(dbs, .makeReal)
  ## Then label the objs with the names
  names(objs) <- dbs
  ## Then filter out any objs that were not listed in pkgs (clip leaf nodes)
  objs <- objs[names(objs) %in% pkgs]
  ## Then return 
  objs
}

## Retrieves the list of DBs "in order of shortest path"
.lookupDbsFromCols <- function(x, cols, keytype){
  ## 1st we want cols ordered so that the one that matches our keytype is FIRST
  ## Also, we must always have keytype be part of cols here
  cols <- unique(c(keytype, cols))
  pkgs <- .lookupDbNamesFromCols(x, cols)
  pkgs <- unique(pkgs)
  ## Now use the graph to decide the path
  .resortDbs(x, pkgs=pkgs, keytype=keytype)
}
## .lookupDbsFromCols(x, cols, keytype)

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
## graph-aware .addAppropriateCols() function.

## I have upgraded .addAppropriateCols() so that it pays attention to the graph.
## IOW, the cols asked for should not just be all of them but only the ones
## that are needed based on the initial values of cols that was passed into
## select.

## SO: Instead of just returning ALL the fkey columns I need to implement the
## following rule:
## IF the cols requested are separated on the graph by too big a distance:
## THEN I need to get the keys from the intermediate nodes.


## helper to get only the linked keys.
.dropUnlinkedKeys <- function(x, fkeys, pkgs, dsts){
  ## run through the sorted final pkgs list, and label each fkey as supported
  ## or not.
  newFkeys <- character()
  for(i in seq_len(length(pkgs)-1)){  ## Or one test per edge
    ## if distance from keytype is unequal: collect keys from this edge.
    if(dsts[[pkgs[i]]]!=dsts[[pkgs[i+1]]]){ 
      key1 <- .mkeys(x, pkgs[i],pkgs[i+1], key = "tbl1")
      key2 <- .mkeys(x, pkgs[i],pkgs[i+1], key = "tbl2")
      newFkeys <- c(newFkeys,key1,key2)
    }
  }
  unique(newFkeys)
}

## This one will actually get the extra cols (foreign keys) for ALL the
## RELEVANT Dbs and then add that information to the keytype and cols.
.addAppropriateCols <- function(x, cols, keytype){
  ## get the graph
  g <- dbGraph(x)
  ## We want to run dijkstras (factor from .resortDbs)
  dst <- .getDistances(x, keytype)
  ## then we need to lookup the DBs we need (based on the cols alone). (
  pkgs <- .lookupDbNamesFromCols(x, cols)
  ## helper function for matching and subsetting
  ## BE CAREFUL! the order matters for match and is reverse of %in%
  ## (iow use the long vector second)
  matchSub <- function(x, m1, m2){
    idx <- match(unlist(m1), unlist(m2))
    idx <- idx[!is.na(idx)]
    x[idx]
  }
  ## nodeWalker recurses to get back to the start node
  nodeWalker <- function(g, dst, pkg, curDist,extraPkgs){ 
    ## get the edges 
    e <- edges(g) 
    enames <- names(e)
    dstNames <- names(dst)

    pkgConts <- matchSub(e, pkg, enames) 
    
    pkgDists <-  matchSub(dst, pkgConts, dstNames) 
    pkg <- names(pkgDists)[pkgDists < curDist]    
    curDist <- matchSub(dst, pkg, dstNames)    
    if(curDist !=0){
      extraPkgs <- c(extraPkgs, pkg)
      nodeWalker(g, dst, pkg, curDist, extraPkgs)
    }else{
      return(extraPkgs)
    }
  }
  ## vector for holding nodes that are "between" start and leaf nodes.
  extraPkgs <- character()
  ## master loop for all leaf pkgs (leaf nodes)
  for(i in seq_len(length(pkgs))){
    pkg <- pkgs[i] 
    curDist <- dst[names(dst) %in% pkg] 
    if(curDist > 1){## then we need to work out how to 
      extraPkgs <- nodeWalker(g, dst, pkg, curDist, extraPkgs) 
    }
  }
  ## Now combine together all the different packages 
  pkgs <- unique(c(.lookupDbNameFromKeytype(x, keytype), pkgs, extraPkgs))
  ## And then order the pkgs so that they are sorted from keytype to leaves
  pkgSubDsts <- dst[match(pkgs, names(dst))]
  sortedPkgs <- pkgs[order(pkgSubDsts)]
  
  ## finally, for each node of pkgs, we need to grab the appropriate fkeys...
  fkeys <- .getDbNameFKeys(x) ## 1st get ALL the fkeys
  ## Then I need to drop keys that point to pkgs that are NOT in the path.
  fkeys <- .dropUnlinkedKeys(x, fkeys, sortedPkgs, pkgSubDsts)
  ## And then add those keys to our cols
  unique(c(keytype, cols, fkeys))
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


.getSelects <- function(x, dbs, keys, cols, keytype){
  res <- list(length(dbs))
  for(i in seq_len(length(dbs))){
    ## in addition to looping over the dbs, the appropriate cols must be
    ## selected for EACH
    dbtype <- names(dbs)[[i]]
    colsLocal <- cols[cols %in% cols(dbs[[i]])]
    if(i==1){
      ## mtype accumulates a history of tables that we have merged so far.
      mtype <- dbtype
      res[[i]] <- select(dbs[[i]], keys, colsLocal, keytype)
    }else{ ## more than one
      mtype <- c(mtype,dbtype)
      ml <- length(mtype)
      keytype <- .mkeys(x, mtype[ml-1], mtype[ml], key="tbl2")
      ## An UGLY exception for GO.db:  (TODO: Is there a more elegant way?)
      if(dbtype=="GODb"){
        keytype="GOID"
      }
      prevKeyType <- .mkeys(x, mtype[ml-1], mtype[ml], key="tbl1")
      keys <- unique(res[[1]][[prevKeyType]])
      res[[i]] <- select(dbs[[i]], keys, colsLocal, keytype)
    }
  }
  names(res) <- names(dbs)
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

.mkeys <- function(x, tbl1, tbl2, key=c("tbl1","tbl2")){
  key <- match.arg(key)
  kf <- keyFrame(x)
  ## process for a double match of tbl1 and tbl2 (in any order)
  ## note: (we should ALWAYS have one when this function is called)
  
  res <- apply(kf[,1:2], MARGIN=1, FUN=.parseCol, tbl1)
  res2 <- apply(kf[,1:2], MARGIN=1, FUN=.parseCol, tbl2)
  res <- res | res2 
  resRowIdx <- res[,1] & res[,2]
  matchRow <- kf[resRowIdx,]
  if(dim(matchRow)[1]<1){stop("No relationship found for ",tbl1," and ",tbl2)}

  ## now the tricky part is that in returning the keys I have to get the
  ## correct keys back to the user...  And this is based on whether tbl1 was
  ## one thing or another.
  if(key=="tbl1"){
    if(grepl(tbl1,matchRow$xDbs)){
      return(as.character(matchRow$xKeys))
    }else{ ## then its reversed of the order in the row...
      return(as.character(matchRow$yKeys))
    }
  }else if(key=="tbl2"){
    if(grepl(tbl2,matchRow$yDbs)){
      return(as.character(matchRow$yKeys))
    }else{ ## and the reverse case
      return(as.character(matchRow$xKeys))
    }
  }
}

.select <- function(x, keys, cols, keytype){
  ## if asked for what they have, just return that.
  if(all(cols %in% keytype)  && length(cols)==1){
    res <- data.frame(keys=keys)
    colnames(res) <- cols
    return(res) }
  
  ## Preserve original cols (we will be adding some to get our results along
  ## the way 
  oriCols <- cols  
  
  ## 1st add any missing foreign key cols (based on what our graph looks like,
  ## and also based on what was asked for)
  cols <- .addAppropriateCols(x, cols, keytype)
  ## Now we only need to get the nodes that *have* those columns.
  dbs <- .lookupDbsFromCols(x, cols, keytype)
  
  ## next we get the data from each.
  sels <- .getSelects(x, dbs, keys, cols, keytype)
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



