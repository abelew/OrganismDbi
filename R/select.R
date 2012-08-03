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

## helper to convert text strings into objects
.makeReal <- function(x){
  eval(parse(text=x))
}

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
  gd <- as.matrix(keyFrame(x))
  ## now give all the keys as a vector, but named by their databases.
  res <- c(gd[,3],gd[,4])
  pkgs <-  c(gd[,1],gd[,2])
  names(res) <- pkgs 
  res
}

## ## this method gives us a vector where names are the types,
## ## and the values are the keys
## .getDbObjFKeys <- function(x){
##   res <- .getDbNameFKeys(x) 
##   objs <- lapply(names(res), .makeReal)
##   names(res) <- lapply(objs, class)
##   res
## }



## this will return the cols, with appropriate things appended for a single Db.
## .getExtraColsForDb <- function(x, db, cols){
##   fkeys <- .getDbObjFKeys(x)
##   ## Now add the fkeys to the cols that go with the db
##   cols <- c(cols, fkeys[names(fkeys) %in% db])
##   unique(cols)
## }


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

## This one will actually get the extra cols for ALL the RELEVANT Dbs.
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
  ## Add extra packages to pkgs
  pkgs <- unique(c(.lookupDbNameFromKeytype(x, keytype), pkgs, extraPkgs))
  ## finally, for each node of pkgs, we need to grab the appropriate fkeys...
  fkeys <- .getDbNameFKeys(x) ## get all the fkeys
  fkeys <- fkeys[names(fkeys) %in% pkgs] ## only keep for pkgs in the path!
  ## And then add those keys to our cols
  unique(c(cols, fkeys))
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


.getSelects <- function(dbnames, keys, cols, keytype){
  res <- list(length(dbs))
  for(i in seq_len(length(dbs))){
    ## in addition to looping over the dbs, the appropriate cols must be
    ## selected for EACH
    dbtype <- as.character(class(dbs[[i]]))
    colsLocal <- cols[cols %in% cols(dbs[[i]])]
    if(i==1){
      ## mtype accumulates a history of tables that we have merged so far.
      mtype <- dbtype
      res[[i]] <- select(dbs[[i]], keys, colsLocal, keytype)
    }else{ ## more than one
      mtype <- c(mtype,dbtype)
      keytype <- mkeys[[paste(mtype,collapse="_")]][2] ## always the 2nd val
      ## An UGLY exception for GO.db:  (TODO: Is there a more elegant way?)
      if(dbtype=="GODb"){
        keytype="GOID"
      }
      prevKeyType <- mkeys[[paste(mtype,collapse="_")]][1]
      keys <- unique(res[[1]][[prevKeyType]])
      res[[i]] <- select(dbs[[i]], keys, colsLocal, keytype)
    }
  }
  names(res) <- sapply(dbs, class)
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
.mergeSelectResults <- function(sels, mkeys){
  for(i in seq_len(length(sels))){
    if(i==1){
      ## mtype starts with table i, and just accumulates a history of tables
      ## that we have merged so far.
      mtype <- names(sels)[i]
      res <- sels[[1]]
    }else{## There is more than one, so we must merge...
        mtype <- c(mtype,names(sels)[i])
        xkey <- mkeys[[paste(mtype,collapse="_")]][1]
        ykey <- mkeys[[paste(mtype,collapse="_")]][2]
        res <- merge(res, sels[[i]],by.x=xkey, by.y=ykey, all=TRUE)
                                        #, suffixes = c("",""))
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
.mkeys <- function(x, tbl1, tbl2, key=c("tbl1","tbl2")){
  key <- match.arg(key)
  kf <- OrganismDbi:::keyFrame(x)
  ## process for a double match of tbl1 and tbl2 (in any order)
  ## note: (we should ALWAYS have one when this function is called)
  parseCol <- function(piece, str){
    grepl(str, piece)
  }
  
  res <- apply(kf[,1:2], MARGIN=1, FUN=parseCol, tbl1)
  res2 <- apply(kf[,1:2], MARGIN=1, FUN=parseCol, tbl2)
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
  ## TODO: investigate setting default keys and cols arguments to NULL.  This
  ## is what the other select statements support, BUT it is more complicated
  ## here because we don't know which things the keys and cols go with at this
  ## point.  (So the check will have to happen later, after we know this
  ## information, and for cols it may be split up...)  
  ## if(is.null(keys)) keys <- keys(x) ## if no keys provided: use them all
  ## if(is.null(cols)) cols <- cols(x) ## if no cols provided: use them all
  ## check that the keytype matches the keys
  ## ktKeys = keys(x, keytype=keytype)
  ## if(!(any(ktKeys %in% keys))){
  ##   stop("keys must be of the same keytype as the actual keytype")
  ## }

  mkeys <- .mkeys()
  
  ## Preserve original cols (we will be adding some to get our results along
  ## the way 
  oriCols <- cols  
  
  ## 1st add any missing foreign key cols (based on what our graph looks like,
  ## and also based on what was asked for)
  cols <- .addAppropriateCols()
  ## Now we only need to get the nodes that *have* those columns.
  dbs <- .lookupDbsFromCols(x, cols, keytype)
  
  ## next we get the data from each.
  sels <- .getSelects(dbs, keys, cols, keytype, mkeys)
  ## Then we need to merge them together using the foreign keys
  res <- .mergeSelectResults(sels, mkeys)
  
  
  ## Then we need to filter out all columns that we didn't ask for.  
  ## Actually that is not quite right, what we want to do is make a blacklist
  ## of columns that were added (in fkeys) and that were NOT requested
  ## (oriCols and keytype).

#  extraKeys <- unique(unlist(mkeys))
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
    ## TODO: reqCols is kind of being faked here.  What we really want is to
    ## better clean up colnames(res) and to pass in a SHORTER list here.  Then
    ## the shorter list will act as a filter to "clean up" the columns by
    ## getting rid of unwanted columns. at least one kind of kruft that ends
    ## up in here is cols that are duplicated like: Ontology.x, Ontology.y
    ## etc.
  }
  unique(res)
}




setMethod("select", "OrganismDb",
    function(x, keys, cols, keytype) {
          .select(x, keys, cols, keytype)
        }
)


## debug(Homo.sapiens:::.select)


## planned usage:
## These all work now:
## library(Homo.sapiens);  cols <- cols(Homo.sapiens)[c(7,10,11,12)]; keys <- head(keys(org.Hs.eg.db, "ENTREZID")); keytype <- "ENTREZID"; res <- select(Homo.sapiens, keys, cols, keytype); head(res); dim(res)

## This 1st example should give me: [1] 51518    11
##(for the dim)



## cols2 <- cols(Homo.sapiens)[c(7,10,11,37)]; res2 <- select(Homo.sapiens, keys, cols2, keytype); head(res2)



## cols5 <- cols(Homo.sapiens)[c(7,8)]; res5 <- select(Homo.sapiens, keys, cols5, keytype); head(res5)



## cols3 <- cols(Homo.sapiens)[c(10,11,37)]; res3 <- select(Homo.sapiens, keys, cols3, keytype); head(res3)

## cols4 <- cols(Homo.sapiens)[c(10,11,12)]; res4 <- select(Homo.sapiens, keys, cols4, keytype); head(res4)



##########################
## Weird finTab issues??? - Resolved with AnnotationDbi patch

## cols6 <- cols(Homo.sapiens)[c(37,38)]; res6 <- select(Homo.sapiens, keys, cols6, keytype); head(res6)

## library(Homo.sapiens); keys <- head(keys(org.Hs.eg.db, "ENTREZID")); keytype <- "ENTREZID"; cols8 <- cols(Homo.sapiens)[c(37)]; res8 <- select(Homo.sapiens, keys, cols8, keytype); head(res8)


##########################
## BLERGH!  These issues are caused by this bug here (patch is avail)
## select(org.Hs.eg.db, keys= head(keys(org.Hs.eg.db, "ENTREZID")), cols = "ENTREZID", keytype="ENTREZID")



## This is the same situation as I patched above in AnnotationDbi, but now
## occuring in the meta-select I have created in the local .select

## cols7 <- cols(Homo.sapiens)[c(10)]; res7 <- select(Homo.sapiens, keys, cols7, keytype); head(res7)






###################################################################
## TODO: more testing starting from other types of IDs!

## library(Homo.sapiens);cols3 <- cols(Homo.sapiens)[c(10,11,37)];cols2 <- cols(Homo.sapiens)[c(7,10,11,37)];cols8 <- cols(Homo.sapiens)[c(37)];

## debug(Homo.sapiens:::.lookupDbsFromCols);
## debug(Homo.sapiens:::.mergeSelectResults);
## debug(Homo.sapiens:::.select);
## debug(Homo.sapiens:::.getSelects);
## debug(Homo.sapiens:::.addAppropriateCols);

## strange bug: GOID should be keytype, but only GO actually works...
## foo = head(keys(GO.db))
## this fails:
## res9 <- select(Homo.sapiens, keys=foo, cols=cols8, keytype="GOID"); head(res9)

## this also bombs:
## res9 <- select(Homo.sapiens, keys=foo, cols=cols3, keytype="GOID"); head(res9)


## But this works? - oops not since I added the GODb exception above to .getSelects...  Damn, that's confusing...
## res10 <- select(Homo.sapiens, keys=foo, cols=cols8, keytype="GO"); head(res10)
## And so does this? Wait, no.  It doesn't really work either...
## res11 <- select(Homo.sapiens, keys=foo, cols=cols3, keytype="GO"); head(res11)




## What about starting with a key from TxDb???
## 

## works
## keys = head(keys(Homo.sapiens, "CDSID")); res12 <- select(Homo.sapiens, keys=keys, cols=cols8, keytype="CDSID"); head(res12)


## But this one doesn't finish in a timely fashion?  (seems to now be an
## efficiency thing)
## keys = head(keys(Homo.sapiens, "TXID")); res13 <- select(Homo.sapiens, keys=keys, cols=cols3, keytype="TXID"); head(res13)





## keys = head(keys(Homo.sapiens, "TXID")); res14 <- select(Homo.sapiens, keys=keys, cols=cols2, keytype="TXID"); head(res14)

## This most recent one is an interesting bug, we are dropping all the results
## simply because we couldn't connect all the way through.  




## Here is what seems to be bombing on these last two:
## system.time(foo <- select(org.Hs.eg.db,keys=keys(org.Hs.eg.db,keytype="ENTREZID"), cols=c("ENTREZID","ACCNUM","GO"), keytype="ENTREZID"))

## so after cols improvements I am trying: 
## foo <- select(org.Hs.eg.db, keys=keys(org.Hs.eg.db,keytype="ENTREZID"), cols=c("ENTREZID","ACCNUM"), keytype="ENTREZID")
## And this is still WAY too slow.  In fact, I really need to speed that up in AnnotationDbi...  But in the meantime lets implement solution 2 here.


## proposed efficiencies:
## I should be able to help things by NOT getting all of the keys?
## length(keys)==42106 is a LOT of keys. (and that is indeed ALL of them)
## I should also be able to save some time by being more selective
## about the cols (don't get GO unless I really must).

### for both? of these things, I should be able to do something smart
### based only on knowing whether or not I am processing an OrgDb
### (ie. whether or not this node as more than one edge) and also based
### on whether or not I have another select to process after this one
### (if not, then we don't need to get another table).

## MAJOR TODO: selectively drop fkeys based on how many selects we will
## need...  There is just no reason to retrieve GO IDs if we don't need them
## (for example)


## SOLUTION 1:
## Actually I think I can do better criteria.  I know I need a GO ID "added"
## only IF 1) there is a GODb in the chain somewhere 2) the chain is longer
## than 1 AND 3) we are on the OrgDb node.  The flip side to this would be the
## case where I only want to add the ENTREZID link, which would happen in the
## complementary situation.  A more general solution is to see that 1) we are
## on an OrgDb node, 2) this is not the 1st link (if it is we don't need to
## add anything), 3) what is the length of dbs (later this would be outdegree
## to current node) and then 4) if it is 2, then we need to choose which one
## using the mkeys (helper function)
## DONE



## Deciding about row-key filtering is much harder.  It means that I
## have to start passing in the keys from .select, and then I have to
## consider repercussions intelligently.  But I think I can do it
## simply because I think that I basically always want an inner join
## here between the selects.  That is, I start with some key in some
## table, and I want to return all the things from the other tables
## that can be matched (inner join).  This is the best I can do.

## SOLUTION 2:
## So I start with the keys in the 1st table (it will be the 1st table
## because I have sorted them).  And then I just use mkeys to get the
## key types (for the previous table) and use that info to get the
## column from the sels[[i-1]][[prevKeyType]].  And THOSE are my new
## keys...
## DONE


## ALSO: I have a problem with things that were NOT in the original query
## ending up in the result.  See res8, and notice how I get Evidence, and
## Ontology included in the results even though I did succeed at filtering out
## "GO".  I have a finite list of exceptions like this that can end up in the
## results without having been in the initial "cols" argument.  I probably
## need a special helper method to scrub these exceptions out whenever "GO" is
## no longer present in the colnames etc.  This isn't my favorite solution,
## but it would probably work very well most of the time since all of my
## exceptions are coming from Org Packages and will be known about in advance.

## Solution 3:
## write a post-process filter-helper to remove these when appropriate.









































## ANOTHER external bug.  Crikey!

## debug(Homo.sapiens:::.getSelects)
## There may be a problem with
## txdb= TxDb.Hsapiens.UCSC.hg19.knownGene; foo = select(txdb, keys=keys(txdb,keytype="GENEID"), cols=c("CDSSTART", "GENEID"),keytype="GENEID"); head(foo); head(foo)
## txdb= TxDb.Hsapiens.UCSC.hg19.knownGene; foo = select(txdb, keys=keys(txdb,keytype="GENEID"), cols=c("TXSTART", "GENEID"),keytype="GENEID"); head(foo); head(foo)


## BUT:
## library(TxDb.Hsapiens.UCSC.hg19.knownGene); x <- txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene;
## cols = c("GENEID","TXID","TXSTART"); keys = head(keys(x, "GENEID")); foo = select(x, keys, cols = cols, keytype="GENEID");head(foo)

## SO THIS is the problem!  (using a lot of keys in a TXDB select() call is apparently unforgivably slow...
##  foo = select(x, keys(x, "GENEID"), cols = cols, keytype="GENEID");head(foo)






## TODO #1: bug is impeding progress because select() for GO.db is not giving us
## what we want.

## test code
## head(toTable(GOTERM))
## head(select(GO.db, keys = head(keys(GO.db)), cols="TERM"))
## BPPARENTS ??? WTH did that come from???
## Adding in more stuff still leaves us with the bizarre rename event...
## colnames(select(GO.db, keys = head(keys(GO.db)), cols=c("TERM","SYNONYM")))
## debug(AnnotationDbi:::.select)

## trouble seems to be coming from here:
## debug(AnnotationDbi:::.renameColumnsWithRepectForExtras)

## digging a bit deeper,
## debug(AnnotationDbi:::.getAllColAbbrs)
## debug(AnnotationDbi:::.makeColAbbrs)
## These methods seem to do what they are supposed to.
## AnnotationDbi:::.getAllColAbbrs(GO.db)
## But it seems I have a logical error.  I cannot use the redundant short
## column names to map back (accurately) to the table abbreviations when they
## are keys used to join across multiple tables...

## For the AnnotationDbi implementation, there is a general problem when I
## have multiple join columns that are used in different mappings...
## For the general mappings: accession can take us to REFSEQ *OR* to ACCNUM,
## and ipi_id can take us to PFAM *OR* PROSITE. These I could almost live
## with, but the really killer one is GO.db, where almost EVERYTHING has go_id
## mapped, and the one you get (the 1st one) is BPPARENTS, which produces that
## WTH momemt from earlier...

## To fix it, I really just need to make some exceptions.  If the column ID is
## go_id, is should return "TERM", and if it is accession or ipi_id, we
## probably need to make a smarter choice about how we are getting our cols
## back out.
## May have to look in the value for cols to get the right decision???
## So for PFAM, and PROSITE, or REFSEQ and ACCNUM this sort of column sub as
## the last step in .renameColumnsWithRepectForExtras() is a clean solution.
## Still testing for GO though...


## colnames(select(GO.db, keys = head(keys(GO.db)), cols=c("TERM")))
## debug(AnnotationDbi:::.renameColumnsWithRepectForExtras)

## colnames(select(org.Hs.eg.db, keys = head(keys(org.Hs.eg.db)), cols=c("REFSEQ", "ACCNUM")))

## colnames(select(org.Hs.eg.db, keys = head(keys(org.Hs.eg.db)), cols=c("REFSEQ", "ACCNUM", "GO")))

## x = GO.db
## AnnotationDbi:::.renameColumnsWithRepectForExtras(x, toTable(GOTERM))
## This too is now FIXED.


## FOR GO, I ALSO have another bug too (asking for multiple things results in
## ambiguous merge() calls). Example:
## colnames(select(GO.db, keys = head(keys(GO.db)), cols=c("TERM","CCPARENTS")))
## This bug may relate to the fact that calling toTable() on some GO Bimaps
## returns more than one column with the same exact name...  EXAMPLE
## colnames(toTable(GOCCPARENTS))
## colnames(toTable(GOPBPARENTS))
## colnames(toTable(GOTERM))
## To fix these kinds of situations, I need to do a pre-filter to remove the
## unwanted go_id columns (they show up after you say toTable - not sure why).
## Once such redundant columns are removed (sometime before we try to merge
## them) we should be able to do a swap like above, but to just replace go_id
## with "GOID" - FIXED (this part at least).





################################################################################
## HUGE TODO: Write a generic way to track which foreign keys are needed to
## stitch indiv. selects together and use that for the merge above.
## One idea is to use graphs objects


## Another huge problem is GO.  Right now we have 'GOID' for GO.db and 'GO' as
## a keytype for org and chip pkgs.  That is ugly, and we can't change that
## without namespacing the cols and keytypes for both groups.  Because then if
## you want to use GO IDs and get from an org package to a transcriptDB
## package it might be confused with using a GO ID and starting from the GO.db
## package.  (a different path)


## THREE constraints for graphs that we could use to sort out the joins on:
## 1) Must be acyclic
## 2) Only can have one edge between any two nodes
##    (where the edge represents a relationship between foreign keys).
## 3) There must be a unique (namespaces?) way of naming IDs that connect nodes
##    (ie. GOID and GO)



## as an input I think I will ask for just a data.frame (from the user).
## Then internally I can call ftM2graphNEL().  That way users won't have to
## call new() to make a graphNEL
## the three collumns can be DB1, DB2, and key to use for merging.
## Paul recommends that I use graphNELs (I can use internally)
## graphNEL

## Paul also recommends that use 
## RBGL for a shortest path algorithm.  dijkstra.sp() ?



