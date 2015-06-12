setGeneric("getTxDbIfAvailable", function(x, ...) standardGeneric("getTxDbIfAvailable"))

## getter/setter set
setGeneric("TxDb",function(x, ...) standardGeneric("TxDb"))
setGeneric("TxDb<-",signature="x",function(x, value) standardGeneric("TxDb<-"))


setGeneric("selectByRanges", signature="x",
           function(x, ranges, columns, overlaps, ignore.strand)
           standardGeneric("selectByRanges"))

setGeneric("selectRangesById", signature="x",
           function(x, keys, columns, keytype, feature)
           standardGeneric("selectRangesById"))
