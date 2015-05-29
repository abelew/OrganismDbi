setGeneric("getTxDbIfAvailable", function(x, ...) standardGeneric("getTxDbIfAvailable"))

setGeneric("selectByRanges", signature="x",
           function(x, ranges, columns, overlaps, ignore.strand)
           standardGeneric("selectByRanges"))

setGeneric("selectRangesById", signature="x",
           function(x, keys, columns, keytype, overlaps)
           standardGeneric("selectRangesById"))
