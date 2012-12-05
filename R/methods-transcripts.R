## This is where I will put methods to overload things like
## transcripts() and exons()...



## new argument: cols here can be any legit value for cols. (not just
## tx_id and tx_name etc.)

## vals will just pass through to the internal transcripts call.

## columns arg is just for b/c support and will just pass through to
## the internal transcripts call


## I need helpers to extract that Txdb, and to remove cols that are
## part of a ranges object.

## 1st helper I need is one to identify which node (if any) is a TranscriptDb...
## if there is one, return it, otherwise error out with message to explain.
.findTranscriptDb <- function(x){
    

}


## 
setMethod("transcripts", "OrganismDb",
    function(x, vals=NULL, columns=c("tx_id", "tx_name")){
                
        ## Here we call the other transcripts method (on the TxDb)

        ## Then we call select on the rest.

        ## Then we can assemble it all together.
        
      }
)

