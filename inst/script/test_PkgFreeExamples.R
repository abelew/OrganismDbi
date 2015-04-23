## Old example
library(OrganismDbi)

## debug(OrganismDbi:::.mungeGraphData)
## Bug fixed. - But ANOTHER one still lurks!

gd <- list(join1 = c(GO.db="GOID", org.Hs.eg.db="GO"),
           join2 = c(org.Hs.eg.db="ENTREZID",
             TxDb.Hsapiens.UCSC.hg19.knownGene="GENEID"))

makeOrganismPackage(pkgname = "Homsaps",
  graphData = gd,
  organism = "Homo sapiens",
  version = "1.0.0",
  maintainer = "Pkg Maintainer<maintainer@somewhere.org>",
  author = "Some Body",
  destDir = ".",
  license = "Artistic-2.0")


## ## This remaining bug causes this to fail: :(
## ## first an example that *will* work (traditional)
## gd <- list(join1 = c(GO.db="GOID", org.Hs.eg.db="GO"))
## makeOrganismPackage(pkgname = "sapiens",
##                     graphData = gd,
##                     organism = "Homo sapiens",
##                     version = "1.0.0",
##                     maintainer = "Pkg Maintainer<maintainer@somewhere.org>",
##                     author = "Some Body",
##                     destDir = ".",
##                     license = "Artistic-2.0")



## This is not something that I want to support? IOW, I think I want to be able to make OrganismDb objects, but NOT to make packages (unless there are real package depenedencies)
## and then one that requires we not depend on packages...
## debug(OrganismDbi:::.initPkg)
## debug(OrganismDbi:::.loadOrganismDbiPkg)
library(OrganismDbi)
library(AnnotationHub)
ah <- AnnotationHub()
puffer <- ah[['AH13961']]
##gd <- list(join1 = c(GO.db="GOID", puffer="GO"))

gd <- list(join1 = c(GO.db="GOID", puffer="GO"),
           join2 = c(puffer="ENTREZID",
             TxDb.Hsapiens.UCSC.hg19.knownGene="GENEID")) ## irrelevant TxDb
makeOrganismPackage(pkgname = "puff",
                    graphData = gd,
                    organism = "Tetraodon nigroviridis",
                    version = "1.0.0",
                    maintainer = "Pkg Maintainer<maintainer@somewhere.org>",
                    author = "Some Body",
                    destDir = ".",
                    license = "Artistic-2.0")

## then install gives me an error because its trying to install it without access to puffer....
install.packages("./puff/", repos=NULL)



####################################################################
## Lets make a simpler test that just uses the constructor:
## We need to be able to make and load this without an entire package...

## Here is what shows up in zzz.R (from .onLoad hook):
## So its just saving a data.frame (nothing special - but I *DO* need it...
load(system.file("data","graphData.rda",package='puff'))
## OrganismDbi:::.loadOrganismDbiPkg(pkgname='puff',
##                                  graphData=graphData)

## So basically: export this constructor here so that it can be called separately.
## Note: dbType argument is not even used...
## debug(OrganismDbi:::OrganismDb)
debug(OrganismDbi:::.initPkg)
debug(OrganismDbi:::.addLocalPkgsToNamespace)
debug(OrganismDbi:::.loadOrganismDbiPkg)

## install.packages("./puff/", repos=NULL)
library("puff")
puff

## puff <- OrganismDbi:::OrganismDb(graphData=graphData)
columns(puff)







########################################################################
## Now lets make an example with an Annotationhub resource, but that does not store the resource locally... (instead looks it up each time from the hub cache).

## But before I can do that, I need to make another change...  Basically, I need to change the internals of the object so that instead of storing the 'names' of the objects (which could change or eventually clash etc.) I instead use dbfile() to get (and store) the locations of the database files that exist for each one.  For those cases where an object does not yet have a file location, THEN we will do what we are doing now (and save it to the package/rely on it being in the path).

## so the plan will be to rely more on something like this:
## Does it have a value for dbfile()? (this covers installed packages AND AnnotHub objects)

## I am not sure if I can support objects that are only on the path (and unsaved):
## Ideally if it is *only* on the path (IOW there is no value for dbfile)? If so, then save it to the package (only can offer this saving option in case where they call makeOrganismPackageFromXXX() family... So users will be able to make an object that can't load the next time (if saved)

## So that would mean that the heuristic would be:
## 1 - is it on the path already? (exists())
## 2 - if not call loadDb on the file.path in the object

## And internally, the object would store either a name OR a character vector (IOW file path)

## Recommendation in that case will be to them to use the constructor only with reliable resources (packages or AnnotationHub objects).


## Alternatively I could just REQUIRE that all objects all be based on objects that have been saveDb/loadDb cycled...  But I think I want to keep the current flexibility.




## WHAT I REALLY WANT:

## all resources in a OrganismDb object must point to saved DBs (file paths only).
## if the end user calls makeOrganismPackageFromXXX() and passes in an object that is not saved to disc, it will be saved to their package
## if the user tries to make an object with something that is not saved to disc (using the constructor) then they should get a warning (but the constructor should basically just call it with exists() etc. (as it basically does now). 


## ??? Also the arguments (even for the general case) need to be a little different since I will now need to pass in actual objects (and not just pairs of strings).  This means I need to modify this to also have the actual objects passed in? Or maybe that's not all that important (since internally they will be object whether the actual object or just a name is passed in)...

## Here is what these are currently:
## gd <- list(join1 = c(GO.db="GOID", puffer="GO"),
##            join2 = c(puffer="ENTREZID",
##              TxDb.Hsapiens.UCSC.hg19.knownGene="GENEID")) ## irrelevant TxDb
