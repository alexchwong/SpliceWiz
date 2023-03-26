#' @describeIn covDataObject-class Constructs a covDataObject object.
#' For internal use only
#' @export
covDataObject <- function(
        args = list(),
        annotations = list(),
        colData = data.frame(),
        covData = list(),
        juncData = list(),
        normData = list()        
) {
    newCovObj <- new("covDataObject",
        args = args,
        annotations = annotations,
        colData = colData,
        covData = covData, juncData = juncData,
        normData = normData
    )
    newCovObj
}