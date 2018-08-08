setClassUnion("charnull", members = c("character", "NULL"))

#' TNS: An S4 class for survival analysis using transcriptional 
#' networks inferred by the RTN package.
#'
#' @slot TNI a previously computed \linkS4class{TNI}-class object.
#' @slot survivalData a data frame containing the survival data for all samples. 
#' Samples must be in the rows and the survival variables in the columns. 
#' Time of last update and event in last update (0 for alive, 1 for deceased).
#' @slot para a list with the parameters used to compute the GSEA2 analysis.
#' @slot results a list with results from TNS methods.
#' @slot status a vector containing the processing status of the TNS object.
#'
#' @method tnsPreprocess 
#' \code{\link{tnsPreprocess}}
#' @method tnsGSEA2 \code{\link{tnsGSEA2}}
#' @method tnsKM \code{\link{tnsKM}}
#' @method tnsPlotKM \code{\link{tnsPlotKM}}
#' @method tnsCox \code{\link{tnsCox}}
#' @method tnsPlotCox \code{\link{tnsPlotCox}}
#' @method tnsGet \code{\link{tnsGet}}
#' @method tnsInteraction \code{\link{tnsInteraction}}
#' @method tnsKmInteraction \code{\link{tnsKmInteraction}}
#' @method tnsPlotKmInteraction \code{\link{tnsPlotKmInteraction}}
#' @method tnsCoxInteraction \code{\link{tnsCoxInteraction}}
#' @method tnsPlotCoxInteraction \code{\link{tnsPlotCoxInteraction}}
#' @method tnsPlotGSEA2 \code{\link{tnsPlotGSEA2}}
#' @aliases TNS
#' 
#' @section Constructor:
#' 
#' see \code{\link{tni2tnsPreprocess}} constructor.
#'
#'
#' @exportClass TNS
#'
## Class TNS (Transcriptional Network - Survival)
setClass("TNS", 
         representation(TNI = "TNI", 
                        survivalData = "data.frame", 
                        para = "list", 
                        results = "list", 
                        status = "character"), 
         prototype = list(TNI = NULL, 
                          survivalData = data.frame(), 
                          para = list(), 
                          results = list(), 
                          status = character()
         )
)



