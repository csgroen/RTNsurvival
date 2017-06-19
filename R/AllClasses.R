setClassUnion("TNInull", members = c("TNI", "NULL"))

#' TNS: An S4 class for survival survival analysis using transcriptional 
#' networks inferred by the RTN package.
#'
#' @slot survivalData a data frame containing the survival data for all samples. 
#' Samples must be in the rows and the survival variables in the columns. 
#' Time of last update and event in last update (0 for alive, 1 for deceased).
#' @slot EScores a list created by tnsGSEA2, containing the enrichment scores 
#' for all samples
#' @slot keycovar a string vector containing the key covariables used to 
#' compute the Cox regression. They must be present in the survivalData table.
#' @slot tni a \linkS4class{TNI}-class object, previously computed. It is added 
#' to the TNS via tnsGSEA2
#' @slot status a vector containing the processing status of the TNS object.
#' @slot para a list with the parameters used to compute the GSEA2 analysis.
#'
#' @method tnsPreprocess 
#' \code{\link[RTNsurvival:tnsPreprocess]{tnsPreprocess}}
#' @method tnsGSEA2 \code{\link[RTNsurvival:tnsGSEA2]{tnsGSEA2}}
#' @method tnsKM \code{\link[RTNsurvival:tnsKM]{tnsKM}}
#' @method tnsCox \code{\link[RTNsurvival:tnsCox]{tnsCox}}
#' @method tnsGet \code{\link[RTNsurvival:tnsGet]{tnsGet}}
#' @aliases TNS
#' 
#' @section Constructor:
#' 
#' tnsPreprocess(tni, survivalData, keycovar, time = 1,
#'               event = 2, samples = NULL)
#' 
#' \itemize{
#' \item {tni} - {a 'TNI' object.}
#' \item {survivalData} - {a data.frame containing at least 2 columns}
#' \item {keycovar} - {A character vector or NULL}
#' \item {time} - {A numeric value or character vector}
#' \item {event} - {A numeric value or character vector}
#' \item {samples} - {A character vector}
#' }
#'
#' @exportClass TNS
#'
## Class TNS (Transcriptional Network - Survival)
setClass("TNS", representation(tni = "TNInull", survivalData = "data.frame", keycovar = "character", 
    para = "list", EScores = "list", status = "character"), prototype = list(tni = NULL, 
    survivalData = data.frame(), keycovar = character(), para = list(), EScores = list(), 
    status = character()))
