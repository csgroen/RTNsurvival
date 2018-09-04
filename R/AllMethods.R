
#' Preprocessing of TNS class objects
#'
#' Creates TNS class onbjects for regulons an survival data.
#'
#' @param tni A \linkS4class{TNI} class, already processed with the same samples
#' listed in the survival data.frame.
#' @param survivalData A named data.frame with samples in rows and survival data 
#' in the columns (this does not need to be provided if avaibale in the 'TNI' object).
#' @param regulatoryElements A character vector specifying which 
#' 'TNI' regulatory elements should be evaluated.
#' @param time A numeric or character value corresponding to the column of the 
#' data.frame where the time of last observation is given.
#' @param event A numeric or character value, corresponding to the columm of 
#' the data.frame where the 'event' information is given.
#' @param endpoint A numeric value. It represents the cut-off point for the 
#' 'time', if any.
#' @param pAdjustMethod A single character value specifying the p-value 
#' adjustment method to be used (see 'p.adjust' function for details).
#' @param keycovar A character vector of 'keycovars' listed in 'survivalData' columns.
#' @param samples An optional character vector listing samples to be analyzed.
#' @param excludeMid A logical value. If TRUE, inconclusive dES values is not
#' consired in the survival analysis.
#' @return A preprocessed \linkS4class{TNS} class
#' @examples
#' # load survival data
#' data(survival.data)
#' 
#' # load TNI-object
#' data(stni, package = "RTN")
#' 
#' # create a new TNS object
#' stns <- tni2tnsPreprocess(stni, survivalData = survival.data, 
#'         keycovar = c('Grade','Age'), time = 1, event = 2)
#'
#' @seealso \code{\link{tni.preprocess}} for similar 
#' preprocessing.
#' @import methods
#' @importFrom RTN upgradeTNI tni.get
#' @docType methods
#' @rdname tni2tnsPreprocess-methods
#' @aliases tni2tnsPreprocess
#' @export
#' 
setMethod("tni2tnsPreprocess", "TNI", 
          function(tni, survivalData = NULL, regulatoryElements = NULL, 
                   time = 1, event = 2, endpoint = NULL, pAdjustMethod = "BH", 
                   keycovar = NULL, samples = NULL, excludeMid = FALSE)
          {
            #-- tni checks
            tni <- upgradeTNI(tni)
            .tns.checks(tni, type="TNI")
            if(!is.null(regulatoryElements)){
              .tns.checks(regulatoryElements, type="regulatoryElements")
              regulatoryElements <- .checkRegel(tni, regulatoryElements)
              tni@regulatoryElements <- regulatoryElements
            }
            
            #-- is.null
            if (is.null(survivalData)){
              res <- try(tni.get(tni, "colAnnotation"), silent = TRUE)
              if (class(res) == "try-error"){
                stop("Must provide a 'survivalData' object.")
              } else if(nrow(res)){
                survivalData <- res
              }
              survivalData <- res
            }
            
            #-- par checks
            .tns.checks(survivalData, type = "survivalData")
            time = .tns.checks(time, survivalData, type = "time")
            event = .tns.checks(event, survivalData, type = "event")
            .tns.checks(endpoint, type = "endpoint")
            .tns.checks(keycovar, survivalData, "Keycovars")
            samples = .tns.checks(samples, survivalData, type = "samples")
            .tns.checks(excludeMid, type = "excludeMid")
            .tns.checks(pAdjustMethod, type = "pAdjustMethod")
            
            #-- other checks
            if (!all(samples %in% colnames(tni@gexp))) {
              stop("NOTE: all sample names listed in 'survivalData' rownames must be available in the 'tni' object!")
            }
            
            #-- reorganize survivalData
            idx <- c(time, event)
            te.data <- survivalData[, idx]
            survivalData <- survivalData[, -idx]
            survivalData <- cbind(te.data, survivalData)
            names(survivalData)[1:2] <- c("time", "event")
            survivalData <- survivalData[samples, ]
            
            #-- set endpoint
            if(is.null(endpoint)){
              endpoint <- max(survivalData$time, na.rm = TRUE)
            }
            
            #-- making TNS object
            para <- list(time=time, event=event, endpoint=endpoint, 
                         keycovar=keycovar, excludeMid=excludeMid,
                         pAdjustMethod=pAdjustMethod)
            object <- new("TNS", TNI = tni, survivalData = survivalData, para = para)
            
            #-- status update
            object <- tns.set(object, what = "status")
            return(object)
          })


#' Compute regulon activity using 2-tailed Gene Set Enrichment Analysis
#'
#' Works as a wrapper for \code{\link{tni.gsea2}}, performing a 
#' 2-tailed GSEA analysis on a \linkS4class{TNI} class object and integrating 
#' the results into the \linkS4class{TNS} class object.
#'
#' @param tns A \linkS4class{TNS} class, which has been preprocessed
#' @param ... Additional parameters passed to \code{\link{tni.gsea2}} function.
#' @return A \linkS4class{TNS} class, with added regulon activity scores.
#' @examples
#' # load survival data
#' data(survival.data)
#' 
#' # load TNI-object
#' data(stni, package = "RTN")
#'
#' stns <- tni2tnsPreprocess(stni, survivalData = survival.data, 
#' keycovar = c('Grade','Age'), time = 1, event = 2)
#' stns <- tnsGSEA2(stns)
#'
#'\dontrun{
#'
#'# parallel version with SNOW package!
#'library(snow)
#'options(cluster=makeCluster(3, "SOCK"))
#'stns <- tnsGSEA2(stns)
#'stopCluster(getOption("cluster"))
#'
#'}
#'
#' @seealso \code{\link{tni.gsea2}} for information on all 
#' parameters.
#' @importClassesFrom RTN TNI
#' @importFrom RTN tni.gsea2
#' @docType methods
#' @rdname tnsGSEA2-methods
#' @aliases tnsGSEA2
#' @export
#'
setMethod("tnsGSEA2", "TNS", function(tns, ...) {
  
  #-- checks
  if (tns@status["Preprocess"] != "[x]") 
    stop("NOTE: TNS object requires preprocessing!")
  
  #-- update para
  para <- tnsGet(tns, what = "para")
  para$regulonActivity <- "gsea2"
  tns <- tns.set(tns, para, "para")
  
  #-- run gsea2 and update TNS
  tni <- tnsGet(tns, what = "TNI")
  regulonActivity <- tni.gsea2(tni, ...=...)
  tns <- tns.set(tns, regulonActivity, "regulonActivity")
  return(tns)
})


#' Compute regulon activity by calling aREA (analytic Rank-based Enrichment Analysis) algorithm
#'
#' Uses \code{\link{tni.area3}} function to compute regulon activity
#' for \linkS4class{TNS} class objects.
#'
#' @param tns A \linkS4class{TNS} class, which has been preprocessed
#' @param ... Additional parameters passed to \code{\link{tni.area3}} function.
#' @return A \linkS4class{TNS} class, with added regulon activity scores.
#' @references Alvarez et al. Functional characterization of somatic mutations in cancer
#' using network-based inference of protein activity. Nature Genetics, 48(8):838-847, 2016.
#' @examples
#' # load survival data
#' data(survival.data)
#' 
#' # load TNI-object
#' data(stni, package = "RTN")
#'
#' stns <- tni2tnsPreprocess(stni, survivalData = survival.data, 
#' keycovar = c('Grade','Age'), time = 1, event = 2)
#' 
#' stns <- tnsAREA3(stns)
#'
#' @seealso \code{\link{tni.area3}} for additional details.
#' @importClassesFrom RTN TNI
#' @importFrom RTN tni.area3
#' @importFrom scales rescale
#' @docType methods
#' @rdname tnsAREA3-methods
#' @aliases tnsAREA3
#' @export
#'
setMethod("tnsAREA3", "TNS", function(tns, ...){
  
  #-- checks
  if (tns@status["Preprocess"] != "[x]") 
    stop("NOTE: TNS object requires preprocessing!")
  
  #-- update para
  para <- tnsGet(tns, what = "para")
  para$regulonActivity <- "aREA"
  tns <- tns.set(tns, para, "para")
  
  #-- run area3
  tni <- tnsGet(tns, what = "TNI")
  regulonActivity <- tni.area3(tni, ...=...)
  
  #-- rescale to fit graphics
  regulonActivity$dif <- apply(regulonActivity$dif, 2, rescale, to=c(-1.8, 1.8))
  
  #-- update TNA
  tns <- tns.set(tns, regulonActivity, "regulonActivity")
  return(tns)
})


#' Kaplan-Meier analysis for TNS class objects
#'
#' Creates survival curves and tests if there is a difference between 
#' curves using 'survfit' and 'survdiff' functions, respectivelly.
#'
#' @param tns A \linkS4class{TNS} object, which must have passed GSEA2 analysis.
#' @param regs An optional string vector listing regulons to be tested.
#' @param nSections A numeric value for sample stratification. The larger
#' the number, the more subdivisions will be created for the Kaplan-Meier 
#' analysis.
#' @param verbose A logical value specifying to display detailed messages 
#' (when verbose=TRUE) or not (when verbose=FALSE).
#' 
#' @return Results from 'survfit' and 'survdiff', including log-rank statistics.
#' @examples
#' # load survival data
#' data(survival.data)
#' 
#' # load TNI-object
#' data(stni, package = "RTN")
#'
#' stns <- tni2tnsPreprocess(stni, survivalData = survival.data, 
#'         keycovar = c('Grade','Age'), time = 1, event = 2)
#' stns <- tnsGSEA2(stns)
#' stns <- tnsKM(stns)
#' tnsGet(stns, "kmTable")
#'
#' @importFrom survival survdiff survfit coxph Surv
#' @docType methods
#' @rdname tnsKM-methods
#' @aliases tnsKM
#' @export
#' 
setMethod("tnsKM", "TNS", 
          function(tns, regs = NULL, nSections = 1, verbose = TRUE){
            #-- checks
            .tns.checks(tns, type = "Activity")
            .tns.checks(regs, type = "regs")
            .tns.checks(nSections, type = "nSections")
            .tns.checks(verbose, type = "verbose")
            
            #-- run stratification
            tns <- tnsStratification(tns, nSections = nSections)
            
            #-- get data and para
            regulonActivity <- tnsGet(tns, what = "regulonActivity")
            survData <- tnsGet(tns, what = "survivalData")
            para <- tnsGet(tns, what = "para")
            endpoint <- para$endpoint
            excludeMid <- para$excludeMid
            pAdjustMethod <- para$pAdjustMethod
            
            #--- update para
            para$nSections <- nSections
            tns <- tns.set(tns, para, "para")
            
            #-- set endpoint
            survData$event[survData$time > endpoint] <- 0
            survData$time[survData$time > endpoint] <- endpoint
            
            #-- making reglist
            reglist <- colnames(regulonActivity$regstatus)
            if(!is.null(regs)) {
              if (!all(regs %in% reglist)) {
                stop("all names in 'regs' should be listed in the slot 'results$regulonActivity' of the 'tns' object!")
              }
              reglist <- regs
            }
            
            idx <- apply(regulonActivity$regstatus, 2, function(es) {
              any(is.na(es))
            })
            validregs <- colnames(regulonActivity$regstatus)[!idx]
            reglist <- reglist[reglist %in% validregs]
            
            if (verbose) {
              cat("Computing survival curves...\n")
              pb <- txtProgressBar(min = 0, max = length(reglist), style = 3)
            }
            
            #--- logrank
            kmFit <- list()
            kmTable <- NULL
            for(reg in reglist){
              res <- .survstats(regulonActivity, survData=survData, reg=reg, excludeMid=excludeMid)
              kmFit[[reg]]$survfit <- res$survfit
              kmFit[[reg]]$survdiff <- res$survdiff
              kmTable <- rbind(kmTable,res$kmTable)
              if(verbose) setTxtProgressBar(pb, which(reglist == reg))
            }
            if(verbose) close(pb)
            rownames(kmTable) <- reglist
            kmTable <- data.frame(Regulons=reglist, kmTable, stringsAsFactors = FALSE)
            #---
            kmTable$Adjusted.Pvalue <- p.adjust(kmTable$Pvalue, method = pAdjustMethod)
            kmTable <- kmTable[sort.list(kmTable[,"Pvalue"]),, drop=FALSE]
            #---
            resKM <- list(Table=kmTable, Fit=kmFit)
            tns <- tns.set(tns, resKM, "KM")
            return(tns)
          })


#' Kaplan-Meier plots for TNS class objects
#'
#' Makes a 2 or 3 panel plot for survival analysis. The first panel shows the
#' differential Enrichment score (dES) for all samples, ranked by dES 
#' in their sections. The second (optional) panel shows the status of other 
#' attributes which may be present in the survival data frame for all samples. 
#' The third panel shows a Kaplan-Meier plot computed for the given survival 
#' data, with a curve for each section.
#'
#' @param tns A \linkS4class{TNS} object, which must have passed GSEA2 analysis.
#' @param regs An optional string vector specifying regulons to make the plot.
#' @param attribs A numeric vector. Contains the columns of the survival 
#' data.frame which will be plotted for the second panel.
#' @param fname A string. The name of the file in which the plot will be saved
#' @param fpath A string. The path to the directory where the plot will be saved
#' @param xlab A string. The label for the x axis on the third panel. This should
#' be the measure of time shown in the survival data frame after the last 
#' check-up.
#' @param ylab A string. The label for the y axis on the third panel
#' @param colorPalette A string, which can be 'red', 'blue', 'redblue', or 'bluered'. 
#' Alternatively, it can be colors or hex values.
#' @param plotpdf A logical value. If TRUE, the plot is saved as a pdf file. 
#' If false, it is plotted in the plotting area.
#' @param plotbatch A logical value. If TRUE, plots for all regulons are saved in 
#' the same file. If FALSE, each plot for each regulon is saved in a different file.
#' @param width A numeric value. Represents the width of the plot.
#' @param height A numeric value. Represents the height of the plot.
#' @param panelWidths A numeric vector of length=3 specifying the relative 
#' width of the internal panels.
#' 
#' @return A plot, showing a graphical analysis for the 'tnsKM' function.
#' @examples
#' # load survival data
#' data(survival.data)
#' 
#' # load TNI-object
#' data(stni, package = "RTN")
#'
#' stns <- tni2tnsPreprocess(stni, survivalData = survival.data, 
#'         keycovar = c('Grade','Age'), time = 1, event = 2)
#' stns <- tnsGSEA2(stns)
#' stns <- tnsKM(stns)
#' tnsPlotKM(stns)
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom survival survdiff survfit coxph Surv
#' @docType methods
#' @rdname tnsPlotKM-methods
#' @aliases tnsPlotKM
#' @export
#' 
setMethod("tnsPlotKM", "TNS", 
          function(tns, regs = NULL, attribs = NULL, 
                   fname = "survplot", fpath = ".", xlab = "Months", 
                   ylab = "Survival probability", colorPalette = "bluered", 
                   plotpdf = FALSE, plotbatch = FALSE, width = 6.3, 
                   height = 3.6, panelWidths = c(3, 2, 4)){
            #-- checks
            .tns.checks(tns, type = "Activity")
            .tns.checks(regs, type = "regs")
            .tns.checks(attribs, tns@survivalData, type = "attribs")
            .tns.checks(fname, type = "fname")
            .tns.checks(fpath, type = "fpath")
            .tns.checks(xlab, type = "xlab")
            .tns.checks(ylab, type = "ylab")
            .tns.checks(plotpdf, type = "plotpdf")
            .tns.checks(plotbatch, type = "plotbatch")
            .tns.checks(width, type = "width")
            .tns.checks(height, type = "height")
            .tns.checks(panelWidths, type = "panelWidths")
            tnstatus <- tnsGet(tns, what = "status")
            if(tnstatus["KM"] != "[x]")
              stop("NOTE: TNS object needs to be evaluated by 'tnsKM'!", 
                   call. = FALSE)
            
            #-- get data and para
            regulonActivity <- tnsGet(tns, what = "regulonActivity")
            survData <- tnsGet(tns, what = "survivalData")
            kmTable <- tnsGet(tns, what = "kmTable")
            kmFit <- tnsGet(tns, what = "kmFit")
            para <- tnsGet(tns, what = "para")
            endpoint <- para$endpoint
            excludeMid <- para$excludeMid
            nSections <-  para$nSections
            
            #-- check colorPalette with nSections
            .tns.checks(colorPalette, nSections, type = "colorPalette")
            
            #-- set endpoint
            survData$event[survData$time > endpoint] <- 0
            survData$time[survData$time > endpoint] <- endpoint
            
            #-- get attribs
            if (!is.null(attribs)) {
              if (is.list(attribs)) {
                groups <- unlist(lapply(attribs, length))
                idx <- unlist(attribs)
                attribs <- as.matrix(survData[, idx])
              } else {
                groups <- NULL
                attribs <- as.matrix(survData[, attribs])
              }
              if (!all(attribs %in% c(0, 1, NA))) 
                stop("'attribs' variables should only include binary values!")
            }
            
            #-- making reglist
            reglist <- kmTable$Regulons
            if(!is.null(regs)) {
              if (!all(regs %in% reglist)) {
                stop("all names in 'regs' should be listed in the slot 'results$regulonActivity' of the 'tns' object!")
              }
              reglist <- regs
            }
            
            #--- remove invalid string
            fname <- gsub(".pdf", '',fname, ignore.case = TRUE)
            
            #---plot
            if(plotbatch & plotpdf){
              fname <- paste(fname, ".pdf", sep = "")
              pdf(file = paste(fpath, "/", fname, sep = ""), width = width, height = height)
              for(reg in reglist){
                kmtb <- kmTable[reg,]
                survdf <- kmFit[[reg]]$survdiff
                survft <- kmFit[[reg]]$survfit
                .survplot(regulonActivity, kmtb=kmtb, survdf=survdf, survft=survft, reg=reg, 
                          endpoint=endpoint, xlab=xlab, ylab=ylab, colorPalette=colorPalette, 
                          panelWidths=panelWidths, excludeMid=excludeMid, attribs=attribs, 
                          groups=groups)
              }
              dev.off()
              if(plotpdf){
                tp1 <- c("NOTE: file '",fname,"' should be available either in the\n")
                tp2 <- c("working directory or in a user's custom directory!\n")
              }
            } else {
              for(reg in reglist){
                if(plotpdf){
                  pdf(file = paste(fpath, "/", fname, "_", reg, ".pdf", sep = ""), 
                      width = width, height = height)
                }
                kmtb <- kmTable[reg,]
                survdf <- kmFit[[reg]]$survdiff
                survft <- kmFit[[reg]]$survfit
                .survplot(regulonActivity, kmtb=kmtb, survdf=survdf, survft=survft, reg=reg, 
                          endpoint=endpoint, xlab=xlab, ylab=ylab, colorPalette=colorPalette, 
                          panelWidths=panelWidths, excludeMid=excludeMid, attribs=attribs, 
                          groups=groups)
                if (plotpdf) dev.off()
              }
              if(plotpdf){
                if(length(reglist)>1){
                  fname <- paste(fname, "...pdf", sep = "")
                  tp1 <- c("NOTE: ",length(reglist)," files named '",fname, 
                           "' should be available either in the\n")
                  tp2 <- c("working directory or in a user's custom directory!\n")
                } else {
                  fname <- paste(fname,"_",reglist, ".pdf", sep = "")
                  tp1 <- c("NOTE: file '",fname,
                           "' should be available either in the\n")
                  tp2 <- c("working directory or in a user's custom directory!\n")
                }
                message(tp1,tp2)
              }
            }
          })


#' Cox regression analysis for TNS class objects
#'
#' Run Cox multivariate regression for regulons and other covariates.
#'
#' @param tns A \linkS4class{TNS} object, which must have passed GSEA2 analysis.
#' @param regs An optional string vector listing regulons to be tested.
#' @param qqkeycovar A logical value. If TRUE, only the samples in the 2nd and 
#' 3rd quartils of 'keycovar' are used in the analysis. If FALSE, all samples
#' are used (see \code{\link{tni2tnsPreprocess}}).
#' @param verbose A logical value specifying to display detailed messages 
#' (when verbose=TRUE) or not (when verbose=FALSE).
#' @return Cox hazard models and statistics.
#' @examples
#' # load survival data
#' data(survival.data)
#' 
#' # load TNI-object
#' data(stni, package = "RTN")
#'
#' stns <- tni2tnsPreprocess(stni, survivalData = survival.data, 
#' keycovar = c('Age','Grade'), time = 1, event = 2)
#' stns <- tnsGSEA2(stns)
#' stns <- tnsCox(stns, regs = c('PTTG1','E2F2','FOXM1'))
#' tnsGet(stns, "coxTable")
#' 
#' @docType methods
#' @importFrom stats p.adjust
#' @rdname tnsCox-methods
#' @aliases tnsCox
#' @export
#' 
setMethod("tnsCox", "TNS", 
          function(tns, regs = NULL, qqkeycovar = FALSE, verbose = TRUE)
          {
            
            #-- checks
            .tns.checks(tns, type = "Activity")
            .tns.checks(regs, type = "regs")
            .tns.checks(qqkeycovar, type = "qqkeycovar")
            .tns.checks(verbose, type = "verbose")
            
            #-- gets
            regulonActivity <- tnsGet(tns, what = "regulonActivity")
            survData <- tnsGet(tns, what = "survivalData")
            .tns.checks(survData, type = "survival_cox")
            para <- tnsGet(tns, what = "para")
            keycovar <- para$keycovar
            excludeMid <- para$excludeMid
            pAdjustMethod <- para$pAdjustMethod
            endpoint <- para$endpoint
            
            #-- set endpoint
            survData$event[survData$time > endpoint] <- 0
            survData$time[survData$time > endpoint] <- endpoint
            
            #-- checks
            dif <- regulonActivity$dif
            if (excludeMid) {
              dif[regulonActivity$regstatus == regulonActivity$center] <- NA
            }
            
            if(!is.null(regs)) {
              if (!all(regs %in% colnames(dif))) {
                stop("Not all 'regs' have regulonActivity!")
              }
              idx <- colnames(dif) %in% regs
              dif <- dif[, idx]
              dif <- dif[, regs]
            }
            
            #---set valid names for a 'formula'
            regs <- colnames(dif)
            xregs <- .namesCorrect(regs)
            colnames(dif) <- xregs
            names(regs) <- xregs
            
            #---combine dif and survivalData 
            dtsumm <- cbind(survData[rownames(dif), ], dif)
            
            #--- set keycovar by quantile
            if(qqkeycovar && !is.null(keycovar)){
              for(kvar in keycovar){
                tp <- dtsumm[[kvar]]
                if(is.numeric(tp)){
                  ql <- quantile(tp, c(0.25, 0.75), na.rm = TRUE)
                  tp[tp < ql[1]] <- NA
                  tp[tp > ql[2]] <- NA
                  dtsumm[[kvar]] <- tp
                } 
              }
            }
            
            #--- filter data
            dtsumm <- dtsumm[, c("time", "event", keycovar, xregs)]
            
            #--- get cox formula
            if(is.null(keycovar)){
              fm1 <- "Surv(time, event) ~ "
            } else {
              fm1 <- paste("Surv(time, event) ~ ", 
                           paste(keycovar, collapse = "+"), 
                           sep = "")
            }
            
            if (verbose) {
              cat("Computing Cox regression models...\n")
              pb <- txtProgressBar(min = 0, max = length(xregs), style = 3)
            }
            
            #--- fit cox regression models
            coxFit <- lapply(xregs, function(rg){
              if(verbose) setTxtProgressBar(pb, which(xregs == rg))
              nas <- is.na(dtsumm[, rg])
              if (sum(nas) > nrow(dtsumm)/2) {
                cxmd <- NA
              } else {
                if(!is.null(keycovar)){
                  fm2 <- formula(paste(fm1, rg, sep = "+"))
                  cxmd <- coxph(fm2, data = dtsumm[!nas, ])
                } else {
                  fm2 <- formula(paste(fm1, rg, sep = ""))
                  cxmd <- coxph(fm2, data = dtsumm[!nas, ])
                }
              }
              cxmd
            })
            if(verbose) close(pb)
            names(coxFit) <- xregs
            
            #--- get probs
            coxprobs <- sapply(xregs, function(rg) {
              cxmd <- coxFit[[rg]]
              if(class(cxmd)=="coxph"){
                res <- summary(cxmd)
                res <- res$coefficients[rg,c("Pr(>|z|)")]
              } else {
                res <- 1
              }
              res
            })
            ci <- p.threshold(pvals=coxprobs, alpha=0.05, method=pAdjustMethod)
            
            #--- get coefs
            coxcoefs <- sapply(xregs, function(rg) {
              cxmd <- coxFit[[rg]]
              if(class(cxmd)=="coxph"){
                res <- summary(cxmd, conf.int = 1-ci)
                res <- res$conf.int[rg,c(1,3:4)]
              } else {
                res <- c(1, 0.99, 1.01)
              }
              res
            })
            coxcoefs <- t(coxcoefs)
            coxTable <- cbind(coxcoefs,coxprobs)
            
            #--- add keycovars to coxTable
            if(!is.null(keycovar)){
              idx <- which.max(coxTable[, 1])
              rg <- rownames(coxTable)[idx]
              cxmd <- coxFit[[rg]]
              resref <- summary(cxmd, conf.int = 1-ci)
              ci <- resref$conf.int[,c(1,3:4)]
              pr <- resref$coefficients[,c("Pr(>|z|)")]
              resref <- cbind(ci,pr)
              coxTable <- rbind(resref[-nrow(resref),], coxTable)
              rownames(coxTable)[1:length(keycovar)] <- keycovar
            }
            colnames(coxTable) <- c("HR","Lower95","Upper95","Pvalue")
            
            #--- add original symbols to rownames
            idx <- match(names(regs), rownames(coxTable))
            rownames(coxTable)[idx] <- regs
            coxTable <- data.frame(Regulons=rownames(coxTable), coxTable)
            
            #--- adjust pvalues and assign significant results
            coxTable$Adjusted.Pvalue <- p.adjust(coxTable$Pvalue, method = pAdjustMethod)
            coxTable <- coxTable[sort.list(coxTable[,"Pvalue"]),, drop=FALSE]
            
            #--- return
            res <- list(Table=coxTable, Fit=coxFit)
            tns <- tns.set(tns, res, "Cox")
            return(tns)
            
          })



#' Cox plots for TNS class objects
#'
#' Run Cox multivariate regression for regulons and key covariables.
#'
#' @param tns A \linkS4class{TNS} object, which must have passed GSEA2 analysis.
#' @param regs An optional string vector specifying regulons to make the plot.
#' @param fname A string. The name of the PDF file which will contain the plot.
#' @param fpath A string. The directory where the file will be saved.
#' @param ylab A string. The label of the y-axis, describing what is represented.
#' @param xlab A string. The label of the x-axis.
#' @param width A numeric value. The width of the plot.
#' @param height A numeric value. The height of the plot.
#' @param xlim A vector with 2 values indicating lowest and highest HR values.
#' @param sortregs A logical value. If TRUE, regulons are sorted from most 
#' negatively associated with hazard to most positively associated with hazard.
#' @param plotpdf A logical value.
#' @return A Cox hazard model plot and statistics.
#' @examples
#' # load survival data
#' data(survival.data)
#' 
#' # load TNI-object
#' data(stni, package = "RTN")
#'
#' stns <- tni2tnsPreprocess(stni, survivalData = survival.data, 
#' keycovar = c('Age','Grade'), time = 1, event = 2)
#' stns <- tnsGSEA2(stns)
#' stns <- tnsCox(stns, regs = c('PTTG1','E2F2','FOXM1'))
#' tnsPlotCox(stns)
#' 
#' @docType methods
#' @importFrom stats p.adjust
#' @rdname tnsPlotCox-methods
#' @aliases tnsPlotCox
#' @export
#' 
setMethod("tnsPlotCox", "TNS", 
          function(tns, regs = NULL, fname = "coxplot", 
                   fpath = ".", ylab = "Regulons and other covariates", 
                   xlab = "Hazard Ratio (95% CI)", width = 5, height = 5, 
                   xlim = c(0.3, 3), sortregs = TRUE, plotpdf = FALSE){
            
            #-- checks
            .tns.checks(tns, type = "Activity")
            .tns.checks(regs, type = "regs")
            .tns.checks(fname, type = "fname")
            .tns.checks(fpath, type = "fpath")
            .tns.checks(ylab, type = "ylab")
            .tns.checks(xlab, type = "xlab")
            .tns.checks(width, type = "width")
            .tns.checks(height, type = "height")
            .tns.checks(xlim, type = "xlim_log")
            .tns.checks(sortregs, type = "sortregs")
            .tns.checks(plotpdf, type = "plotpdf")
            tnstatus <- tnsGet(tns, what = "status")
            if(tnstatus["Cox"] != "[x]")
              stop("NOTE: TNS object needs to be evaluated by 'tnsCox'!", 
                   call. = FALSE)
            
            #-- gets
            coxTable <- tnsGet(tns, what = "coxTable")
            para <- tnsGet(tns, what = "para")
            keycovar <- para$keycovar
            
            if(!is.null(regs)){
              if (any(regs %in% keycovar)) {
                stop("'regs' should not be listed in 'keycovar'!")
              }
              if (!all(regs %in% coxTable$Regulons)) {
                stop("All 'regs' should be listed in the 'coxTable', please see 'tnsGet' function!")
              }
              idx <- coxTable$Regulons %in% c(keycovar,regs)
              coxTable <- coxTable[idx, ]
            } else {
              regs <- setdiff(coxTable$Regulons,keycovar)
            }
            
            #--- sort coxTable
            coxTable <- coxTable[sort.list(coxTable$HR, decreasing = T),]
            idx <- which(coxTable$Regulons%in%keycovar)
            if(length(idx)>0)
              coxTable <- rbind(coxTable[-idx,],coxTable[idx,])
            
            #--- plot
            .plotCox(coxTable, regs = regs, keycovar = keycovar, fpath = fpath, 
                     fname = fname, width = width, height = height, xlim = xlim, 
                     xlab = xlab, ylab = ylab, plotpdf = plotpdf)
            
          })


setMethod("show", "TNS", function(object) {
  message("a TNS (Transcriptional Network - Survival) object:\n")
  message("--status:")
  print(object@status, quote = FALSE)
})

#' Get information from slots in a TNS object
#'
#'Get information from individual slots in a TNS object and any
#'available results from a previous analysis.
#'
#' @param tns A \linkS4class{TNS} object.
#' @param what A string specifying what should be retrieved from the object.
#' Options: 'status','survivalData', 'regulonActivity', 'TNI', 'para', 'kmTable', 'kmFit', 
#' 'coxTable', 'coxFit', 'kmInteractionTable', 'kmInteractionFit', 
#' 'coxInteractionTable', 'coxInteractionFit', and 'regulatoryElements'.
#' @return Content from slots in the \linkS4class{TNS} object.
#' @examples
#' # load survival data
#' data(survival.data)
#' 
#' # load TNI-object
#' data(stni, package = "RTN")
#'
#' stns <- tni2tnsPreprocess(stni, survivalData = survival.data, 
#' keycovar = c('Grade','Age'), time = 1, event = 2)
#' stns <- tnsGSEA2(stns)
#' regulonActivity <- tnsGet(stns, 'regulonActivity')
#'
#' @docType methods
#' @rdname tnsGet-methods
#' @aliases tnsGet
#' @export

setMethod("tnsGet", "TNS", function(tns, what)
{
  if (what == "survivalData"){
    return(tns@survivalData) 
  } else if(what == "regulonActivity"){
    return(tns@results$regulonActivity) 
  } else if(what == "TNI"){
    return(tns@TNI)
  } else if(what == "kmTable"){
    return(tns@results$KM$Table)
  } else if (what == "kmFit"){
    return(tns@results$KM$Fit)
  } else if(what == "coxTable"){
    return(tns@results$Cox$Table)
  } else if(what == "coxFit"){
    return(tns@results$Cox$Fit)
  } else if(what == "kmInteractionTable"){
    return(tns@results$KmInt$Table)
  } else if(what == "kmInteractionFit"){
    return(tns@results$KmInt$Fit)
  } else if(what == "coxInteractionTable"){
    return(tns@results$CoxInt$Table)
  } else if(what == "coxInteractionFit"){
    return(tns@results$CoxInt$Fit)
  } else if (what == "para"){
    return(tns@para)
  } else if(what == "regulatoryElements"){
    return(tns@TNI@regulatoryElements)
  } else if(what == "status"){
    return(tns@status)
  } else {
    stop("'what' must be one of:\n",
         "'status', 'survivalData', 'regulonActivity', 'TNI', 'para','kmTable', 'kmFit',\n",
         "'coxTable', 'coxFit', 'kmInteractionTable', 'kmInteractionFit',\n",
         "'coxInteractionTable', 'coxInteractionFit', and 'regulatoryElements'.")
  }
})

#' Survival analysis for dual regulons
#'
#' A generic call to 'tnsCoxInteraction' and 'tnsKmInteraction' functions.
#' 
#' @param tns A \linkS4class{TNS} object, which must have passed GSEA2 analysis.
#' @param ... Parameters passed to \code{\link{tnsKmInteraction}} and 
#' \code{\link{tnsCoxInteraction}} functions.
#' @param verbose A logical value specifying to display detailed messages 
#' (when verbose=TRUE) or not (when verbose=FALSE).
#' 
#' @return A \linkS4class{TNS} object evaluated by the 'tnsKmInteraction' and 
#' 'tnsCoxInteraction' functions.
#' @examples
#' # load survival data
#' data(survival.data)
#' 
#' # load TNI-object
#' data(stni, package = "RTN")
#'
#' stns <- tni2tnsPreprocess(stni, survivalData = survival.data, 
#' keycovar = c('Grade','Age'), time = 1, event = 2)
#' stns <- tnsGSEA2(stns)
#'
#' # infer dual regulons
#' smbr <- tni2mbrPreprocess(stni)
#' smbr <- mbrAssociation(smbr)
#'
#' # survival analysis for dual regulons
#' # stns <- tnsInteraction(stns, smbr)
#' 
#' @importClassesFrom RTNduals MBR
#' @importFrom RTNduals mbrGet
#' @importFrom survival survdiff survfit coxph Surv
#' @importFrom stats complete.cases
#' @docType methods
#' @rdname tnsInteraction-methods
#' @aliases tnsInteraction
#' @export
#' 
setMethod("tnsInteraction", "TNS", 
          function(tns, ..., verbose=TRUE){
            if (verbose)
              cat("Assessing interaction between regulons...\n")
            tns <- tnsKmInteraction(tns, ...=...)
            tns <- tnsCoxInteraction(tns, ...=...)
            return(tns)
          })


#' Kaplan-Meier analysis for dual regulons
#'
#' Kaplan-Meier analysis for dual regulons, assessing the interaction between regulons.
#' 
#' @param tns A \linkS4class{TNS} object, which must have passed GSEA2 analysis.
#' @param mbr An \linkS4class{MBR} object computed from the same 'TNI' object used
#' to make the 'TNS' object.
#' @param stepFilter A single logical value specifying to use a step-filter 
#' algorithm, testing dual regulons that have at least one significant predictor
#' in the 'tnsKM' method (when stepFilter=TRUE) or not (when stepFilter=FALSE).
#' @param pValueCutoff An numeric value. The p-value cutoff applied to the results
#' from the previous steps of the analysis pipeline (when stepFilter=TRUE).
#' @param verbose A logical value specifying to display detailed messages 
#' (when verbose=TRUE) or not (when verbose=FALSE).
#' 
#' @return Results from 'survfit' and 'survdiff', including log-rank statistics.
#' @examples
#' # load survival data
#' data(survival.data)
#' 
#' # load TNI-object
#' data(stni, package = "RTN")
#'
#' stns <- tni2tnsPreprocess(stni, survivalData = survival.data, 
#' keycovar = c('Grade','Age'), time = 1, event = 2)
#' stns <- tnsGSEA2(stns)
#'
#' # infer dual regulons
#' smbr <- tni2mbrPreprocess(stni)
#' smbr <- mbrAssociation(smbr)
#'
#' # KM analysis for dual regulons
#' # stns <- tnsKmInteraction(stns, smbr, stepFilter = FALSE)
#' # tnsGet(stns, "kmInteractionTable")
#' 
#' @importClassesFrom RTNduals MBR
#' @importFrom RTNduals mbrGet
#' @importFrom survival survdiff survfit coxph Surv
#' @importFrom stats complete.cases
#' @docType methods
#' @rdname tnsKmInteraction-methods
#' @aliases tnsKmInteraction
#' @export
#' 
setMethod("tnsKmInteraction", "TNS", 
          function(tns, mbr, stepFilter = TRUE, pValueCutoff = 0.05, verbose = TRUE){
            
            #-- checks
            .tns.checks(tns, type = "Activity")
            .tns.checks(mbr, type = "MBR")
            .tns.checks(stepFilter, type = "stepFilter")
            .tns.checks(pValueCutoff, type = "pValueCutoff")
            .tns.checks(verbose, type = "verbose")
            
            #-- check TNIs
            if(!identical(tnsGet(tns, what = "TNI"), mbrGet(mbr, what = "TNI"))){
              stop("NOTE: slots 'TNI' must be identical between 'tns' and 'mbr' objects!")
            }
            
            #-- stratification
            tns <- tnsStratification(tns, nSections = 1, center = TRUE)
            
            #-- get data
            regulonActivity <- tnsGet(tns, what = "regulonActivity")
            survData <- tnsGet(tns, what = "survivalData")
            para <- tnsGet(tns, what = "para")
            endpoint <- para$endpoint
            excludeMid <- para$excludeMid
            pAdjustMethod <- para$pAdjustMethod
            
            #-- set endpoint
            survData$event[survData$time > endpoint] <- 0
            survData$time[survData$time > endpoint] <- endpoint
            
            #-- get dualregs
            dualregs <- mbrGet(mbr, what = "dualRegulons")
            dualtb <- t(sapply(dualregs, function(dual){
              unlist(strsplit(dual, "~"))
            }))
            colnames(dualtb) <- c("reg1","reg2")
            dualtb <- data.frame(dualtb, stringsAsFactors = FALSE)
            
            #-- apply stepFilter from previous methods
            if(stepFilter){
              tnstatus <- tnsGet(tns, what = "status")
              if(tnstatus["KM"] != "[x]"){
                stop("NOTE: when 'stepFilter=TRUE', the TNS object needs to be evaluated by 'tnsKM'!", 
                     call. = FALSE)
              }
              mktb <- tnsGet(tns, what = "kmTable")
              mktb <- mktb[mktb$Adjusted.Pvalue<pValueCutoff, ]
              idx <- dualtb$reg1%in%rownames(mktb) | dualtb$reg2%in%rownames(mktb)
              dualtb <- dualtb[idx,]
              if(nrow(dualtb)==0){
                message("NOTE: using 'stepFilter=TRUE': no significant regulon from the 'tnsKM' method!",
                        call. = FALSE)
                return(tns)
              }
              dualregs <- rownames(dualtb)
            }
            
            if (verbose) {
              cat("Computing survival curves...\n")
              pb <- txtProgressBar(min = 0, max = length(dualregs), style = 3)
            }
            
            #--- logrank
            kmFit <- list()
            kmTable <- NULL
            for(dual in dualregs){
              regs <- as.character(dualtb[dual,])
              res <- .survstatsDuals(regulonActivity, survData=survData, regs=regs, excludeMid=excludeMid)
              kmFit[[dual]]$survfit <- res$survfit
              kmFit[[dual]]$survdiff <- res$survdiff
              kmTable <- rbind(kmTable,res$kmTable)
              if(verbose) setTxtProgressBar(pb, which(dualregs == dual))
            }
            if(verbose) close(pb)
            kmTable <- data.frame(Regulon1=dualtb$reg1, Regulon2=dualtb$reg2,
                                  kmTable, stringsAsFactors = FALSE)
            rownames(kmTable) <- dualregs
            #---
            kmTable$Adjusted.Pvalue <- p.adjust(kmTable$Pvalue, method = pAdjustMethod)
            kmTable <- kmTable[sort.list(kmTable[,"Pvalue"]),, drop=FALSE]
            #---
            resKM <- list(Table=kmTable, Fit=kmFit)
            tns <- tns.set(tns, resKM, "KmInt")
            return(tns)
            
          })


#' Plot results from Kaplan-Meier analysis for dual regulons
#'
#'
#' @param tns A \linkS4class{TNS} object, which must have passed GSEA2 analysis.
#' @param dualreg A character string with the name of a dual regulon.
#' @param fname A string. The name of the file in which the plot will be saved
#' @param fpath A string. The path to the directory where the plot will be saved
#' @param xlab A string. The label for the x axis on the third panel. This should
#' be the measure of time shown in the survival data.frame after the last 
#' check-up.
#' @param ylab A string. The label for the y axis on the third panel
#' @param colorPalette A string, which can be 'red', 'blue', 'redblue', or 'bluered'. 
#' Alternatively, it can be a vector of five colors or hex values.
#' @param width A numeric value. Represents the width of the plot.
#' @param height A numeric value. Represents the height of the plot.
#' @param plotpdf A logical value. If TRUE, the plot is saved as a pdf file. 
#' If false, it is plotted in the plotting area.
#' 
#' @return  A plot, showing a graphical analysis for the 'tnsKmInteraction' function.
#' @examples
#' # load survival data
#' data(survival.data)
#' 
#' # load TNI-object
#' data(stni, package = "RTN")
#'
#' stns <- tni2tnsPreprocess(stni, survivalData = survival.data, 
#' keycovar = c('Grade','Age'), time = 1, event = 2)
#' stns <- tnsGSEA2(stns)
#' 
#' # infer dual regulons
#' smbr <- tni2mbrPreprocess(stni)
#' smbr <- mbrAssociation(smbr)
#' 
#' # KM analysis for dual regulons
#' # stns <- tnsKmInteraction(stns, smbr, stepFilter=FALSE)
#' # tnsPlotKmInteraction(stns, dualreg = "FOXM1~PTTG1")
#' 
#' @importFrom RColorBrewer brewer.pal
#' @importFrom survival survdiff survfit coxph Surv
#' @importFrom stats complete.cases
#' @docType methods
#' @rdname tnsPlotKmInteraction-methods
#' @aliases tnsPlotKmInteraction
#' @export
#' 
setMethod("tnsPlotKmInteraction", "TNS", 
          function(tns, dualreg, fname = "kmInteraction", 
                   fpath = ".", xlab = "Months", 
                   ylab = "Survival probability", 
                   colorPalette = "bluered", width = 4, 
                   height = 4, plotpdf = FALSE){
            
            #-- checks
            .tns.checks(tns, type = "Activity")
            .tns.checks(dualreg, type = "dualreg")
            .tns.checks(fname, type = "fname")
            .tns.checks(fpath, type = "fpath")
            .tns.checks(xlab, type = "xlab")
            .tns.checks(ylab, type = "ylab")
            .tns.checks(colorPalette, 2, type = "colorPalette")
            .tns.checks(width, type = "width")
            .tns.checks(height, type = "height")
            .tns.checks(plotpdf, type = "plotpdf")
            
            #-- stratification
            tns <- tnsStratification(tns, nSections = 1, center = TRUE)
            
            #-- get data
            regulonActivity <- tnsGet(tns, what = "regulonActivity")
            survData <- tnsGet(tns, what = "survivalData")
            para <- tnsGet(tns, what = "para")
            endpoint <- para$endpoint
            excludeMid <- para$excludeMid
            
            #-- set endpoint
            survData$event[survData$time > endpoint] <- 0
            survData$time[survData$time > endpoint] <- endpoint
            
            #-- get regs
            regs <- unlist(strsplit(dualreg, split = "~", fixed=TRUE))
            
            #-- making reglist
            if (!all(regs %in% colnames(regulonActivity$regstatus))) {
              stop("all names in 'dualreg' should be listed in the slot 'results$regulonActivity' of the 'tns' object!")
            }
            
            #---plot
            if(plotpdf){
              fname <- gsub(".pdf", '',fname, ignore.case = TRUE)
              fname <- paste(fname,"_",paste(regs,collapse = "_"),".pdf", sep = "")
              pdf(file = paste(fpath, "/", fname, sep = ""), width = width, height = height)
            }
            .survplotDuals(regulonActivity, survData=survData, regs=regs, 
                           endpoint=endpoint, excludeMid=excludeMid, 
                           ylab=ylab, xlab=xlab, colorPalette=colorPalette)
            if(plotpdf){
              dev.off()
              tp1 <- c("NOTE: file '",fname,
                       "' should be available either in the\n")
              tp2 <- c("working directory or in a user's custom directory!\n")
              message(tp1,tp2)
            }
            
          })


#' Cox regression analysis for dual regulons
#'
#' Cox regression analysis for dual regulons, including the interaction term.
#'
#' @param tns A \linkS4class{TNS} object with regulons used to compute the dual regulons.
#' @param mbr An \linkS4class{MBR} object computed from the same 'TNI' object used
#' to make the 'TNS' object.
#' @param stepFilter A single logical value specifying to use a step-filter 
#' algorithm, testing dual regulons that have at least one significant predictor
#' in the 'tnsCox' method (when stepFilter=TRUE) or not (when stepFilter=FALSE).
#' @param pValueCutoff An numeric value. The p-value cutoff applied to the results
#' from the previous steps of the analysis pipeline (when stepFilter=TRUE).
#' @param verbose A logical value specifying to display detailed messages 
#' (when verbose=TRUE) or not (when verbose=FALSE).
#' 
#' @return Cox hazard models and statistics.
#' @examples
#' # load survival data
#' data(survival.data)
#' 
#' # load TNI-object
#' data(stni, package = "RTN")
#' 
#' # perform survival analysis for regulons
#' stns <- tni2tnsPreprocess(stni, survivalData = survival.data, 
#' time = 1, event = 2)
#' stns <- tnsGSEA2(stns)
#' 
#' # infer dual regulons
#' smbr <- tni2mbrPreprocess(stni)
#' smbr <- mbrAssociation(smbr)
#'
#' # run Cox regression for dual regulons
#' stns <- tnsCoxInteraction(stns, smbr, stepFilter = FALSE)
#' tnsGet(stns, "coxInteractionTable")
#'
#' @seealso \code{\link{tni2mbrPreprocess}} for all plot parameters
#' @return An updated TNS-class object containing Cox regression models 
#' for all given duals
#' @importClassesFrom RTNduals MBR
#' @importFrom RTNduals mbrGet
#' @docType methods
#' @rdname tnsCoxInteraction-methods
#' @aliases tnsCoxInteraction
#' @export
#' 
setMethod("tnsCoxInteraction", "TNS",
          function(tns, mbr, stepFilter = TRUE, pValueCutoff = 0.05, 
                   verbose = TRUE)
          {
            
            #-- checks
            .tns.checks(tns, type = "Activity")
            .tns.checks(mbr, type = "MBR")
            .tns.checks(stepFilter, type = "stepFilter")
            .tns.checks(pValueCutoff, type = "pValueCutoff")
            .tns.checks(verbose, type = "verbose")
            
            #-- check TNIs
            if(!identical(tnsGet(tns, what = "TNI"), mbrGet(mbr, what = "TNI"))){
              stop("NOTE: slots 'TNI' must be identical between 'tns' and 'mbr' objects!")
            }
            
            #-- get dualregs
            dualregs <- mbrGet(mbr, what = "dualRegulons")
            dualtb <- t(sapply(dualregs, function(dual){
              unlist(strsplit(dual, "~"))
            }))
            colnames(dualtb) <- c("reg1","reg2")
            dualtb <- data.frame(dualtb, stringsAsFactors = FALSE)
            
            #-- apply stepFilter from previous methods
            if(stepFilter){
              tnstatus <- tnsGet(tns, what = "status")
              if(tnstatus["Cox"] != "[x]"){
                stop("NOTE: when 'stepFilter=TRUE', the TNS object needs to be evaluated by 'tnsCox'!", 
                     call. = FALSE)
              }
              coxtb <- tnsGet(tns, what = "coxTable")
              coxtb <- coxtb[coxtb$Adjusted.Pvalue<=pValueCutoff, ]
              idx <- dualtb$reg1%in%rownames(coxtb) | dualtb$reg2%in%rownames(coxtb)
              dualtb <- dualtb[idx,]
              if(nrow(dualtb)==0){
                message("NOTE: using 'stepFilter=TRUE': no significant regulon from the 'tnsCox' method!",
                        call. = FALSE)
                return(tns)
              }
              dualregs <- rownames(dualtb)
            }
            
            #-- gets
            regulonActivity <- tnsGet(tns, what = "regulonActivity")
            survData <- tnsGet(tns, what = "survivalData")
            .tns.checks(survData, type = "survival_cox")
            para <- tnsGet(tns, what = "para")
            excludeMid <- para$excludeMid
            pAdjustMethod <- para$pAdjustMethod
            endpoint <- para$endpoint
            
            #-- set endpoint
            survData$event[survData$time > endpoint] <- 0
            survData$time[survData$time > endpoint] <- endpoint
            
            #-- mbr vs tns checks
            tns.regs <- colnames(regulonActivity$dif)
            if (!all(dualtb$reg1 %in% tns.regs) | !all(dualtb$reg2 %in% tns.regs)) {
              idx1 <- which(dualtb$reg1 %in% tns.regs)
              idx2 <- which(dualtb$reg2 %in% tns.regs)
              dualtb <- dualtb[intersect(idx1, idx2),,drop=FALSE]
              dualregs <- rownames(dualtb)
              tp <- paste("Not all regulons have had enrichment scores computed.\n",
                          "This is possibly due to 'minRegulonSize'. Regression being computed for only ", 
                          length(dualregs), " regulon pairs.", sep="")
              warning(tp, call. = FALSE)
            }
            
            #-- checks
            dif <- regulonActivity$dif
            if (excludeMid) {
              dif1[regulonActivity$regstatus == regulonActivity$center] <- NA
            }
            
            #-- correct names for a 'formula'
            colnames(dif) <- .namesCorrect(colnames(dif))
            
            #-- combine data frames
            dtsumm <- cbind(survData, dif)
            
            if (verbose) {
              cat("Computing Cox regression models...\n")
              pb <- txtProgressBar(min = 0, max = length(dualregs), style = 3)
            }
            
            #--- fit cox regression model
            coxFit <- lapply(dualregs, function(dual){
              if(verbose) setTxtProgressBar(pb, which(dualregs == dual))
              regs <- unlist(strsplit(dual, "~"))
              regs <- .namesCorrect(regs)
              nas1 <- is.na(dif[, regs[1]])
              nas2 <- is.na(dif[, regs[2]])
              if (sum(nas1) > nrow(dif)/2 | sum(nas2) > nrow(dif)/2){
                cxmd <- NA
              } else {
                coxdual <- paste(regs[1], regs[2], sep = "*")
                fm <- formula(paste("Surv(time, event) ~ ", coxdual, sep = ""))
                cxmd <- coxph(fm, data = dtsumm[!(nas1&nas2), ])
              }
              cxmd
            })
            if(verbose) close(pb)
            names(coxFit) <- dualregs
            
            #--- get probs
            coxprobs <- sapply(dualregs, function(dual){
              cxmd <- coxFit[[dual]]
              if(class(cxmd)=="coxph"){
                res <- summary(cxmd)
                res <- res$coefficients[3,c("Pr(>|z|)")]
              } else {
                res <- 1
              }
              res
            })
            ci <- p.threshold(pvals=coxprobs, alpha=0.05, method=pAdjustMethod)
            
            #--- get coefs
            coxcoefs <- sapply(dualregs, function(dual){
              cxmd <- coxFit[[dual]]
              if(class(cxmd)=="coxph"){
                res <- summary(cxmd, conf.int = 1-ci)
                res <- res$conf.int[3,c(1,3:4)]
              } else {
                res <- c(1, 0.99, 1.01)
              }
              res
            })
            coxcoefs <- t(coxcoefs)
            coxTable <- cbind(coxcoefs,coxprobs)
            #---
            nms <- t(sapply(dualregs, function(dual){
              unlist(strsplit(dual, "~"))
            }))
            #---
            coxTable <- data.frame(nms, coxTable, stringsAsFactors = FALSE)
            colnames(coxTable) <- c("Regulon1","Regulon2", "HR", "Lower95", "Upper95", "Pvalue")
            coxTable$Adjusted.Pvalue <- p.adjust(coxTable$Pvalue, method = pAdjustMethod)
            coxTable <- coxTable[sort.list(coxTable[,"Pvalue"]),, drop=FALSE]
            res <- list(Table=coxTable, Fit=coxFit)
            tns <- tns.set(tns, res, "CoxInt")
            return(tns)
            
          })

#' Plot results from Cox regression analysis for dual regulons
#'
#' @param tns A `TNS` object with regulons used to compute the dual regulon.
#' @param dualreg A character string with the name of a dual regulon.
#' @param xlim A numeric vector of length 2, i.e. xlim = c(x1, x2),
#' indicating the limits of the plot for the first member of the dual regulon. 
#' If xlim = NULL, it will be derevided from the observed data ranges. 
#' Values must be in the range [-2,2].
#' @param ylim A numeric vector of length 2, i.e. ylim = c(y1, y2),
#' indicating the limits of the plot for the second member of the dual regulon. 
#' If ylim = NULL, it will be derevided from the observed data ranges. 
#' Values must be in the range [-2,2]. If plotype='2D', ylim represents the 
#' two fixed values for the second member of the dual regulon.
#' @param hlim A numeric vector of length 2, i.e. hlim = c(h1, h2), 
#' indicating the limits of the plot for the Hazard Ratio (HR).
#' If hlim = NULL, it will be derevided from the observed data ranges.
#' If plotype='2D', HR is represented in the y-axis.
#' @param hcols A vector of length 2 indicating a diverging color scheme for 
#' the Hazard Ratio (HR).
#' @param showdata A logical value indicating whether to show the original data 
#' used to fit linear model.
#' @param colorPalette A string, which can be 'red', 'blue', 'redblue', or 'bluered'. 
#' Alternatively, it can be a vector of five colors or hex values.
#' @param fname A string. The name of the PDF file (when plotpdf=TRUE).
#' @param fpath A string. The directory where the file will be saved.
#' @param width A numeric value. The width of the plot.
#' @param height A numeric value. The height of the plot.
#' @param plotype A string indicating '2D' of '3D' plot type. If plotype = '2D', 
#' the Hazard Ratio is represented in the y-axis.
#' @param plotpdf A logical value.
#' @return A Cox hazard model plot and statistics.
#' @examples
#' # load survival data
#' data(survival.data)
#' 
#' # load TNI-object
#' data(stni, package = "RTN")
#' 
#' # perform survival analysis for regulons
#' stns <- tni2tnsPreprocess(stni, survivalData = survival.data, time = 1, event = 2)
#' stns <- tnsGSEA2(stns, verbose=FALSE)
#' 
#' # infer dual regulons
#' smbr <- tni2mbrPreprocess(stni)
#' smbr <- mbrAssociation(smbr)
#'
#' # run Cox regression for dual regulons
#' # stns <- tnsCoxInteraction(stns, smbr, stepFilter = FALSE)
#' # tnsPlotCoxInteraction(stns, dualreg = "FOXM1~PTTG1")

#' @seealso \code{\link{tnsKM}},
#' \code{\link{tnsCox}}
#' @return A 3D heatmap plot.
#' @importClassesFrom RTNduals MBR 
#' @importFrom RTNduals mbrPlotInteraction
#' @importFrom survival coxph
#' @importFrom grDevices col2rgb
#' @importFrom graphics title
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats model.frame
#' @docType methods
#' @rdname tnsPlotCoxInteraction-methods
#' @aliases tnsPlotCoxInteraction
#' @export
#'
setMethod("tnsPlotCoxInteraction", "TNS", 
          function(tns, dualreg, xlim=NULL, ylim=NULL, hlim=NULL, 
                   hcols = c("#008080ff","#d45500ff"), showdata = TRUE, 
                   colorPalette = "bluered", fname = "coxInteraction", 
                   fpath = ".", width = 4.5, height = 4, plotype = "3D", 
                   plotpdf = FALSE){
            
            #-- checks
            .tns.checks(tns, type = "Activity")
            .tns.checks(dualreg, type = "dualreg")
            if(!is.null(xlim)) .tns.checks(xlim, type = "xlim_reg")
            if(!is.null(ylim)) .tns.checks(ylim, type = "ylim_reg")
            if(!is.null(hlim)) .tns.checks(hlim, type = "hlim_log")
            .tns.checks(hcols, type = "hcols")
            .tns.checks(showdata, type = "showdata")
            .tns.checks(colorPalette, 2, type = "colorPalette")
            .tns.checks(fname, type = "fname")
            .tns.checks(fpath, type = "fpath")
            .tns.checks(width, type = "width")
            .tns.checks(height, type = "height")
            .tns.checks(plotype, type = "plotype")
            .tns.checks(plotpdf, type = "plotpdf")
            tnstatus <- tnsGet(tns, what = "status")
            if(tnstatus["CoxInt"] != "[x]")
              stop("NOTE: TNS object needs to be evaluated by 'tnsCoxInteraction'!", 
                   call. = FALSE)
            
            #-- get regs
            regs <- unlist(strsplit(dualreg, split = "~", fixed=TRUE))
            
            #-- get cox model
            coxInteractionFit <- tnsGet(tns, "coxInteractionFit")
            coxInteractionTable <- tnsGet(tns, "coxInteractionTable")
            if(dualreg%in%names(coxInteractionFit)){
              model <- coxInteractionFit[[dualreg]]
              model$pAdjustInteraction <- coxInteractionTable[dualreg,"Adjusted.Pvalue"]
            } else {
              tp <- paste(regs[2:1], collapse = "~")
              if(!tp%in%names(coxInteractionFit)){
                stop("NOTE: 'dualreg' should be listed in 'coxInteractionFit'!\nsee 'tnsGet' function.\n", 
                     call.=FALSE)
              }
              model <- coxInteractionFit[[tp]]
              model$pAdjustInteraction <- coxInteractionTable[tp,"Adjusted.Pvalue"]
            }
            
            #--- get labs
            if(plotype=="3D"){
              xlab <- paste("Regulon activity of",regs[1])
              ylab <- paste("Regulon activity of",regs[2])
              zlab <- "HR"
            } else {
              xlab <- paste("Regulon activity of",regs[1])
              ylab <- paste("Regulon activity\nof",regs[2],"(dES)")
              zlab <- "Hazard Ratio (HR)"
            }
            
            #-- correct names for a 'formula'
            xregs <- .namesCorrect(regs)
            names(xregs) <- regs
            
            #--- set colorPalette
            if(showdata){
              tns <- tnsStratification(tns, nSections = 1, center = TRUE)
              regulonActivity <- tnsGet(tns, "regulonActivity")
              excludeMid <- tnsGet(tns, what = "para")$excludeMid
              regstatus <- .getSurvplotCols(regulonActivity, regs, excludeMid, colorPalette)
              obdata <- stats::model.frame(model, drop.unused.levels=TRUE)
              if(all(rownames(obdata) %in% names(regstatus$cols))){
                datacols <- regstatus$cols[rownames(obdata)]
              } else {
                stop("...unanticipated error occurred while processing this call!")
              }
              if(is.null(xlim)) 
                xlim <- range(obdata[,xregs[1]])*1.04
              xlim <- max(abs(xlim))
              xlim <- c(-xlim,xlim)
              if(is.null(ylim))
                ylim <- range(obdata[,xregs[2]])*1.04
              ylim <- max(abs(ylim))
              ylim <- c(-ylim,ylim)
            } else {
              datacols <- NA
            }
            
            #----- Coxplot
            if(plotpdf){
              fname <- gsub(".pdf", '',fname, ignore.case = TRUE)
              fname <- paste(fname,"_",paste(regs,collapse = "_"),sep = "")
            }
            mbrPlotInteraction(model=model, vars=xregs, xlim=xlim, ylim=ylim, zlim=hlim, 
                               xlab=xlab, ylab=ylab, zlab=zlab, zlog=TRUE,
                               zcols=hcols, showdata=showdata, datacols=datacols,
                               fname=fname, fpath=fpath, width=width, height=height,
                               plotpdf=plotpdf, plotype = plotype)
          })

#' Plot 2-tailed GSEA for a sample from a TNS
#'
#' Makes a 2-tailed GSEA plot for a certain phenotype (sample)
#' present in a TNS. A wrapper of \code{\link{tna.plot.gsea2}}
#'
#' @param tns A \linkS4class{TNS} object
#' @param aSample A string specifying a given sample number present in the 
#' 'survivalData' table.
#' @param regs An optional string vector specifying regulons to make the plot.
#' @param checklog A logical value. If TRUE, expression values are transformed 
#' into log space.
#' @param verbose A logical value specifying to display detailed messages 
#' (when verbose=TRUE) or not (when verbose=FALSE).
#' @param ntop An optional integer value. The number of regulons for which the 
#' GSEA2 will be plotted.
#' @param pValueCutoff An numeric value. The p-value cutoff for the analysis.
#' @param pAdjustMethod A character. Specifies the adjustment method for the 
#' pvalue.
#' See \code{\link{p.adjust}}
#' @param refsamp A character vector.
#' @param plotpdf A single logical value.
#' @param ... parameters which will be passed to \code{\link{tna.plot.gsea2}},
#' such as ylimPanels, heightPanels, width, height, ylabPanels, xlab...
#' @return A plot containing the 2-tailed GSEA analysis for a phenotype.
#' @examples
#' # load survival data
#' data(survival.data)
#' 
#' # load TNI-object
#' data(stni, package = "RTN")
#'
#' stns <- tni2tnsPreprocess(stni, survivalData = survival.data, 
#'         keycovar = c('Grade','Age'), time = 1, event = 2)
#' stns <- tnsGSEA2(stns, verbose=FALSE)
#' tnsPlotGSEA2(stns, 'MB-5115', regs = 'FOXM1', plotpdf = FALSE)
#'
#' @seealso \code{\link{tna.plot.gsea2}} for all plot parameters
#' @importFrom RTN tni.gsea2 tna.plot.gsea2 tna.gsea2 tni2tna.preprocess tni.get
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom graphics abline axis barplot image layout legend lines mtext 
#' par plot plot.new points segments
#' @importFrom stats formula pchisq quantile
#' @docType methods
#' @rdname tnsPlotGSEA2-methods
#' @aliases tnsPlotGSEA2
#' @export
#' 
setMethod("tnsPlotGSEA2", "TNS",
          function(tns, aSample, regs = NULL, refsamp = NULL, checklog = FALSE, 
                   ntop = NULL, pValueCutoff = 0.05, pAdjustMethod = "BH", 
                   verbose = TRUE, plotpdf = FALSE, ...)
          {
            #-- checks
            .tns.checks(tns, type = "TNSpreprocess")
            .tns.checks(aSample, object2 = tns@survivalData, type = "aSample")
            .tns.checks(regs, type = "regs")
            .tns.checks(refsamp, object2 = tns@survivalData, type = "refsamp")
            .tns.checks(checklog, type = "checklog")
            .tns.checks(ntop, type = "ntop")
            .tns.checks(pValueCutoff, type = "pValueCutoff")
            .tns.checks(pAdjustMethod, type = "pAdjustMethod")
            .tns.checks(verbose, type = "verbose")
            .tns.checks(plotpdf, type = "plotpdf")
            
            if (verbose) 
              cat("-Transforming TNS object to inherit methods... \n")
            
            tni <- tnsGet(tns, what = "TNI")
            gexp <- tni.get(tni, what="gexp")
            regulatoryElements <- tni.get(tni, what="regulatoryElements")
            if (!is.null(regs)){
              nin <- length(regs)
              idx1 <- which(regulatoryElements%in%regs)
              idx2 <- which(names(regulatoryElements)%in%regs)
              if(length(idx1)>length(idx2)){
                regs <- regulatoryElements[idx1]
              } else {
                regs <- regulatoryElements[idx2]
              }
              if (length(regs)<nin) 
                stop("NOTE: one or more 'regs' not listedin the TNS object")
            } else {
              regs <- regulatoryElements
            }
            
            ##------ compute reference gx vec
            if (is.null(refsamp)){
              gxref <- apply(gexp, 1, mean)
            } else {
              idx <- colnames(gexp) %in% refsamp
              if (!all(sum(idx) %in% length(refsamp))){
                stop("'refsamp' should list only valid sample names!")
              }
              gxref <- apply(gexp[, refsamp], 1, mean)
            }
            
            ##-- check log space/transform into log space
            qx <- as.numeric(quantile(gexp, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
            LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) || 
              (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
            if (checklog || LogC){
              dt <- log2(1 + gexp) - log2(1 + gxref)
              if (LogC) 
                warning("'gexp' values seem not to be in log space! ..log2 transformation has been applied!")
            } else {
              dt <- gexp - gxref
            }
            
            #-- get phenotype vector
            pheno <- dt[, aSample]
            names(pheno) <- rownames(dt)
            
            #-- make tna from scratch
            rtna <- tni2tna.preprocess(tni, phenotype = pheno, verbose = verbose)
            rtna <- tna.gsea2(rtna, stepFilter = FALSE, pValueCutoff, 
                              pAdjustMethod, tfs = regs, 
                              verbose = verbose)
            
            #-- plot
            tna.plot.gsea2(rtna, labPheno = aSample, tfs = regs, plotpdf = plotpdf, ...=...)
            if(verbose & plotpdf){
              tp1 <- c("NOTE: 'PDF' file for '",aSample,"' should be available either in the\n")
              tp2 <- c("working directory or in a user's custom directory!\n")
              message(tp1,tp2)
            }
            
          })


