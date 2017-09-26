
#' Preprocessing of TNS class objects.
#'
#' Creates TNS class onbjects for regulons an survival data.
#'
#' @param tni A \linkS4class{TNI} class, already processed with the same samples
#' listed in the survival data.frame.
#' @param survivalData A named data.frame with samples in rows and survival data 
#' in the columns.
#' @param keycovar A character vector of the 'keycovars' listed in the 
#' data.frame columns.
#' @param time A numeric or character value corresponding to the column of the 
#' data.frame where the time of last observation is given.
#' @param event A numeric or character value, corresponding to the columm of 
#' the data.frame where the 'event' information is given.
#' @param samples An optional character vector listing samples to be analyzed.
#' @return A preprocessed \linkS4class{TNS} class
#' @examples
#' # load survival data
#' data(survival.data)
#' 
#' # load TNI-object
#' data(stni, package = "RTN")
#' 
#' # create a new TNS object
#' stns <- tnsPreprocess(stni, survival.data, keycovar = c('Grade','Age'), 
#' time = 1, event = 2)
#'
#' @seealso \code{\link[RTN:tni.preprocess]{tni.preprocess}} for similar 
#' preprocessing.
#' @import methods
#' @docType methods
#' @rdname tnsPreprocess-methods
#' @aliases tnsPreprocess
#' @export
#' 
setMethod("tnsPreprocess", "TNI", function(tni, survivalData, keycovar = NULL, time = 1, 
    event = 2, samples = NULL) {
    #-- tni checks
    if (tni@status["Preprocess"] != "[x]") 
        stop("NOTE: TNI object requires preprocessing in the 
                       RTN package!")
    if (tni@status["Permutation"] != "[x]") 
        stop("NOTE: TNI object requires permutation/bootstrap and 
                       DPI filter in the RTN package!")
    if (tni@status["DPI.filter"] != "[x]") 
        stop("NOTE: TNI object requires DPI filter in the RTN 
                       package!")
    
    #-- missing
    if (missing(survivalData)) {
        res <- try(tni.get(tni, "colAnnotation"), silent = TRUE)
        if (class(res) == "try-error")
            stop("Must provide a 'survivalData' object.")
        else
            survivalData <- res
    }
    
    #-- par checks
    .tns.checks(survivalData, type = "survivalData")
    time = .tns.checks(time, survivalData, type = "Time")
    event = .tns.checks(event, survivalData, type = "Event")
    .tns.checks(keycovar, survivalData, "Keycovars")
    samples = .tns.checks(samples, survivalData, type = "Samples")
    
    #-- other checks
    if (!all(samples %in% colnames(tni@gexp))) {
        stop("all samples listed in 'survivalData' rownames must be 
available in the 'tni' object!")
    }
    
    #-- reorganize survivalData
    idx <- c(time, event)
    te.data <- survivalData[, idx]
    survivalData <- survivalData[, -idx]
    survivalData <- cbind(te.data, survivalData)
    names(survivalData)[1:2] <- c("time", "event")
    survivalData <- survivalData[samples, ]
    
    #-- making TNS object
    object <- new("TNS", tni = tni, survivalData = survivalData, keycovar = keycovar)
    
    #-- status update
    object <- tns.set(object, what = "status-1")
    
    object
})


#' 2-tailed Gene Set Enrichment Analysis on Transcriptional Networks.
#'
#' Works as a wrapper for \code{\link[RTN:tni.gsea2]{tni.gsea2}}, performing a 
#' 2-tailed GSEA analysis on a \linkS4class{TNI} class object and integrating 
#' the results into the \linkS4class{TNS} class object.
#'
#' @param tns A \linkS4class{TNS} class, which has been preprocessed
#' @param ... Parameters passed to the \code{\link[RTN:tni.gsea2]{tni.gsea2}} 
#' function.
#' @return A \linkS4class{TNS} class, with added Enrichment Scores.
#' @examples
#' # load survival data
#' data(survival.data)
#' 
#' # load TNI-object
#' data(stni, package = "RTN")
#'
#' stns <- tnsPreprocess(stni, survival.data, keycovar = c('Grade','Age'), time = 1, event = 2)
#' stns <- tnsGSEA2(stns, verbose=FALSE)
#'
#' @seealso \code{\link[RTN:tni.gsea2]{tni.gsea2}} for information on all 
#' parameters.
#' @importClassesFrom RTN TNI
#' @docType methods
#' @rdname tnsGSEA2-methods
#' @aliases tnsGSEA2
#' @export
#'
setMethod("tnsGSEA2", "TNS", function(tns, ...) {
    
    #-- checks
    if (tns@status["Preprocess"] != "[x]") 
        stop("NOTE: TNS object requires preprocessing!")
    
    #-- run gsea2 and update TNS
    tni <- tnsGet(tns, what = "TNI")
    EScores <- tni.gsea2(tni, ... = ...)
    tns <- tns.set(tns, EScores, "EScores")
    tns <- tns.set(tns, what = "status-2")
    
    return(tns)
})


#' Kaplan-Meier analysis for TNS class objects.
#'
#' Makes a 2 or 3 panel plot for survival analysis. The first panel shows the
#' differential Enrichment score (dES) for all samples, ranked by expression 
#' in their sections. The second (optional) panel shows the status of other 
#' attributes which may be present in the survival data.frame for all samples. 
#' The third panel shows a Kaplan-Meier plot computed for the given survival 
#' data, with a curve for each section.
#'
#' @param tns a \linkS4class{TNS} object, which must have passed GSEA2 analysis.
#' @param regs a string vector. Contains all the regulons which are going to be 
#' plotted.
#' @param attribs a numeric vector. Contains the columns of the survival 
#' data.frame which will be plotted for the second panel.
#' @param nSections A numeric value for the stratification of the sample. The 
#' larger the number, the more subdivisions will be created for the Kaplan-Meier 
#' analysis.
#' @param endpoint a numeric value. It represents the cut-off point for the 
#' 'time', if any.
#' @param fname a string. The name of the file in which the plot will be saved
#' @param fpath a string. The path to the directory where the plot will be saved
#' @param ylab a string. The label for the y axis on the third panel
#' @param xlab a string. The label for the x axis on the third panel. This should
#' be the measure of time shown in the survival data.frame after the last 
#' check-up.
#' @param pal a string, which can be 'red', 'blue' or 'redblue'. Represents the 
#' colors used in the first and third panels. Alternatively, it can 
#' contains the hex values.
#' @param excludeMid a logical value. If TRUE, inconclusive dES values will not
#'  be consired in the survival analysis.
#' @param flipcols a logical value. If TRUE, flips the order of the samples to 
#' lowest expression on top, highest on the bottom.
#' @param plotpdf a logical value. If TRUE, the plot is saved as a pdf file. 
#' If false, it is plotted in the plotting area.
#' @param plotbatch a logical value. If TRUE, plots for all regs are saved in 
#' the same file.
#' If FALSE, each plot for each reg is saved in a different file.
#' @param width a numeric value. Represents the width of the plot.
#' @param height a numeric value. Represents the height of the plot.
#' @param panelWidths a numeric vector of length=3 specifying the relative 
#' width of the internal panels.
#' @param dES.ylab a string. The label for the y axis of the first panel.
#' @param show.KMlegend a logical value. If TRUE, shows the sample stratification 
#' information on the third panel.
#' @param KMlegend.pos a string. Provides the location of the sample 
#' stratification legend on the third panel. One of: 'bottomright', 
#' 'bottom', 'bottomleft', 'left', 
#' 'topleft', 'top', 'topright', 'right' and 'center'.
#' @param KMlegend.cex a numeric value. Provides the character expansion factor 
#' for the stratification legend, which alters the size and spacing of the font 
#' in the legend.
#' @param show.pval a logical value. If TRUE, shows the Logrank P-value for the 
#' Kaplan-Meier plot on the third panel.
#' @param pval.cex a numeric value. Provides the character expansion factor for 
#' the pvalue legend, which alters the size and spacing of the font in the 
#' legend.
#' @param pval.pos a string. Provides the location of the Logrank P-value on 
#' the third panel. One of: 'bottomright', 'bottom', 'bottomleft', 'left', 
#' 'topleft', 'top', 'topright', 'right' and 'center'.
#' 
#' @return A plot, showing the graphical analysis of provided survival data.
#' @examples
#' # load survival data
#' data(survival.data)
#' 
#' # load TNI-object
#' data(stni, package = "RTN")
#'
#' stns <- tnsPreprocess(stni, survival.data, keycovar = c('Grade','Age'), 
#' time = 1, event = 2)
#' stns <- tnsGSEA2(stns, verbose=FALSE)
#' tnsKM(stns, regs='FOXM1', attribs = list(c('ER+','ER-'),c('G1','G2','G3')), 
#' plotpdf = FALSE)
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom survival survdiff survfit coxph Surv
#' @docType methods
#' @rdname tnsKM-methods
#' @aliases tnsKM
#' @export
#' 
setMethod("tnsKM", "TNS", function(tns, regs = NULL, attribs = NULL, nSections = 2, 
    endpoint = 60, fname = "survplot", fpath = ".", ylab = "Survival probability", 
    xlab = "Months", pal = "redblue", excludeMid = FALSE, flipcols = FALSE, plotpdf = TRUE, 
    plotbatch = FALSE, width = 6.3, height = 3.6, panelWidths = c(3, 2, 4), dES.ylab = "Samples", 
    show.KMlegend = TRUE, KMlegend.pos = "bottomleft", KMlegend.cex = 1, show.pval = TRUE, 
    pval.cex = 1, pval.pos = "topright") {
    #-- checks
    .tns.checks(tns, type = "status")
    .tns.checks(nSections, type = "nSec")
    .tns.checks(fname, type = "Fname")
    .tns.checks(fpath, type = "Path")
    .tns.checks(ylab, type = "Ylab")
    .tns.checks(xlab, type = "Xlab")
    .tns.checks(regs, type = "Regs")
    .tns.checks(attribs, tns@survivalData, type = "Attribs")
    .tns.checks(pal, tns@para$strat, type = "Pal")
    .tns.checks(excludeMid, type = "ExcludeMid")
    .tns.checks(flipcols, type = "FlipCols")
    .tns.checks(plotpdf, type = "PlotPDF")
    .tns.checks(plotbatch, type = "PlotBatch")
    .tns.checks(width, height, type = "WidthHeight")
    .tns.checks(endpoint, type = "EndPoint")
    .tns.checks(panelWidths, type = "panelWidths")
    .tns.checks(dES.ylab, type = "dES.ylab")
    .tns.checks(show.KMlegend, type = "showKMlegend")
    .tns.checks(KMlegend.pos, type = "KMlegend.pos")
    .tns.checks(KMlegend.cex, type = "KMlegend.cex")
    .tns.checks(show.pval, type = "show.pval")
    .tns.checks(pval.cex, type = "pval.cex")
    .tns.checks(pval.pos, type = "pval.pos")
    
    #-- stratification
    tns <- .tns.stratification(tns, nSections = nSections)
    
    #-- organize data
    survData <- tnsGet(tns, what = "survivalData")
    EScores <- tnsGet(tns, what = "EScores")
    #-- 
    survData$event[survData$time > endpoint] <- 0
    survData$time[survData$time > endpoint] <- endpoint
    #-- 
    sp1 <- rownames(survData)
    sp2 <- rownames(EScores$regstatus)
    survData <- survData[sp1 %in% sp2, , drop = FALSE]
    EScores$pos <- EScores$pos[sp2 %in% sp1, , drop = FALSE]
    EScores$neg <- EScores$neg[sp2 %in% sp1, , drop = FALSE]
    EScores$dif <- EScores$dif[sp2 %in% sp1, , drop = FALSE]
    #-- 
    tns <- tns.set(tns, EScores, what = "EScores")
    tns <- tns.set(tns, survData, what = "survivalData")
    
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
            stop("NOTE: 'attribs' variables should only include binary values!")
    }
    
    #-- making reglist
    reglist <- colnames(tns@EScores$regstatus)
    if (!is.null(regs)) {
        if (!all(regs %in% reglist)) {
            stop("NOTE: all names in 'regs' should be listed 
                     in the slot 'EScores' of the 'tns' object!")
        }
        reglist <- regs
    }
    
    idx <- apply(EScores$regstatus, 2, function(es) {
        any(is.na(es))
    })
    validregs <- colnames(EScores$regstatus)[!idx]
    reglist <- reglist[reglist %in% validregs]
    
    #---plot
    if (plotbatch & plotpdf) {
        pdf(file = paste(fpath, "/", fname, ".pdf", sep = ""), width = width, height = height)
        for (reg in reglist) {
            .survplot(EScores, survData, reg, fname, fpath, ylab, xlab, pal, panelWidths, 
                plotpdf, excludeMid, flipcols, attribs, groups, endpoint, dES.ylab = dES.ylab, 
                show.KMlegend = show.KMlegend, KMlegend.pos = KMlegend.pos, KMlegend.cex = KMlegend.cex, 
                show.pval = show.pval, pval.cex = pval.cex, pval.pos = pval.pos)
        }
        dev.off()
    } else {
        for (reg in reglist) {
            if (plotpdf) {
                pdf(file = paste(fpath, "/", reg, fname, ".pdf", sep = ""), width = width, 
                  height = height)
            }
            .survplot(EScores, survData, reg, fname, fpath, ylab, xlab, pal, panelWidths, 
                plotpdf, excludeMid, flipcols, attribs, groups, endpoint, dES.ylab = dES.ylab, 
                show.KMlegend = show.KMlegend, KMlegend.pos = KMlegend.pos, KMlegend.cex = KMlegend.cex, 
                show.pval = show.pval, pval.cex = pval.cex, pval.pos = pval.pos)
            if (plotpdf) {
                message("NOTE: 'PDF' file was generated")
                dev.off()
            }
        }
    }
    
    invisible(list(EScores = EScores, survivalData = survData))
})


#' Cox regression analysis for TNS class objects.
#'
#' Run Cox multivariate regression for regulons and key covariables.
#'
#' @param tns a \linkS4class{TNS} object, which must have passed GSEA2 analysis.
#' @param regs a string vector. Contains the regulons which will be used to 
#' compute the Cox multivariate model. If left NULL, all regulons will be used.
#' @param endpoint a numeric value. The final point in time for the samples. All
#' time values larger than endpoint will be set at endpoint.
#' @param fname a string. The name of the PDF file which will contain the plot.
#' @param fpath a string. The directory where the file will be saved.
#' @param ylab a string. The label of the y-axis, describing what is represented.
#' @param xlab a string. The label of the x-axis.
#' @param qqkeycovar a logical value. If TRUE, only the samples in the 2nd and 
#' 3rd quarters of dES are used to compute. If FALSE, all samples are used.
#' @param excludeMid a logical value. If TRUE, inconclusive dES values will not be
#' consired in the survival analysis.
#' @param width a numeric value. The width of the plot.
#' @param height a numeric value. The height of the plot.
#' @param xlim a vector with 2 values. The first value represents the lowest
#' value in the x-axis, the second value is the highest.
#' @param sortregs a logical value. If TRUE, regulons are sorted from most 
#' negatively associated with hazard to most positively associated with hazard.
#' @param plotpdf a logical value.
#' @return A Cox hazard model plot. If TRUE, generates a pdf plot.
#' @examples
#' # load survival data
#' data(survival.data)
#' 
#' # load TNI-object
#' data(stni, package = "RTN")
#'
#' stns <- tnsPreprocess(stni, survival.data, keycovar = c('Grade','Age'), 
#' time = 1, event = 2)
#' stns <- tnsGSEA2(stns, verbose=FALSE)
#' tnsCox(stns, regs = c('PTTG1','E2F2','FOXM1'), sortregs = TRUE, 
#' plotpdf = FALSE)
#' 
#' @docType methods
#' @rdname tnsCox-methods
#' @aliases tnsCox
#' @export
#' 
setMethod("tnsCox", "TNS", function(tns, regs = NULL, endpoint = 60, fname = "coxplot", 
    fpath = ".", ylab = "Regulons and key covariates", xlab = "Hazard Ratio (95% CI)", 
    qqkeycovar = FALSE, excludeMid = FALSE, width = 5, height = 5, xlim = c(0.2, 
        10), sortregs = TRUE, plotpdf = TRUE) {
    
    #-- checks
    .tns.checks(tns, type = "status")
    .tns.checks(regs, type = "Regs")
    .tns.checks(fname, type = "Fname")
    .tns.checks(fpath, type = "Path")
    .tns.checks(ylab, type = "Ylab")
    .tns.checks(xlab, type = "Xlab")
    .tns.checks(qqkeycovar, type = "QQCovar")
    .tns.checks(endpoint, type = "EndPoint")
    .tns.checks(excludeMid, type = "ExcludeMid")
    .tns.checks(width, height, type = "WidthHeight")
    .tns.checks(xlim, type = "Xlim")
    .tns.checks(sortregs, type = "SortRegs")
    .tns.checks(plotpdf, type = "PlotPDF")
    .tns.checks(tns@survivalData, type = "survival_cox")
    
    #-- gets
    EScores <- tnsGet(tns, what = "EScores")
    survData <- tnsGet(tns, what = "survivalData")
    keycovar <- tnsGet(tns, what = "keycovar")
    
    #-- checks
    dif <- EScores$dif
    if (excludeMid) {
        dif[EScores$regstatus == EScores$mid] <- NA
    }
    
    if (!is.null(regs)) {
        if (!all(regs %in% colnames(dif))) {
            stop("Not all 'regs' have EScores!")
        }
        idx <- colnames(dif) %in% regs
        dif <- dif[, idx]
        dif <- dif[, regs]
    }
    
    #---set names to a valid format
    regs <- colnames(dif)
    xregs <- gsub("-|\\+|\\.", "_", regs)
    xregs <- gsub("\\s", "", xregs)
    colnames(dif) <- xregs
    names(regs) <- xregs
    
    #---combine dif and survivalData 
    summary <- cbind(survData[rownames(dif), ], dif)
    
    #--- set keycovar by quantile
    kvarlist <- list()
    if (qqkeycovar) {
        for (kvar in tns@keycovar) {
            tp <- summary[[kvar]]
            ql <- quantile(tp, c(0.25, 0.75), na.rm = TRUE)
            tp[tp < ql[1]] <- NA
            tp[tp > ql[2]] <- NA
            summary[[kvar]] <- tp
            kvarlist[[kvar]] <- ql
        }
    }
    
    #--- filter data
    summary <- summary[, c("time", "event", keycovar, xregs)]
    
    #--- get cox formula
    if (is.null(tns@keycovar)) {
        fm1 <- "Surv(time, event) ~ "
    } else {
        fm1 <- paste("Surv(time, event) ~ ", paste(keycovar, collapse = "+"), sep = "")
    }
    
    #--- fit cox regression model
    resall <- sapply(xregs, function(rg) {
        nas <- is.na(summary[, rg])
        if (sum(nas) > nrow(summary)/2) {
            c(1, 1, 0.99, 1.01)
        } else {
            fm2 <- formula(paste(fm1, rg, sep = "+"))
            summary(coxph(fm2, data = summary[!nas, ]))$conf.int[rg, , drop = FALSE]
        }
    })
    resall <- t(resall)
    dimnames(resall) = list(xregs, 
                    c("exp(coef)", "exp(-coef)", "lower .95", "upper .95"))
    if (sortregs) {
        resall <- resall[sort.list(resall[, 1]), ]
    }
    
    #--- fit cox model for keycovars and adding to resall
    idx <- which.max(resall[, "exp(coef)"])
    fm2 <- formula(paste(fm1, rownames(resall)[idx], sep = "+"))
    resref <- summary(coxph(fm2, data = summary))$conf.int
    resall <- rbind(resref[-nrow(resref), ], resall)
    rownames(resall)[1:length(keycovar)] <- tnsGet(tns, "keycovar")
    resall <- resall[nrow(resall):1, ]
    
    #--- add symbols to rownames
    idx <- match(names(regs), rownames(resall))
    rownames(resall)[idx] <- regs
    
    #--- plot
    filen = paste(fpath, "/", fname, ".pdf", sep = "")
    .plotCox(resall, regs = regs, keycovar = keycovar, filen = filen, width = width, 
        height = height, xlim = xlim, xlab = xlab, ylab = ylab, plotpdf = plotpdf)
    
    #--- return
    invisible(list(resall = resall, kvarlist = kvarlist))
    
})

setMethod("show", "TNS", function(object) {
    message("a TNS (Transcriptional Network - Survival) object:\n")
    message("--status:")
    print(object@status, quote = FALSE)
})


#' Get information from slots in a TNS object
#'
#'Get information from individual slots in a TNS object and 
#'any available results from a previous analysis.
#'
#' @param object a TNS object
#' @param what a character vector specifying what should be retrieved from the
#' object. Options: 'survivalData', 'EScores', 'TNI', 'keycovar'
#' 
#' @return A plot containing the 2-tailed GSEA analysis for a phenotype.
#' @examples
#' # load survival data
#' data(survival.data)
#' 
#' # load TNI-object
#' data(stni, package = "RTN")
#'
#' stns <- tnsPreprocess(stni, survival.data, keycovar = c('Grade','Age'), 
#' time = 1, event = 2)
#' stns <- tnsGSEA2(stns, verbose=FALSE)
#' enrichmentScores <- tnsGet(stns, 'EScores')
#'
#' @docType methods
#' @rdname tnsGet-methods
#' @aliases tnsGet
#' @export

setMethod("tnsGet", "TNS", function(object, what)
{
    if (what == "survivalData") {
        return(object@survivalData) 
    } 
    else if (what == "EScores") {
        return(object@EScores) 
    } 
    else if (what == "TNI") {
        return(object@tni) 
    } 
    else if (what == "keycovar") {
        return(object@keycovar)
    }
    else if (what == "regulatoryElements") {
        return(object@tni@regulatoryElements)
    }
    else {
        stop("'what' must be one of: 'survivalData', 'EScores', 'TNI',
             'keycovar' and 'regulatoryElements'.")
    }
})

#' Panel plot of regulons in a inferred dual
#'
#' Returns a folder with 5 plots that form a panel: 
#' [1] sample stratification based on activity of regulon 1, 
#' [2] sample stratification based on activity of regulon 2, 
#' [3] a scatter plot of samples's rank in both stratifications based on 
#' regulon activity
#' [4] Kaplan-Meier curves based on both regulon stratifications and interaction 
#' of regulons in dual,
#' [5] Cox regression for both regulons and the interaction between the activity
#' information for the dual. 
#' The function also returns a PNG panel with all plots aligned.
#'
#' @param mbr an \linkS4class{MBR} object computed from
#' the same \linkS4class{TNI} objects that were used to make the
#' \linkS4class{TNS} objects.
#' @param tns1 a `TNS` object with regulons used to compute the duals.
#' @param tns2 another `TNS`  object computed with the other regulons used to 
#' compute the duals. It's optional if all duals are from the same 'TNI' object.
#' @param dual a character string with the name of a dual
#' @param attribs a numeric vector. Contains the columns of the survival 
#' data.frame which will be plotted as sample status in first and second plots.
#' @param endpoint a numeric value. It represents the cut-off point for the 
#' 'time', if any.
#' @param path a string. The path to the directory where the plot will be saved.
#' @param nSections numeric value for the stratification of the sample. The 
#' larger the number, the more subdivisions will be created for the 
#' Kaplan-Meier and sample stratification analysis.
#' @param pal a palette name to be passed to \code{\link[RColorBrewer:brewer.pal]{brewer.pal}}.
#' @param excludeMid a logical value. if TRUE, doesn't use the samples that fall in mid section 
#' of the stratification for computations.
#' @param sectionsLegend a string vector. It must contain names for the sections
#' of the stratification on the interaction Kaplan-Meier. If NULL, the sections
#' are simply numbered in the legend.
#' @param panelWidths a numeric vector of length 2. It must contain a proportion
#' for the two panels in plots [1] and [2].
#' @param png.res a numeric value. It represents the resolution in ppi for the
#' png panel.
#' 
#' @examples
#' # load survival data
#' data(survival.data)
#' 
#' # load TNI-object
#' data(stni, package = "RTN")
#' stni <- upgradeTNI(stni)
#' 
#' # perform survival analysis
#' stns <- tnsPreprocess(stni, survival.data, keycovar = c('Grade','Age'),
#'                       time = 1, event = 2)
#' stns <- tnsGSEA2(stns, verbose=FALSE)
#' 
#' # create MBR-object using TF-TF duals
#' library(RTNduals)
#' tni <- tnsGet(stns, "TNI")
#' mbr <- tni2mbrPreprocess(tni, tni, verbose = FALSE)
#' mbr <- mbrAssociation(mbr, prob = 0.75)
#' mbr <- mbrDuals(mbr)
#' duals <- mbrGet(mbr, what="dualRegulons")
#' 
#' # create panel plot
#' 
#' dualSurvivalPanel(mbr, stns, dual = duals[1], 
#'                   attribs = c("ER+", "ER-", "PR+", "PR-"))
#'
#' @seealso \code{\link[RTNsurvival:tnsKM]{tnsKM}},
#' \code{\link[RTNsurvival:tnsCox]{tnsCox}}
#' @return A folder containing 5 pdf survival plots and one png survival panel
#' with all plots
#' @importFrom RTN tni.get
#' @importClassesFrom RTNduals MBR
#' @importFrom png readPNG
#' @importFrom graphics rasterImage box
#' @importFrom grDevices png
#' @importFrom stats lm cor
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @docType methods
#' @rdname dualSurvivalPanel-methods
#' @aliases dualSurvivalPanel
#' @export
#' 
setMethod("dualSurvivalPanel", "MBR", function(mbr, tns1, tns2 = NULL, dual, attribs = NULL, endpoint = 60,
                   path = NULL,  nSections = 2, 
                   pal = "BrBG", excludeMid = FALSE, sectionsLegend = NULL,
                   panelWidths = c(2,3), png.res = 150)
{
    
    .tns.checks(mbr, type = "MBR")
    .tns.checks(tns1, type = "status")
    if (!is.null(tns2))
        .tns.checks(tns2, type = "status")
    .tns.checks(dual, mbr, type = "Duals")
    .tns.checks(attribs, tnsGet(tns1, "survivalData"), type = "attribs")
    .tns.checks(endpoint, type = "EndPoint")
    .tns.checks(path, type = "Path2")
    .tns.checks(pal, type = "CBpal")
    .tns.checks(nSections, type = "NSections")
    .tns.checks(panelWidths, type = "PanelWidths2")
    .tns.checks(png.res, type = "png.res")

    #-- check sections legend
    if(!is.null(sectionsLegend)) {
        if(excludeMid) {
            if (length(sectionsLegend) != nSections*2) {
                stop("`sectionsLegend` length doesn't match number of sections.",
                     .call = FALSE)
            }
        }
        else {
            if (length(sectionsLegend) != nSections*2+1) {
                stop("`sectionsLegend` length doesn't match number of sections.",
                     .call = FALSE)
            }
            }
    }
    regs <- unlist(strsplit(dual, "~"))
    
    #--- creating folder (if not existent)
    if (is.null(path)) {
        if(!dir.exists(dual))
            dir.create(dual)
        path <- dual
    }
    else if(!dir.exists(path)) {
        dir.create(path)
    }
    
    #-- data wrangling
    tns1.reg.el <- tnsGet(tns1, "regulatoryElements")
    if (is.null(tns2)) {
        all.regs <- unlist(strsplit(dual, "~"))
        if(!all(all.regs %in% tns1.reg.el) & !all(all.regs %in% names(tns1.reg.el)))
        {
            stop("If tns2 is not given, all regulons used to compute duals must be
                 present in `regulatoryElements` of tns1.")
        }
        tns2 <- tns1
    }
    else {
        tns2.reg.el <- tnsGet(tns2, "regulatoryElements")
        if (!(regs[1] %in% tns1.reg.el) & !(regs[1] %in% names(tns1.reg.el))){
            if(!(regs[2] %in% tns1.reg.el) & !(regs[2] %in% names(tns1.reg.el))) {
                stop("`tns1` doesn't contain any useful information.")
            }
            else {
                tmp <- tns2
                tns2 <- tns1
                tns1 <- tmp
                rm(tmp)
                if (!(regs[1] %in% tns2.reg.el) & !(regs[1] %in% names(tns2.reg.el))) {
                    stop("`tns2` doesn't contain any useful information.")
                }    
            }
            
        }
        else {
            if (!(regs[2] %in% tns2.reg.el) & !(regs[2] %in% names(tns2.reg.el))) {
                stop("`tns2` doesn't contain any useful information.") 
            }
        }
    }
    #-- check if tns are for same gexp
    if(!identical(tni.get(tnsGet(tns1, "TNI"), "gexp"), 
              tni.get(tnsGet(tns2, "TNI"), "gexp"))) {
        stop("The expression matrix (gexp) used to compute the TNIs must be the same.")
    }
    
    #--- Coxtable
    dualCoxTable <- dualCoxTable(mbr, tns1, tns2, duals = dual, 
                                 excludeMid = excludeMid, verbose = FALSE)
    
    #----- dES plots
    fname <- paste(path, "/", "1.dESregPlot_", regs[1], ".pdf", sep = "")
    pdf(file = fname, width = 3.897, height = 4.1335)
    dESregPlot(tns1, regs = regs[1], attribs = attribs, 
               panelWidths = panelWidths, flipGraph = TRUE, nSections = nSections)
    dev.off()
    
    fname <- paste(path, "/", ".1.dESregPlot_", regs[1], ".png", sep = "")
    png(filename = fname, width = 3.897*png.res, height = 4.1335*png.res, res = png.res)
    dESregPlot(tns1, regs = regs[1], attribs = attribs, 
               panelWidths = panelWidths, flipGraph = TRUE, nSections = nSections,
               xname = "a")
    dev.off()
    
    fname <- paste(path, "/", "2.dESregPlot_", regs[2], ".pdf", sep = "")
    pdf(file = fname, width = 3.897, height = 4.1335)
    dESregPlot(tns2, regs = regs[2], attribs = attribs, 
               panelWidths = panelWidths, flipGraph = FALSE, nSections = nSections)
    dev.off()
    
    fname <- paste(path, "/", ".2.dESregPlot_", regs[2], ".png", sep = "")
    png(filename = fname, width = 3.897*png.res, height = 4.1335*png.res, res = png.res)
    dESregPlot(tns2, regs = regs[2], attribs = attribs, 
               panelWidths = panelWidths, flipGraph = FALSE, nSections = nSections,
               xname = "b")
    dev.off()
    
    #----- KM plots (separate regs)
    fname <- paste(path, "/", "4.KMplot_", dual, ".pdf", sep = "")
    pdf(file = fname, width = 3.897, height = 6.3)
    layout(matrix(c(1,2,3,3), ncol = 2, byrow = TRUE), widths = c(4,3.3), heights = c(3,5))
    par(mar = c(1,1,1,1))
    dES1 <- KMregPlot(tns1, reg = regs[1], endpoint = endpoint, xlab = "", show.KMlegend = FALSE, pval.cex = 0.9,
                      pval.pos = "bottomright", ylab.cex = 1.2, xlab.cex = 1.2, nSections = nSections,
                      excludeMid = excludeMid)
    dES2 <- KMregPlot(tns2, reg = regs[2], endpoint = endpoint, xlab = "", show.KMlegend = FALSE, pval.cex = 0.9,
                      pval.pos = "bottomright", y.axis = TRUE, ylab.cex = 1.2, xlab.cex = 1.2, ylab = "", 
                      nSections = nSections, excludeMid = excludeMid)
    dES2$regstatus <- dES2$regstatus[rownames(dES1$regstatus),]
    samples <- rownames(dES1$regstatus)[dES1$regstatus[,regs[1]] == dES2$regstatus[,regs[2]]]
    pal.cols <- brewer.pal((nSections*4), pal)[c(2:(nSections+1), 
                                                 (nSections+3):(nSections*4))]
    if (!excludeMid)
        pal.cols <- c(pal.cols[1:nSections], "grey80", pal.cols[nSections+2:length(pal.cols)])
    #-- if samples are disagreeing: use opposites for KM
    samples_cor <- cor(dES1$regstatus[,regs[1]], dES2$regstatus[,regs[2]])
    
    if (samples_cor < 0) {
        samples <- NULL
        nStrat <- nSections*2+1
        up <- 1:nStrat
        down <- nStrat:1
        up <- up[-(nSections+1)]
        down <- down[-(nSections+1)]
        
        for (i in 1:length(up)) {
            current_samps <- rownames(dES1$regstatus)[dES1$regstatus[,regs[1]] == up[i]]
            current_samps[current_samps %in% rownames(dES2$regstatus)[dES2$regstatus[,regs[2]] == down[i]]]
            samples <- c(samples, current_samps)
        }
        
        KMregPlot(tns1, reg = regs[1], endpoint = endpoint, xlab = "Months", show.KMlegend = TRUE, pval.cex = 1,
                  pval.pos = "bottomright", title = "Interaction", samples = samples, KMlegend.cex = 0.7, 
                  ylab.cex = 1.2, xlab.cex = 1.2, title.cex = 1.2, nSections = nSections,
                  excludeMid = excludeMid, pal = pal.cols)
        mode <- "disagreement"
        
    }
    else {
        KMregPlot(tns1, reg = regs[1], endpoint = endpoint, xlab = "Months", show.KMlegend = TRUE, pval.cex = 1,
                  title = "Interaction", ylab.cex = 1.2, xlab.cex = 1.2, title.cex = 1.2, nSections = nSections, pal = pal.cols,
                  excludeMid = excludeMid)
        mode <- "agreement"
    }
    dev.off()
    
    fname <- paste(path, "/", ".4.KMplot_", dual, ".png", sep = "")
    png(filename = fname, width = 3.897*png.res, height = 6.3*png.res, res = png.res)
    layout(matrix(c(1,2,3,3), ncol = 2, byrow = TRUE), widths = c(4,3.3), heights = c(3,5))
    par(mar = c(1,1,1,1))
    KMregPlot(tns1, reg = regs[1], endpoint = endpoint, xlab = "", show.KMlegend = FALSE, pval.cex = 0.9,
              pval.pos = "bottomright", ylab.cex = 1.2, xlab.cex = 1.2, nSections = nSections,
              excludeMid = excludeMid, xname = "c")
    KMregPlot(tns2, reg = regs[2], endpoint = endpoint, xlab = "", show.KMlegend = FALSE, pval.cex = 0.9,
              pval.pos = "bottomright", y.axis = TRUE, ylab.cex = 1.2, xlab.cex = 1.2, ylab = "", 
              nSections = nSections, excludeMid = excludeMid)
    
    if (mode == "disagreement") {
        KMregPlot(tns1, reg = regs[1], endpoint = endpoint, xlab = "Months", show.KMlegend = TRUE, pval.cex = 1,
                  pval.pos = "bottomright", title = "Interaction", samples = samples, KMlegend.cex = 0.7, 
                  ylab.cex = 1.2, xlab.cex = 1.2, title.cex = 1.2, nSections = nSections,
                  excludeMid = excludeMid, pal = pal.cols)
    }
    else {
        KMregPlot(tns1, reg = regs[1], endpoint = endpoint, xlab = "Months", show.KMlegend = TRUE, pval.cex = 1,
                  pval.pos = "bottomright", title = "Interaction", samples = samples, KMlegend.cex = 0.8, 
                  ylab.cex = 1.2, xlab.cex = 1.2, title.cex = 1.2, nSections = nSections, pal = pal.cols,
                  excludeMid = excludeMid)
    }
    
    dev.off()
    
    #----- Scatter plot
    
    fname <- paste(path, "/",  "3.rankScatter_", dual, ".pdf", sep = "")
    pdf(fname, width = 3.3, height = 3.5)
    rankScatter(tns1, tns2, dual, dES1, dES2, nSections = nSections, mode = mode,
                pal = pal)
    dev.off()
    
    fname <- paste(path, "/",  ".3.rankScatter_", dual, ".png", sep = "")
    png(fname, width = 3.3*png.res, height = 3.5*png.res, res = png.res)
    rankScatter(tns1, tns2, dual, dES1, dES2, nSections = nSections, mode = mode,
                pal = pal)
    dev.off()
    
    
    #----- Cox plot
    
    fname <- paste(path, "/",  "5.coxplot_", dual, ".pdf", sep = "")
    pdf(file = fname, width = 4.25, height = 2.5)
    dualCoxPlot(dual, dualCoxTable)
    dev.off()
    
    fname <- paste(path, "/",  ".5.coxplot_", dual, ".png", sep = "")
    png(filename = fname, width = 4.25*png.res, height = 2.5*png.res, res = png.res)
    dualCoxPlot(dual, dualCoxTable, xname = "d")
    dev.off()
    
    pngPanel(dual, path, png.res)
    
    msg <- paste0("Plots successfull. ", "Files can be found in ", path, " directory. ", 
                  "More information on the interpretation of the plot can be found in Figure 5 in vignette('RTNsurvival')")
    print(msg)
    
})

