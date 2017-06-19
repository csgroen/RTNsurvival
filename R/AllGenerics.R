
setGeneric("tnsPreprocess", function(tni, survivalData, keycovar, time = 1, event = 2, 
    samples = NULL) standardGeneric("tnsPreprocess"), package = "RTNsurvival")

setGeneric("tnsGSEA2", function(tns, ...) standardGeneric("tnsGSEA2"), package = "RTNsurvival")

setGeneric("tnsKM", function(tns, regs = NULL, attribs = NULL, nSections = 2, endpoint = 60, 
    fname = "survplot", fpath = ".", ylab = "Survival probability", xlab = "Months", 
    pal = "redblue", excludeMid = FALSE, flipcols = FALSE, plotpdf = TRUE, plotbatch = FALSE, 
    width = 6.3, height = 3.6, panelWidths = c(3, 2, 4), dES.ylab = "Samples", show.KMlegend = TRUE, 
    KMlegend.pos = "bottomleft", KMlegend.cex = 1, show.pval = TRUE, pval.cex = 1, 
    pval.pos = "topright") standardGeneric("tnsKM"), package = "RTNsurvival")

setGeneric("tnsCox", function(tns, regs = NULL, endpoint = 60, fname = "coxplot", 
    fpath = ".", ylab = "Regulons and key covariates", xlab = "Hazard Ratio (95% CI)", 
    qqkeycovar = FALSE, excludeMid = FALSE, width = 5, height = 5, xlim = c(0.2, 
        10), sortregs = TRUE, plotpdf = TRUE) standardGeneric("tnsCox"), package = "RTNsurvival")

setGeneric("tnsGet", function(object, what) standardGeneric("tnsGet"), 
           package = "RTNsurvival")