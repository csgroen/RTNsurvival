
setGeneric("tni2tnsPreprocess", 
           function(tni, survivalData = NULL, regulatoryElements = NULL, 
                    time = 1, event = 2, endpoint = NULL, pAdjustMethod = "BH", 
                    keycovar = NULL, samples = NULL, excludeMid = FALSE) 
             standardGeneric("tni2tnsPreprocess"), package = "RTNsurvival")

setGeneric("tnsGSEA2", 
           function(tns, ...)
             standardGeneric("tnsGSEA2"), package = "RTNsurvival")

setGeneric("tnsAREA3", 
           function(tns, ...)
             standardGeneric("tnsAREA3"), package = "RTNsurvival")

setGeneric("tnsKM", 
           function(tns, regs = NULL, nSections = 1, verbose = TRUE) 
             standardGeneric("tnsKM"), package = "RTNsurvival")

setGeneric("tnsPlotKM", 
           function(tns, regs = NULL, attribs = NULL, fname = "survplot", 
                    fpath = ".", xlab = "Months", ylab = "Survival probability", 
                    colorPalette = "bluered", plotpdf = FALSE, plotbatch = FALSE, 
                    width = 6.3, height = 3.6, panelWidths = c(3, 2, 4)) 
             standardGeneric("tnsPlotKM"), package = "RTNsurvival")

setGeneric("tnsCox", 
           function(tns, regs = NULL, qqkeycovar = FALSE, verbose = TRUE) 
             standardGeneric("tnsCox"), package = "RTNsurvival")

setGeneric("tnsPlotCox", 
           function(tns, regs = NULL, fname = "coxplot", 
                    fpath = ".", ylab = "Regulons and other covariates", 
                    xlab = "Hazard Ratio (95% CI)", width = 5, 
                    height = 5, xlim = c(0.3, 3), 
                    sortregs = TRUE, plotpdf = FALSE) 
             standardGeneric("tnsPlotCox"), package = "RTNsurvival")

setGeneric("tnsGet", 
           function(tns, what) standardGeneric("tnsGet"), 
           package = "RTNsurvival")

setGeneric("tnsInteraction",
           function(tns, ..., verbose = TRUE)
             standardGeneric("tnsInteraction"), package = "RTNsurvival")

setGeneric("tnsKmInteraction",
           function(tns, mbr, stepFilter = TRUE, pValueCutoff = 0.05, verbose = TRUE)
             standardGeneric("tnsKmInteraction"), package = "RTNsurvival")

setGeneric("tnsPlotKmInteraction", 
           function(tns, dualreg = NULL, fname = "kmInteraction", 
                    fpath = ".", xlab = "Months", 
                    ylab = "Survival probability", colorPalette = "bluered", 
                    width = 4, height = 4, plotpdf = FALSE) 
             standardGeneric("tnsPlotKmInteraction"), package = "RTNsurvival")

setGeneric("tnsCoxInteraction",
           function(tns, mbr, stepFilter = TRUE, pValueCutoff = 0.05, verbose = TRUE) 
             standardGeneric("tnsCoxInteraction"), package = "RTNsurvival")

setGeneric("tnsPlotCoxInteraction",
           function(tns, dualreg, xlim = NULL, ylim = NULL, hlim = NULL, 
                    hcols = c("#008080ff","#d45500ff"), showdata = TRUE, 
                    colorPalette = "bluered", fname = "coxInteraction", 
                    fpath = ".", width = 5, height = 4, plotype = "3D", plotpdf = FALSE) 
             standardGeneric("tnsPlotCoxInteraction"), package = "RTNsurvival")

setGeneric("tnsPlotCoxInteraction",
           function(tns, dualreg, xlim = NULL, ylim = NULL, hlim = NULL, 
                    hcols = c("#008080ff","#d45500ff"), showdata = TRUE, 
                    colorPalette = "bluered", fname = "coxInteraction", 
                    fpath = ".", width = 5, height = 4, plotype = "3D", plotpdf = FALSE) 
               standardGeneric("tnsPlotCoxInteraction"), package = "RTNsurvival")

setGeneric("tnsPlotGSEA2",
           function(tns, aSample, regs = NULL, refsamp = NULL, checklog = FALSE, 
                    ntop = NULL, pValueCutoff = 0.05, pAdjustMethod = "BH", 
                    verbose = TRUE, plotpdf = FALSE, ...) 
             standardGeneric("tnsPlotGSEA2"), package = "RTNsurvival")

setGeneric("tnsPlotCovariates",
           function(tns, regs = NULL, attribs = NULL, fname = "covarplot", 
                    fpath = ".", plotpdf = FALSE, plotbatch = FALSE,
                    panelHeights = c(1,1), width = 5.3, height = 4,
                    dummyEncode = TRUE, divs = NULL)
               standardGeneric("tnsPlotCovariates"), package = "RTNsurvival")

setGeneric("tnsSRE",
           function(tns, subgroup, regs = NULL, pValueCutoff = 0.05, 
                    pAdjustMethod = "BH")
               standardGeneric("tnsSRE"), package = "RTNsurvival")

setGeneric("tnsPlotSRE",
           function(tns, subgroup = NULL, by = "nGroups",
                    nGroupsEnriched = 1, nTopEnriched = 10, 
                    breaks = seq(-1.5, 1.5, 0.1),
                    markEnriched = FALSE, ...)
               standardGeneric("tnsPlotSRE"), package = "RTNsurvival")

