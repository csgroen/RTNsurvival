#' Cox regression analysis for dual regulons
#'
#' Returns a data.frame with Cox regression results for the interaction between 
#' regulons identified as duals. 
#'
#' @param mbr an \linkS4class{MBR} object computed from
#' the same \linkS4class{TNI} objects that were used to make the
#' \linkS4class{TNS} objects.
#' @param tns1 a 'TNS' object with regulons used to compute the duals.
#' @param tns2 another 'TNS' object computed with the other regulons used to 
#' compute the duals. It's optional if all duals are from the same 'TNI' object.
#' @param duals an optional vector containing duals to create the table. If NULL,
#' the regression is performed for all duals.
#' @param verbose a logical value. If TRUE, prints status of function while 
#' executing.
#' @param excludeMid if TRUE, doesn't use the samples that fall in mid section 
#' of the stratification for computations.
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
#' results <- mbrGet(mbr, what="dualsInformation")
#' 
#' # create dual Cox regression
#' dualCox <- dualCoxTable(mbr, stns, verbose = FALSE)
#'
#' @seealso \code{\link[RTNduals:tni2mbrPreprocess]{tni2mbrPreprocess}} for all 
#' plot parameters
#' @return A matrix containing the Cox regression results for all given duals
#' @importFrom RTNduals tni2mbrPreprocess mbrAssociation mbrDuals mbrGet
#' @export

dualCoxTable <- function(mbr, tns1, tns2 = NULL, duals = NULL,
                         verbose = TRUE, excludeMid = FALSE) {
    
    #-- checks
    .tns.checks(mbr, type =  "MBR")
    .tns.checks(tns1, type =  "status")
    if (!is.null(tns2)) 
        .tns.checks(tns2, type =  "status")
    .tns.checks(duals, mbr, type = "Duals")
    .tns.checks(verbose, type = "Verbose")
    .tns.checks(excludeMid, type = "ExcludeMid")
    
    #-- get duals
    if (is.null(duals))
        duals <- mbrGet(mbr, what = "dualRegulons")
    
    
    tns1.reg.el <- tnsGet(tns1, "regulatoryElements")
    #-- data wrangling
    if (is.null(tns2)) {
        all.regs <- unique(unlist(strsplit(duals, "~")))
        if(!all(all.regs %in% tns1.reg.el) & !all(all.regs %in% names(tns1.reg.el)))
        {
            stop("If tns2 is not given, all regulons used to compute duals must be
                 present in `regulatoryElements` of tns1.")
        }
        tns2 <- tns1
    }
    else {
        regs1 <- unlist(strsplit(duals, "~"))[c(TRUE, FALSE)]
        regs2 <- unlist(strsplit(duals, "~"))[c(FALSE, TRUE)]
        tns2.reg.el <- tnsGet(tns2, "regulatoryElements")
        if (!all(regs1 %in% tns1.reg.el) & !all(regs1 %in% names(tns1.reg.el))){
            if(!all(regs2 %in% tns1.reg.el) & !all(regs2 %in% names(tns1.reg.el))) {
                stop("`tns1` doesn't contain any useful information.")
            }
            else {
                tmp <- tns2
                tns2 <- tns1
                tns1 <- tmp
                rm(tmp)
                if (!all(regs1 %in% tns2.reg.el) & !all(regs1 %in% names(tns2.reg.el))) {
                    stop("`tns2` doesn't contain any useful information.")
                }    
            }
            
        }
        else {
            if (!all(regs2 %in% tns2.reg.el) & !all(regs2 %in% names(tns2.reg.el))) {
                stop("`tns2` doesn't contain any useful information.") 
            }
        }
    }
    
    #-- gets
    EScores1 <- tnsGet(tns1, what = "EScores")
    EScores2 <- tnsGet(tns2, what = "EScores")
    survData <- tnsGet(tns1, what = "survivalData")
    keycovar <- tnsGet(tns1, what = "keycovar")
    
    #-- mbr vs tns checks
    tns.regs1 <- colnames(EScores1$dif)
    tns.regs2 <- colnames(EScores2$dif)
    
    mbr.regs1 <- unlist(strsplit(duals, "~"))[seq(1, length(duals)*2, 2)]
    mbr.regs2 <- unlist(strsplit(duals, "~"))[seq(2, length(duals)*2, 2)]
    
    if (!all(mbr.regs1 %in% tns.regs1) | !all(mbr.regs2 %in% tns.regs2)) {
        
        idx1 <- which(mbr.regs1 %in% tns.regs1)
        idx2 <- which(mbr.regs2 %in% tns.regs2)
        duals <- duals[intersect(idx1, idx2)]
        
        warning(paste("Not all regulons in the duals have had enrichment scores computed. 
             This is possibly due to minRegulonSize. Regression being computed for only", 
                      length(duals), "duals."), call. = FALSE)
    }
    
    #-- checks
    dif1 <- EScores1$dif
    if (excludeMid) {
        dif1[EScores1$regstatus == EScores1$mid] <- NA
    }
    
    dif2 <- EScores2$dif
    if (excludeMid) {
        dif2[EScores2$regstatus == EScores2$mid] <- NA
    }
    
    dif1 <- dif1[rownames(dif2),]
    
    #-- correct names
    colnames(dif1) <- .namesCorrect(colnames(dif1))
    colnames(dif2) <- .namesCorrect(colnames(dif2))
    
    #--- get cox formula
    if (is.null(tns1@keycovar)) {
        fm1 <- "Surv(time, event) ~ "
    } else {
        fm1 <- paste("Surv(time, event) ~ ", paste(keycovar, collapse = "+"), sep = "")
    }
    
    #-- make summary of data
    summary <- cbind(survData, dif1, dif2)
    
    if (verbose) {
        cat("Calculating cox regression for duals...\n")
        pb <- txtProgressBar(min = 0, max = length(duals), style = 3)
    }
    #--- fit cox regression model
    resall <- sapply(duals, function(dual) {
        regs <- unlist(strsplit(dual, "~"))
        regs <- .namesCorrect(regs)
        nas1 <- is.na(dif1[, regs[1]])
        nas2 <- is.na(dif2[, regs[2]])
        if(verbose)
            setTxtProgressBar(pb, which(duals == dual))
        if (sum(nas1) > nrow(dif1)/2 | sum(nas2) > nrow(dif2)/2) {
            hz <- c(1, 1, 0.99, 1.01,
                    1, 1, 0.99, 1.01,
                    1, 1, 0.99, 1.01)
        }
        
        else {
            coxdual <- paste(regs[1], regs[2], sep = "*")
            fm2 <- formula(paste(fm1, coxdual, sep = "+"))
            if (is.null(keycovar)) {
                p.res <- summary(coxph(fm2, data = summary[!(nas1&nas2), ]))$conf.int[1:3,]
                
            }
            else {
                p.res <- summary(coxph(fm2, data = summary[!(nas1&nas2), ]))$conf.int[(length(keycovar)+1):(length(keycovar)+3),]
                
            }
            p.res <- p.res[,c(1,3:4)]
            res <- as.vector(t(p.res))
        }
    })
    if(verbose)
        close(pb)
    resall <- t(resall)
    colnames(resall) <- c("HR - reg1", "lower .95 - reg1", "upper .95 - reg1", "HR - reg2", 
                    "lower .95 - reg2", "upper .95 - reg2", "HR - dual", "lower .95 - dual", 
                    "upper .95 - dual")
    
    resall <- resall[order(abs(resall[,7]), decreasing = TRUE),]
    return(resall)
}

.namesCorrect <- function(regs) {
    xregs <- gsub("-|\\+|\\.", "_", regs)
    xregs <- gsub("\\s", "", xregs)
    xregs
}

