
#' Plot 2-tailed GSEA for a sample from a TNS
#'
#' Makes a 2-tailed GSEA plot for a certain phenotype (sample)
#' present in a TNS. A wrapper of \code{\link[RTN:tna.plot.gsea2]{tna.plot.gsea2}}
#'
#' @param object a TNS object
#' @param aSample a string specifying a given sample number present in the 
#' 'survivalData' table.
#' @param regs an optional string vector specifying regulons to make the plot.
#' @param log a logical value. If TRUE, gexp values are transformed into log 
#' space.
#' @param verbose a logical value. If TRUE, prints status of function while 
#' executing.
#' @param ntop an optional integer value. The number of regulons for which the 
#' GSEA2 will be plotted.
#' @param pValueCutoff an integer. The p cutoff value for the analysis.
#' @param pAdjustMethod a character. Specifies the adjustment method for the 
#' pvalue.
#' See \code{\link[stats:p.adjust]{p.adjust}}
#' @param refsamp a character vector.
#' @param plotpdf a single logical value.
#' @param ... parameters which will be passed to 
#' \code{\link[RTN:tna.plot.gsea2]{tna.plot.gsea2}},
#' such as ylimPanels, heightPanels, width, height, ylabPanels, xlab...
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
#' tnsPlotGSEA2(stns, 'MB-5115', regs = 'FOXM1')
#'
#' @seealso \code{\link[RTN:tna.plot.gsea2]{tna.plot.gsea2}} for all plot parameters
#' @importFrom RTN tni.gsea2 tna.plot.gsea2 tna.gsea2 tni2tna.preprocess
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom graphics abline axis barplot image layout legend lines mtext par 
#' plot plot.new points segments
#' @importFrom stats formula pchisq quantile
#' @export
#' 
tnsPlotGSEA2 <- function(object, aSample, regs = NULL, refsamp = NULL, log = FALSE, 
    ntop = NULL, pValueCutoff = 0.05, pAdjustMethod = "BH", verbose = TRUE, plotpdf = TRUE, 
    ...)
    {
    #-- checks
    if (class(object) != "TNS") 
        stop("NOTE: 'tnsPlotGSEA2' requires a 'TNS' class object!")
    if (object@status["Preprocess"] != "[x]") 
        stop("NOTE: TNS object requires preprocessing!")
    if (object@status["GSEA2"] != "[x]") 
        stop("NOTE: TNS object needs to be evaluated by 'tnsGSEA2'!")
    .tns.checks(aSample, object2 = object@survivalData, type = "aSample")
    .tns.checks(regs, type = "Regs")
    .tns.checks(refsamp, object2 = object@survivalData, type = "Refsamp")
    .tns.checks(log, type = "Log")
    .tns.checks(ntop, type = "Ntop")
    .tns.checks(pValueCutoff, type = "pValueCutoff")
    .tns.checks(pAdjustMethod, type = "pAdjustMethod")
    .tns.checks(verbose, type = "Verbose")
    .tns.checks(plotpdf, type = "PlotPDF")
    
    if (verbose) 
        message("Transforming into TNA object\n")
    
    if (!is.null(regs))
    {
        if (!all(regs %in% colnames(object@EScores$dif))) 
            stop("'regs' are not present in Enrichment Scores in the TNS object")
    } else
    {
        all.regs <- object@EScores$dif[aSample, ]
        all.regs <- sort(abs(all.regs), decreasing = TRUE)
        if (is.null(ntop))
        {
            regs <- names(all.regs)
        } else
        {
            if (ntop > length(all.regs)) 
                ntop = length(all.regs)
            regs <- names(all.regs[seq_len(ntop)])
        }
    }
    gexp <- object@tni@gexp
    
    ##------ compute reference gx vec
    if (is.null(refsamp))
    {
        gxref <- apply(gexp, 1, mean)
    } else
    {
        idx <- colnames(gexp) %in% refsamp
        if (!all(sum(idx) %in% length(refsamp)))
        {
            stop("NOTE: 'refsamp' should list only valid sample names!")
        }
        gxref <- apply(gexp[, refsamp], 1, mean)
    }
    
    ##-- check log space/transform into log space
    qx <- as.numeric(quantile(gexp, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
    LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 
        1 && qx[4] > 1 && qx[4] < 2)
    if (log || LogC)
    {
        dt <- log2(1 + gexp) - log2(1 + gxref)
        if (LogC) 
            warning("NOTE:'gexp' values seem not to be in log space! ..log2 transformation has been applied!")
    } else
    {
        dt <- gexp - gxref
    }
    
    #-- get phenotype vector
    pheno <- dt[, aSample]
    names(pheno) <- rownames(dt)
    
    #-- make tna from scratch
    rtna <- tni2tna.preprocess(object@tni, phenotype = pheno, verbose = verbose)
    
    #-- do gsea2 analysis for TNA
    if (verbose) 
        message("- Recomputing two-tailed GSEA analysis\n")
    
    rtna <- tna.gsea2(rtna, stepFilter = FALSE, pValueCutoff, pAdjustMethod, tfs = regs, 
        verbose = verbose)
    
    #-- plot
    tna.plot.gsea2(rtna, labPheno = aSample, tfs = regs, plotpdf = plotpdf, ...)
    if (verbose) 
        message("- GSEA plot(s) for the selected sample and regulons should be 
        available at the working directory!\n")
    
}

.tns.stratification <- function(object, nSections = 2)
{
    regstatus <- sign(object@EScores$dif)
    for (reg in colnames(regstatus))
    {
        sq <- c(seq_len(nSections))
        pos <- object@EScores$pos[, reg]
        neg <- object@EScores$neg[, reg]
        dif <- object@EScores$dif[, reg]
        #---
        regstatus[sign(pos) == sign(neg), reg] <- 0
        tp <- regstatus[, reg]
        #---
        tp1 <- sort(dif[tp > 0], decreasing = TRUE)
        tp1[] <- rep(sq, each = ceiling(length(tp1)/nSections), length.out = length(tp1))
        regstatus[names(tp1), reg] <- tp1
        #---
        tp2 <- sort(dif[tp < 0], decreasing = TRUE)
        tp2[] <- rep(sq + nSections + 1, each = ceiling(length(tp2)/nSections), length.out = length(tp2))
        regstatus[names(tp2), reg] <- tp2
    }
    regstatus[regstatus == 0] <- nSections + 1
    object@EScores$regstatus <- regstatus
    object@EScores$mid <- nSections + 1
    
    #-- update
    object@para$strat <- c(nSections = nSections)
    return(object)
    
}

.survplot <- function(EScores, dt, reg, fname, fpath, ylab, xlab, pal, widths, plotpdf, 
    excludeMid, flipcols, attribs, groups, endpoint, dES.ylab, show.KMlegend, KMlegend.pos, 
    KMlegend.cex, show.pval, pval.cex, pval.pos)
    {
    
    #-- organzing data
    tumours <- rev(sort(EScores$dif[, reg], decreasing = TRUE))
    regstatus <- EScores$regstatus[names(tumours), reg]
    nclass <- length(unique(regstatus))
    
    #-- organzing plot colors
    if (pal %in% c("red", "blue", "redblue"))
    {
        if (pal == "red")
        {
            cols <- pal1(nclass)
        } else if (pal == "blue")
        {
            cols <- pal2(nclass)
        } else if (pal == "redblue")
        {
            cols <- pal3(nclass)
        }
        if (flipcols) 
            cols <- rev(cols)
    } else (cols <- pal)
    #--- adjusting for right number of colors
    if (nclass%%2 == 0)
    {
        rmc <- (nclass/2) + 1
        cols <- cols[-rmc]
    }
    #--- adjusting graphical parameters
    op <- par(no.readonly = TRUE)
    np <- length(tumours)
    nms <- pretty(c(1, np), eps.correct = 1) + 1
    dp <- nms[2] - nms[1]
    pp <- length(nms)
    
    if ((abs(nms[pp] - np)/dp) > 0.6) 
        nms <- nms[-pp]
    
    nms[length(nms)] <- np
    
    panels <- c(TRUE, !is.null(attribs), TRUE)
    layout(matrix(seq_len(sum(panels)), 1, sum(panels)), widths = widths[panels])
    par(mgp = c(2.5, 0.4, 0), mar = c(6.5, 5, 3, 0.7))
    xlim <- range(tumours) + c(-0.5, 0.5)
    
    #--- first panel plot (Sample statification)
    barplot(tumours, space = 0, xlim = c(-2, 2), axes = FALSE, cex.lab = 1.2, col = cols[as.factor(regstatus)], 
        hor
 = TRUE, border = NA, axisnames = FALSE, ylab = dES.ylab, xlab = "", 
        beside = TRUE, lwd = 1)
    mtext("Enrichment score\n( dES )", 1, adj = 0.5, line = 3, cex = 0.8)
    mtext(reg, 3, adj = 0.1, line = -0.5, cex = 0.8)
    axis(2, at = nms, labels = nms, tcl = -0.2, las = 2, lwd = 1.8, cex.axis = 1.2)
    axis(1, tcl = -0.2, lwd = 1.8, cex.axis = 1.2)
    #--- second panel plot (Covariable status, optional)
    if (!is.null(attribs))
    {
        attribs <- attribs[names(regstatus), ]
        par(mar = c(7.1, 0, 3.8, 0))
        image(t(attribs), col = c("grey95", "black"), axes = FALSE)
        labs <- colnames(attribs)
        axis(1, at = seq(0, 1, length.out = length(labs)), labels = labs, tcl = -0.2, 
            las = 2, lwd = 1.8, cex.axis = 0.8)
        if (!is.null(groups))
        {
            lanes <- cumsum(groups)[-length(groups)]
            pos <- seq(0, 1, length.out = length(labs))
            pos <- (pos[-length(pos)] + (pos[2:length(pos)]))/2
            par(xpd = TRUE)
            for (i in lanes)
            {
                lines(x = c(pos[i], pos[i]), y = c(-0.2, 1), col = "grey40", lwd = 1.2, 
                  lty = "11", lend = 1)
            }
            par(xpd = FALSE)
        }
    }
    #--- third panel plot (Kaplan-Meier)
    if (excludeMid && nclass%%2 != 0 && nclass > 1)
    {
        rmc <- (nclass + 1)/2
        cols <- cols[-rmc]
        idx <- regstatus != rmc
        regstatus <- regstatus[idx]
        tumours <- tumours[idx]
        nclass <- nclass - 1
    }
    sections <- sort(unique(regstatus))
    if (length(sections) < length(cols)) 
        cols <- cols[-((length(cols) + 1)/2)]
    ddt <- dt[names(regstatus), ]
    ddt$class <- regstatus
    #-- survival analysis
    res1 <- survfit(Surv(time, event) ~ class, data = ddt)
    par(mar = c(6.5, 5, 3, 1))
    plot(res1, col = cols, lwd = 1.8, axes = FALSE, cex.lab = 1.2, cex = 0.5, mark.time = TRUE, 
        ylab = ylab, xlab = "")
    mtext(xlab, 1, adj = 0.5, line = 2, cex = 0.8)
    labs <- as.integer(seq(0, endpoint, length.out = 4))
    if (!endpoint %in% labs) 
        labs <- pretty(c(0, endpoint))
    axis(1, at = labs, labels = labs, tcl = -0.2, las = 1, lwd = 1.8, cex.axis = 1.2)
    axis(2, tcl = -0.2, las = 2, lwd = 1.8, cex.axis = 1.2)
    #---log-rank test
    if (nclass > 1)
    {
        res2 <- survdiff(Surv(time, event) ~ class, data = ddt)
        pval <- 1 - pchisq(res2$chisq, length(res2$n) - 1)
        #---legends
        if (nclass == 2)
        {
            legs <- paste(c("Positive dES", "Negative dES")[1:length(sections)], 
                ": ", res2$n, " (", res2$obs, ")", sep = "")
        } else if (nclass == 3)
        {
            legs <- paste(c("Positive dES", "undetermined", "Negative dES")[1:length(sections)], 
                ": ", res2$n, " (", res2$obs, ")", sep = "")
        } else
        {
            legs <- paste("Section ", 1:length(sections), ": ", res2$n, "(", res2$obs, 
                ")", sep = "")
        }
        pval <- paste("Logrank P: ", format(pval, digits = 3, scientific = TRUE))
        if (show.KMlegend)
        {
            legend(KMlegend.pos, legend = legs, col = cols, bty = "n", pch = 15, 
                cex = KMlegend.cex, pt.cex = 1.5)
        }
        if (show.pval)
        {
            legend(pval.pos, cex = pval.cex, legend = pval, bty = "n", adj = c(0, 
                -0.5))
        }
    }
    par(op)
}
pal1 <- function(nclass)
{
    pt <- rev(colorRampPalette(brewer.pal(9, "Reds"))(11))
    if (nclass == 1)
    {
        cols <- "grey"
    } else if (nclass <= 3)
    {
        cols <- c(pt[c(1)], "grey", pt[c(5)])
    } else if (nclass <= 5)
    {
        cols <- c(pt[c(1, 4)], "grey", pt[c(7, 9)])
    } else if (nclass <= 7)
    {
        cols <- c(pt[c(1, 3, 5)], "grey", pt[c(7, 9, 10)])
    } else
    {
        warning("NOTE: please, provide up to 3 classes (7 sections) for 
                stratification!")
        cols <- "grey"
    }
    cols
}
pal2 <- function(nclass)
{
    pt <- rev(colorRampPalette(brewer.pal(9, "Blues"))(11))
    if (nclass == 1)
    {
        cols <- "grey"
    } else if (nclass <= 3)
    {
        cols <- c(pt[c(1)], "grey", pt[c(5)])
    } else if (nclass <= 5)
    {
        cols <- c(pt[c(1, 4)], "grey", pt[c(7, 9)])
    } else if (nclass <= 7)
    {
        cols <- c(pt[c(1, 3, 5)], "grey", pt[c(7, 9, 10)])
    } else
    {
        warning("NOTE: please, provide up to 3 classes (7 sections) for 
                stratification!")
        cols <- "grey"
    }
    cols
}
pal3 <- function(nclass)
{
    ptreds <- rev(colorRampPalette(brewer.pal(9, "Reds"))(11))
    ptblues <- rev(colorRampPalette(brewer.pal(9, "Blues"))(11))
    if (nclass == 1)
    {
        cols <- "grey"
    } else if (nclass <= 3)
    {
        cols <- c(ptreds[c(4)], "grey", rev(ptblues[c(5)]))
    } else if (nclass <= 5)
    {
        cols <- c(ptreds[c(3, 6)], "grey", rev(ptblues[c(3, 6)]))
    } else if (nclass <= 7)
    {
        cols <- c(ptreds[c(2, 5, 8)], "grey", rev(ptblues[c(2, 5, 8)]))
    } else
    {
        warning("NOTE: please, provide up to 3 classes (7 sections) for 
                stratification!")
        cols <- "grey"
    }
    cols
}

.plotCox = function(resall, regs, keycovar, filen, width, height, xlim, xlab, ylab, 
    plotpdf)
    {
    #--- get colors
    # ptreds<-rev(colorRampPalette(brewer.pal(9,'Reds'))(11))
    # ptblues<-rev(colorRampPalette(brewer.pal(9,'Blues'))(11))
    # pal<-c('black',ptblues[5],'grey60',ptreds[4])
    pal <- c("black", "#008080ff", "grey60", "#d45500ff")
    
    #--- assign colors to line representation based on significance
    cols <- rep(NA, nrow(resall))
    names(cols) <- rownames(resall)
    cols[keycovar] <- 1
    idx1 <- resall[regs, 1] <= 1
    idx2 <- (resall[regs, 3] < 1 & resall[regs, 4] > 1) | (resall[regs, 3] > 1 & 
        resall[regs, 4] < 1)
    cols[names(idx1)][idx1 & !idx2] <- 2
    cols[names(idx1)][idx2] <- 3
    cols[names(idx1)][!idx1 & !idx2] <- 4
    resall <- cbind(resall, cols = cols)
    cols[] <- pal[cols]
    
    #--- get name labels for the graph
    labs <- rownames(resall)
    mxchar <- max(nchar(labs))
    len <- max(mxchar/10, 1)
    
    #---
    urd <- function(d, x)
    {
        lxd <- log10(x/d)
        rlxd <- unique(c(floor(lxd), ceiling(lxd)))
        d * 10^rlxd
    }
    xlabs <- .prettylog(xlim)
    xlim <- range(xlim)
    
    #--- plot graph
    nIN = 2
    nOUT = nrow(resall) + 1
    op <- par(no.readonly = TRUE)
    ylim <- c(0, nOUT + 1)
    if (plotpdf) 
        pdf(file = filen, width = width, height = height)
    par(mai = c(0.4, 1.2 * len, 1, 0.7), mgp = c(3, 0.5, 0), yaxs = "i", xaxs = "i")
    plot(NA, log = "x", xlim = xlim, ylim = ylim, axes = FALSE, ylab = "", xlab = "")
    segments(xlim[1], nIN:nOUT, resall[, 3], col = "grey85", lwd = 1.5, lty = "21", 
        lend = 2)
    # abline(v=1, lty=2,lwd=1.5,col='grey60', lend=2)
    lines(x = c(1, 1), y = c(1, nOUT + 1), lwd = 1.5, col = "grey60", lty = "21", 
        lend = 2)
    segments(resall[, 3], nIN:nOUT, resall[, 4], col = cols, lwd = 1.5)
    points(y = nIN:nOUT, x = resall[, 1], pch = 18, cex = 1.5, lwd = 1, col = cols)
    axis(3, lwd = 2, cex.axis = 1.2, tck = -0.02, labels = xlabs$labs, at = xlabs$at)
    mtext(xlab, side = 3, line = 2, cex = 1.2)
    mtext(ylab, side = 2, line = 1 + 3.5 * len, cex = 1.2)
    labs <- rownames(resall)
    labs[labs == "LN"] <- "LN+"
    mtext(text = labs, side = 2, at = nIN:nOUT, las = 2, cex = 0.8)
    par(xpd = TRUE)
    legend(y = 0, x = 1, ncol = 2, legend = c("covariates", "associated, HR<1", "not associated", 
        "associated, HR>1"), col = pal, pch = 18, lwd = 1.2, xjust = 0.5, yjust = 0.5, 
        bty = "n", cex = 0.65, pt.cex = 0.9)
    
    if (plotpdf)
    {
        message("NOTE: 'PDF' file was generated")
        dev.off()
    }
    par(op)
    invisible(resall)
    
}
.prettylog <- function(xlim)
{
    urd <- function(d, x)
    {
        lxd <- log10(x/d)
        rlxd <- unique(c(floor(lxd), ceiling(lxd)))
        d * 10^rlxd
    }
    xlim <- sort(xlim)
    if (length(xlim) == 2)
    {
        tp <- seq(xlim[1], xlim[2], by = xlim[1])
        tp <- sort(unique(c(1, xlim, tp)))
        xlabs <- sort(urd(1, tp))
        #---
        tp <- unlist(sapply(1:(length(xlabs) - 1), function(i)
        {
            pretty(xlabs[i:(i + 1)], n = 10)
        }))
        tp <- unique(as.numeric(tp))
        tp <- tp[findInterval(tp, xlim) == 1]
        xat <- sort(unique(c(1, xlim, tp)))
        tp <- range(xat)
        if (tp[1] > min(xlabs)) 
            xlabs <- c(xlabs, tp[1])
        if (tp[2] < max(xlabs)) 
            xlabs <- c(xlabs, tp[2])
        xlabs <- sort(xlabs)
    } else
    {
        xat <- xlabs <- xlim
    }
    idx <- xat >= 10
    tp <- c(format(xat[!idx], scientific = FALSE), format(xat[idx], nsmall = 0))
    tp[!xat %in% xlabs] <- NA
    tp[xat == 1] <- "1.0"
    list(labs = tp, at = xat)
}

tns.set <- function(object, para = NULL, what)
{
    if (what == "status-1")
    {
        object@status["Preprocess"] <- "[x]"
        object@status["GSEA2"] <- "[ ]"
        object
    } else if (what == "status-2")
    {
        object@status["GSEA2"] <- "[x]"
        object
    } else if (what == "EScores")
    {
        object@EScores <- para
        object
    } else if (what == "survivalData")
    {
        object@survivalData <- para
        object
    }
    
    
}


