
############################
### Plot PNG panel
############################

pngPanel <- function(dual, path, res) {
    #-- regulon names
    regs <- unlist(strsplit(dual, "~"))
    
    #-- import png images
    png.names <- c(paste0(".1.dESregPlot_", regs[1], ".png"), 
                   paste0(".2.dESregPlot_", regs[2], ".png"),
                   paste0(".4.KMplot_", dual, ".png"),
                   paste0(".3.rankScatter_", dual, ".png"),
                   paste0(".5.coxplot_", dual, ".png"))
    img.loc <- paste0(path, "/", png.names)
    img1 <- readPNG(img.loc[1])
    img2 <- readPNG(img.loc[2])
    img3 <- readPNG(img.loc[3])
    img4 <- readPNG(img.loc[4])
    img5 <- readPNG(img.loc[5])
    
    #-- layout for the panel
    lmat <- matrix(c(0,1,1,0,3,0,
                     0,1,1,0,3,0,
                     2,2,4,0,3,0,
                     2,2,4,0,5,0,
                     2,2,0,0,5,0), ncol = 6, byrow = TRUE)
    
    #-- parameters for the plot
    png(filename = paste0(path, "/", dual, ".png"), width = 11.7*res, height = 8.26*res)
    layout(lmat, heights = c(1,0.15,0.5,0.5,0.15), widths = c(1,0.2,1,0.05,1.1,0.1))
    par(xaxs="i", yaxs="i", mar = c(0.5,1,1.5,1.5))
    
    #-- dES1
    plot(seq(0,dim(img1)[1], length.out = 5), seq(0,dim(img1)[2], length.out = 5),
         type = "n", ylab = "", xlab = "", axes = FALSE)
    rasterImage(img1, 0,0, dim(img1)[1], dim(img1)[2])
    
    #-- dES2
    plot(seq(0,dim(img2)[1], length.out = 5), seq(0,dim(img2)[2], length.out = 5),
         type = "n", ylab = "", xlab = "", axes = FALSE)
    rasterImage(img2, 0,0, dim(img2)[1], dim(img2)[2])
    
    #-- KM
    plot(seq(0,dim(img3)[1], length.out = 5), seq(0,dim(img3)[2], length.out = 5),
         type = "n", ylab = "", xlab = "", axes = FALSE)
    rasterImage(img3, 0,0, dim(img3)[1], dim(img3)[2])
    
    #-- rankScatter
    plot(seq(0,dim(img4)[1], length.out = 5), seq(0,dim(img4)[2], length.out = 5),
         type = "n", ylab = "", xlab = "", axes = FALSE)
    rasterImage(img4, 0,0, dim(img4)[1], dim(img4)[2])
    
    #-- Cox
    plot(seq(0,dim(img5)[1], length.out = 5), seq(0,dim(img5)[2], length.out = 5),
         type = "n", ylab = "", xlab = "", axes = FALSE)
    rasterImage(img5, 0,0, dim(img5)[1], dim(img5)[2])
    
    dev.off()
    
    #-- delete the pngs
    sapply(img.loc, file.remove)
}

################ dESregPlot
########### An adaptation of tnsKM plot, for the first 2 panels

dESregPlot <- function (tns, regs = NULL, attribs = NULL, nSections = 2,
                        pal = "redblue", excludeMid = FALSE, flipcols = FALSE, 
                        panelWidths = c(2, 3), flipGraph = FALSE,
                        xname = NULL, attribs.cex = 1) {
    
    dES.ylab = "Sample ranking"
    if (is.null(xname))
        xname <- regs
    
    #-- stratification
    tns <- .tns.stratification(tns, nSections = nSections)
    
    #-- organize data
    survData <- tnsGet(tns, what = "survivalData")
    EScores <- tnsGet(tns, what = "EScores")
    #-- 
    sp1 <- rownames(survData)
    sp2 <- rownames(EScores$regstatus)
    survData <- survData[sp1 %in% sp2, , drop = FALSE]
    EScores$pos <- EScores$pos[sp2 %in% sp1, , drop = FALSE]
    EScores$neg <- EScores$neg[sp2 %in% sp1, , drop = FALSE]
    EScores$dif <- EScores$dif[sp2 %in% sp1, , drop = FALSE]
    
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
    
    #---plot
    .dESplot(EScores, survData, regs, attribs, pal, panelWidths, excludeMid, 
             flipcols, groups, dES.ylab = dES.ylab, flipGraph = flipGraph, xname,
             attribs.cex)
    par(mfrow=c(1,1))
    invisible(EScores)
}

.dESplot <- function(EScores, dt, reg, attribs, pal, widths,
                     excludeMid, flipcols, groups, dES.ylab, flipGraph,
                     xname, attribs.cex)
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
    np <- length(tumours)
    nms <- pretty(c(1, np), eps.correct = 1) + 1
    dp <- nms[2] - nms[1]
    pp <- length(nms)
    
    if ((abs(nms[pp] - np)/dp) > 0.6) 
        nms <- nms[-pp]
    
    nms[length(nms)] <- np
    
    panels <- c(!is.null(attribs), TRUE)
    if(flipGraph)
        layout(matrix(seq_len(sum(panels))), heights = widths[panels])
    else
        layout(matrix(seq_len(sum(panels)), 1, sum(panels)), widths = widths[panels])
    
    
    #-- covar plot
    if (!is.null(attribs))
    {
        attribs <- attribs[names(regstatus), ]
        labs <- colnames(attribs)
        if (flipGraph){
            par(mar = c(0, 5, 3, 1))
            image(attribs, col = c("grey95", "black"), axes = FALSE)
            axis(2, at = seq(0, 1, length.out = length(labs)), labels = labs, tcl = -0.2, 
                 las = 2, lwd = 1.8, cex.axis = attribs.cex, mgp = c(3,0.3,0))
            axis(3, at = seq(0,1,by = 1/(length(nms) - 1)), labels = nms, tcl = -0.2, las = 2, lwd = 1.8, cex.axis = 0.8,
                 mgp = c(3,0.1,0), las = 1)
            mtext(dES.ylab, 3, line = 1.5, cex = 1.2)
            if(xname != reg) {
                mtext(xname, 3, adj = -0.1, line = 1.5, cex = 2, font = 2)
            }
            else {
                mtext(xname, 3, adj = -0.45, line = 1.5, cex = 1.7, font = 2)
            }
            
            
        }
        else {
            par(mar = c(5, 3.5, 1.9, 0))
            image(t(attribs), col = c("grey95", "black"), axes = FALSE)
            axis(1, at = seq(0, 1, length.out = length(labs)), labels = labs, tcl = -0.2, 
                 las = 2, lwd = 1.8, cex.axis = attribs.cex, mgp = c(3,0.3,0))
            axis(2, at = seq(0,1,by = 1/(length(nms) - 1)), labels = nms, tcl = -0.2, las = 2, lwd = 1.8, cex.axis = 0.8,
                 mgp = c(3,0.3,0))
            mtext(dES.ylab, 2, line = 1.8, cex = 1.2)
            if(xname != reg) {
                mtext(xname, 3, adj = -0.4, line = 0.3, cex = 2, font = 2)
            } else {
                mtext(xname, 3, adj = 0.2, line = 0.3, cex = 1.7, font = 2)
            }
            
        }
        
        if (!is.null(groups))
        {
            lanes <- cumsum(groups)[-length(groups)]
            pos <- seq(0, 1, length.out = length(labs))
            pos <- (pos[-length(pos)] + (pos[2:length(pos)]))/2
            par(xpd = TRUE)
            for (i in lanes)
            {
                if (flipGraph)
                    lines(y = c(pos[i], pos[i]), x = c(-0.2, 1), col = "grey40", lwd = 1.2, 
                          lty = "11", lend = 1)
                else
                    lines(x = c(pos[i], pos[i]), y = c(-0.2, 1), col = "grey40", lwd = 1.2, 
                          lty = "11", lend = 1)
            }
            par(xpd = FALSE)
        }
    }
    
    #--- sample stratification
    if (flipGraph){
        par(mgp = c(2, 0.4, 0), mar = c(1, 5, 1, 1), xaxs="i")
        barplot(tumours, space = 0, axes = FALSE, cex.lab = 1, col = cols[as.factor(regstatus)], 
                border = NA, axisnames = FALSE, ylab = "", xlab = "", beside = TRUE, lwd = 1, ylim = c(-2, 2))
        mtext("Enrichment score", 2, line = 1.5, adj = 0.5, cex = 1.2)
        axis(2, at = seq(-2, 2, 1),tcl = -0.2, lwd = 1.8, cex.axis = 1, las = 2)
    }
    else {
        par(mgp = c(2.5, 0.4, 0), mar = c(5, 1, 1.9, 0.7), xaxs="i", yaxs = "i")
        barplot(tumours, space = 0, xlim = c(-2, 2), axes = FALSE, cex.lab = 1, col = cols[as.factor(regstatus)], 
                hor = TRUE, border = NA, axisnames = FALSE, ylab = "", xlab = "", beside = TRUE, lwd = 1)
        mtext("Enrichment score", 1, adj = 0.5, line = 1.5, cex = 1.2)
        axis(1, tcl = -0.2, lwd = 1.8, cex.axis = 1)
    }
}

#--- Dual cox plot function
dualCoxPlot <- function(dual, dualCoxTable, ...) {
    res.table <- matrix(dualCoxTable, ncol = 3, byrow = TRUE)
    res.table <- cbind(res.table[,1], NA, res.table[,2:3])
    
    regs <- unlist(strsplit(dual, "~"))
    
    dimnames(res.table) = list(c(regs[1], regs[2],dual), 
                               c("exp(coef)", "exp(-coef)", "lower .95", 
                                 "upper .95"))
    
    .plotCox_mod(res.table, xlim = c(0.3, 6), plotpdf = FALSE,
                 xlab = "Hazard ratio", ...)
}

.plotCox_mod = function(resall, regs, keycovar, filen, width, height, xlim, xlab,
                        plotpdf, xname = NULL)
{
    
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
    labs[3] <- "Interaction"
    xlabs <- .prettylog(xlim)
    xlim <- range(xlim)
    
    #--- plot graph
    nIN = 2
    nOUT = nrow(resall) + 1
    op <- par(no.readonly = TRUE)
    ylim <- c(0, nOUT + 1)
    if (plotpdf) 
        pdf(file = filen, width = width, height = height)
    par(mai = c(0.4, 1.5, 0.8, 0.7), mgp = c(3, 0.5, 0), yaxs = "i", xaxs = "i")
    plot(NA, log = "x", xlim = xlim, ylim = ylim, axes = FALSE, ylab = "", xlab = "")
    segments(xlim[1], nIN:nOUT, resall[, 3], col = "grey85", lwd = 1.5, lty = "21", 
             lend = 2)
    lines(x = c(1, 1), y = c(1, nOUT + 1), lwd = 1.5, col = "grey60", lty = "21", 
          lend = 2)
    segments(resall[, 3], nIN:nOUT, resall[, 4], col = cols, lwd = 1.5)
    points(y = nIN:nOUT, x = resall[, 1], pch = 18, cex = 1.2, lwd = 1, col = cols)
    axis(3, lwd = 2, cex.axis = 1, tck = -0.02, labels = xlabs$labs, at = xlabs$at)
    mtext(xlab, side = 3, line = 1.5, cex = 1.2)
    if (!is.null(xname)) {
        mtext(xname, 3, line = 2.5, adj = -0.4, cex = 2, font = 2)
    }
    labs[labs == "LN"] <- "LN+"
    mtext(text = labs, side = 2, at = nIN:nOUT, las = 2, cex = 1)
    par(xpd = TRUE)
    legend(y = 0, x = 1, ncol = 2, legend = c("covariates", "associated, HR<1", "not associated", 
                                              "associated, HR>1"), col = pal, pch = 18, lwd = 1.2, xjust = 0.5, yjust = 0.5, 
           bty = "n", cex = 0.65, pt.cex = 0.9)
    
    if (plotpdf)
    {
        message("NOTE: a 'PDF' file should be available at the working directory!\n")
        dev.off()
    }
    par(op)
    invisible(resall)
    
}

KMregPlot <- function (tns, reg, nSections = 2, endpoint = 60, ylab = "Survival probability", 
                       pal = "redblue", xlab = "Months", excludeMid = FALSE,  show.KMlegend = TRUE, 
                       KMlegend.pos = "bottomleft", KMlegend.cex = 1, show.pval = TRUE, 
                       pval.cex = 1, pval.pos = "topright", flipcols = FALSE, y.axis = TRUE, samples = NULL, title = NA,
                       ylab.cex = 1, xlab.cex = 1, title.cex = 1, sectionsLegend = NULL,
                       xname = NULL) {
    tns <- .tns.stratification(tns, nSections = nSections)
    
    #-- title
    if(is.na(title))
        title <- reg
    
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
    
    if(!is.null(samples))
        survData <- survData[samples,]
    
    #---plot
    .KMplot(survData, reg, EScores, endpoint, ylab, 
            xlab, excludeMid,  show.KMlegend, 
            KMlegend.pos, KMlegend.cex, show.pval, 
            pval.cex, pval.pos, pal, flipcols, y.axis, samples, title, ylab.cex, 
            xlab.cex, title.cex, sectionsLegend, xname)
    invisible(EScores)
}


.KMplot <- function (dt, reg, EScores, endpoint, ylab, 
                     xlab, excludeMid,  show.KMlegend, 
                     KMlegend.pos, KMlegend.cex, show.pval, 
                     pval.cex, pval.pos, pal, flipcols, y.axis, samples, title,
                     xlab.cex, ylab.cex, title.cex, sectionsLegend, xname) {
    #-- organzing data
    tumours <- rev(sort(EScores$dif[, reg], decreasing = TRUE))
    regstatus <- EScores$regstatus[names(tumours), reg]
    nclass <- length(unique(regstatus))
    
    if (any(pal %in% c("red", "blue", "redblue")))
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
    
    #--- Start of KM plot
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
    if (y.axis & ylab != "") {
        if (title.cex > 1) {
            par(mar = c(1.5, 3, 0, 0), mgp = c(2, 0.3, 0), mai = c(0.7,0.7,0.2,0.1))
            plot(res1, col = cols, lwd = 1.8, axes = FALSE, cex.lab = ylab.cex, cex = 1, mark.time = TRUE, 
                 ylab = ylab, xlab = "")
        }
        else {
            par(mar = c(1.5, 3, 4, 0), mgp = c(2, 0.3, 0), mai = c(0.4,0.7,0.4,0.1))
            plot(res1, col = cols, lwd = 1.8, axes = FALSE, cex.lab = ylab.cex, cex = 1, mark.time = TRUE, 
                 ylab = ylab, xlab = "")
        }
    }
    else {
        par(mar = c(1.5, 0, 4, 0), mgp = c(2, 0.3, 0), mai = c(0.4,0.3,0.4,0.1))
        plot(res1, col = cols, lwd = 1.8, axes = FALSE, cex = 1, mark.time = TRUE, 
             ylab = "", xlab = "")
    }
    if (title.cex > 1)
        mtext(title, 3, cex = title.cex, adj = 1, line = 0)
    else
        mtext(title, 3, cex = title.cex, adj = 1, line = 0)
    mtext(xlab, 1, adj = 0.5, line = 1.2, cex = xlab.cex)
    if(!is.null(xname)) {
        mtext(xname, 3, adj = -0.1, line = 0.7, cex = 2, font = 2)
    }
    
    
    labs <- as.integer(seq(0, endpoint, length.out = 4))
    if (!endpoint %in% labs) 
        labs <- pretty(c(0, endpoint))
    axis(1, at = labs, labels = labs, tcl = -0.2, las = 1, lwd = 1.8, cex.axis = 1)
    if (y.axis)
        axis(2, tcl = -0.2, las = 2, lwd = 1.8, cex.axis = 1)
    #---log-rank test
    if (nclass > 1)
    {
        res2 <- survdiff(Surv(time, event) ~ class, data = ddt)
        pval <- 1 - pchisq(res2$chisq, length(res2$n) - 1)
        #---legends
        if (is.null(sectionsLegend))
        {
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
        }
        else {
            legs <- paste(sectionsLegend, ": ", res2$n, "(", res2$obs, ")", sep = "")
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
}

rankScatter <- function(tns1, tns2, dual, dES1, dES2, pal = "BrBG", nSections, 
                        mode = "agreement") {
    regs <- unlist(strsplit(dual, "~"))
    
    #-- get regulon activity
    dES_reg1 <- dES1$dif[,regs[1]]
    dES_reg2 <- dES2$dif[,regs[2]]
    
    #-- re-arrange samples
    dES_reg2 <- dES_reg2[names(dES_reg1)]
    
    #-- assign colors to points
    samples <- names(dES_reg1)
    if (nSections > 2)
        cols <- RColorBrewer::brewer.pal(nSections*4-1, pal)
    else
        cols <- RColorBrewer::brewer.pal(nSections*4, pal)
    
    if (mode == "agreement"){
        
        bgcol <- sapply(samples, function(sample) {
            #-- get sample stratums
            strat1 <- dES1$regstatus[sample,regs[1]]
            strat2 <- dES2$regstatus[sample, regs[2]]
            
            if(strat1 == strat2) {
                if (strat1 <= nSections) {
                    reg.col <- cols[strat1+1]
                    b.col <- cols[1]
                }
                else if (strat1 == nSections+1) {
                    reg.col <- "grey80"
                    b.col <- "grey60"
                }
                else {
                    if(nSections == 1)
                        reg.col <- cols[strat1]
                    else
                        reg.col <- cols[strat1+2]
                    b.col <- cols[length(cols)]
                }
            }
            else {
                reg.col <- "white"
                b.col <- "grey60"
            }
            c(reg.col, b.col)
            
        })
    }
    else if (mode == "disagreement") {
        cols <- cols
        bgcol <- sapply(samples, function(sample) {
            #-- get sample stratums
            strat1 <- dES1$regstatus[sample,regs[1]]
            strat2 <- dES2$regstatus[sample, regs[2]]
            
            up_strats <- 1:(nSections*2+1)
            down_strats <- (nSections*2+1):1
            
            if(down_strats[strat1] == strat2) {
                if (strat1 <= nSections) {
                    reg.col <- cols[strat1+1]
                    b.col <- cols[1]
                }
                else if (strat1 == nSections+1) {
                    reg.col <- "grey80"
                    b.col <- "grey60"
                }
                else {
                    if(nSections == 1)
                        reg.col <- cols[strat1]
                    else
                        reg.col <- cols[strat1+2]
                    b.col <- cols[length(cols)]
                }
            }
            
            
            else {
                reg.col <- "white"
                b.col <- "grey60"
            }
            c(reg.col, b.col)
        })
    }
    
    #--- sample position in plot
    samples.pos <- sapply(samples, function(sample){
        reg1ac <- dES_reg1[sample]
        reg2ac <- dES_reg2[sample]
        
        #-- get sample place in samples
        s.order1 <- names(sort(dES_reg1))
        x <- which(sample == s.order1)
        
        s.order2 <- names(sort(dES_reg2))
        y <- which(sample == s.order2)
        c(x,y)
    })
    
    #-- axis
    np <- length(samples)
    nms <- pretty(c(1, np), eps.correct = 1) + 1
    dp <- nms[2] - nms[1]
    pp <- length(nms)
    
    if ((abs(nms[pp] - np)/dp) > 0.6) 
        nms <- nms[-pp]
    
    nms[length(nms)] <- np
    
    #--- layout
    layout(matrix(c(0,1,2,
                    3,0,0), nrow = 3, byrow = TRUE), heights = c(1.4,10,1.3),
           widths = c(1.3,10))
    limit <- length(samples)*0.02
    
    #--- upper color bar
    par(mar = c(0,0,2,1.5))
    if (any(dES1$regstatus[,regs[1]] == nSections+1))
        stratcols <- rev(as.matrix(RColorBrewer::brewer.pal(nSections*2+1,"RdBu")))
    else 
        stratcols <- rev(as.matrix(RColorBrewer::brewer.pal(nSections*2,"RdBu")))
    sections <- sort(unique(dES1$regstatus[, regs[1]]), decreasing = TRUE)
    x.pos <- rep(0, length(sections)-1)
    j <- 1
    for (i in sections[1:(length(sections)-1)]) {
        x.pos[j] <- line.pos(dES1, dES_reg1, i, regs[1])
        j <- j+1
    }
    
    col.pos <- rep(0, length(sections))
    
    col.pos[1] <- x.pos[1]+limit
    col.pos[length(col.pos)] <- ncol(samples.pos)-x.pos[length(x.pos)]+limit
    
    for (i in 2:(length(col.pos)-1)) {
        col.pos[i] <- x.pos[i] - x.pos[i-1]
    }
    
    
    x <- 1:ncol(samples.pos)
    
    col.mat <- matrix(rep(stratcols, col.pos), ncol = 1)
    
    z <- matrix(x, ncol = 1)
    
    image(x, y = 1, z, col = col.mat, axes=FALSE,xlab="",ylab="")
    mtext(regs[1], 3, cex = 1.5, font = 2, line = 0.2)
    
    #--- side color bar
    if (any(dES2$regstatus[,regs[2]] == nSections+1))
        stratcols <- rev(as.matrix(RColorBrewer::brewer.pal(nSections*2+1,"RdBu")))
    else 
        stratcols <- rev(as.matrix(RColorBrewer::brewer.pal(nSections*2,"RdBu")))
    
    par(mar = c(0.5,2,0,0))
    sections <- sort(unique(dES2$regstatus[, regs[2]]), decreasing = TRUE)
    y.pos <- rep(0, length(sections)-1)
    j <- 1
    for (i in sections[1:length(sections)-1]) {
        y.pos[j] <- line.pos(dES2, dES_reg2, i, regs[2])
        j <- j+1
    }
    
    col.pos <- rep(0, length(sections))
    
    col.pos[1] <- y.pos[1]+limit
    col.pos[length(col.pos)] <- ncol(samples.pos)-y.pos[length(y.pos)]+limit
    
    for (i in 2:(length(col.pos)-1)) {
        col.pos[i] <- y.pos[i] - y.pos[i-1]
    }
    
    col.mat <- matrix(rep(stratcols, col.pos), nrow = 1)
    
    y <- 1:ncol(samples.pos)
    z <- matrix(1:ncol(samples.pos), nrow = 1)
    
    image(x = 1, y, z, col = col.mat, axes=FALSE,xlab="",ylab="")
    mtext(regs[2], 2, cex = 1.5, font = 2, line = 0.2)
    
    
    #-- scatter plot
    par(mar = c(0.5,0,0,1.5))
    plot(samples.pos[1,], samples.pos[2,], axes = FALSE, 
         xlab = NA, ylab = NA, xaxs = "i", yaxs = "i", xlim = c(-limit, length(samples)+limit),
         ylim = c(-limit, length(samples)+limit))
    box(lty = "21", col = "grey60", lwd = 1)
    mtext("Interaction", 1, adj = 0, cex = 1.5, font = 2, line = 1)
    
    #-- lines
    for (i in x.pos) {
        abline(v = i, col = "grey60", lty = "21", lwd = 1)
    }
    
    for (i in y.pos) {
        abline(h = i, col = "grey60", lty = "21", lwd = 1) 
    }
    
    #-- points
    points(samples.pos[1,], samples.pos[2,], pch = 21, bg = bgcol[1,], 
           col = bgcol[2,])
    
    #-- legend bar
    # if (mode == "ranking") {
    #     col.mat <- matrix(rep(cols[c((nSections*4-1):(nSections*4-nSections),
    #                                  (nSections+1):2)], nSections*2+1), ncol = (nSections*2+1))
    #     col.mat <- rbind(col.mat[1:nSections,], "grey60", col.mat[(nSections+1):(nSections*2),])
    #     par(mar = c(1.5,3,0,3))
    #     z=matrix(1:((nSections*2+1)), ncol = 1)
    #     x=1:(nSections*2+1)
    #     y=1
    #     image(x,y,z,col=col.mat, axes=FALSE,xlab="",ylab="")
    #     mtext("Activated", 1, adj = 1, cex = 1)
    #     mtext("Repressed",1, adj = 0, cex = 1)
    # }
    # else {
    #     col.mat <- matrix(rep(cols[2:6], ncol = 5))
    #     par(mar = c(2,3,0,3))
    #     z=matrix(1:((nSections*2+1)), ncol = 1)
    #     x=1:(nSections*2+1)
    #     y=1
    #     image(x,y,z,col=col.mat, axes=FALSE,xlab="",ylab="")
    #     mtext("", 1, adj = 1, cex = 1)
    #     mtext("Disagreement", 1, adj = 1, line = 0.5, cex = 1)
    #     mtext("Agreement",1, adj = 0, line = 0.5, cex = 1)
    # }
    
    
    par(mfrow=c(1,1))
    
    invisible(lm(samples.pos[1,] ~ samples.pos[2,]))
}

#--- auxiliary functions
line.pos <- function(dESv, dESregv, strat.pos, reg) {
    strat <- names(which(dESv$regstatus[, reg] == strat.pos))
    sample <- names(which.max(dESregv[strat]))
    s.order <- names(sort(dESregv))
    return(which(sample == s.order))
}






