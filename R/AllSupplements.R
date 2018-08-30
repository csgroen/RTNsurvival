
#' Sample stratification for a TNS object
#'
#' Internal function, used for sample stratification.
#' 
#' @param tns a \linkS4class{TNS} object, which must have passed GSEA2 analysis.
#' @param nSections A numeric value for the stratification of the sample. The 
#' larger the number, the more subdivisions will be created for the Kaplan-Meier 
#' analysis.
#' @param center a logical value. If TRUE, numbers assigned to each group is
#' centralized on regulon activity scale.
#' 
#' @return An updated \linkS4class{TNS} object.
#' @examples
#' # see tnsKM method.
#'
#' @seealso \code{\link{tnsKM}}
#' @aliases tnsStratification
#' @export
#' 
tnsStratification <- function(tns, nSections = 1, center = FALSE){
  
  #-- checks
  .tns.checks(tns, type = "TNS")
  .tns.checks(nSections, type = "nSections")
  .tns.checks(center, type = "center")
  
  #-- stratification
  para <- tnsGet(tns, what = "para")
  if(para$regulonActivity=="gsea2"){
    tns <- .tns.stratification.gsea2(tns, nSections = nSections, center = center)
  } else {
    tns <- .tns.stratification.area(tns, nSections = nSections, center = center)
  }
  
  return(tns)
  
}

##------------------------------------------------------------------------------
.tns.stratification.gsea2 <- function(tns, nSections=1, center=TRUE){
  regstatus <- sign(tns@results$regulonActivity$dif)
  for (reg in colnames(regstatus)){
    sq <- c(seq_len(nSections))
    pos <- tns@results$regulonActivity$pos[, reg]
    neg <- tns@results$regulonActivity$neg[, reg]
    dif <- tns@results$regulonActivity$dif[, reg]
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
  mid <- nSections + 1
  regstatus[regstatus == 0] <- mid
  #---
  if(center){
    regstatus <- -1 * (regstatus - mid)
    mid <- 0
  }
  #---
  tns@results$regulonActivity$regstatus <- regstatus
  tns@results$regulonActivity$nSections <- nSections
  tns@results$regulonActivity$center <- mid
  return(tns)
}
##------------------------------------------------------------------------------
.tns.stratification.area <- function(tns, nSections=1, center=FALSE){
  regstatus <- sign(tns@results$regulonActivity$dif)
  for (reg in colnames(regstatus)){
    sq <- c(seq_len(nSections))
    dif <- tns@results$regulonActivity$dif[, reg]
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
  #--- obs. this stratification generates a 'midle' group 
  #--- only when there are samples with regulon activity assigned with 0 or NA
  mid <- nSections + 1
  regstatus[regstatus == 0 | is.na(regstatus)] <- mid
  #---
  if(center){
    regstatus <- -1 * (regstatus - mid)
    mid <- 0
  }
  #---
  tns@results$regulonActivity$regstatus <- regstatus
  tns@results$regulonActivity$nSections <- nSections
  tns@results$regulonActivity$center <- mid
  return(tns)
}

##------------------------------------------------------------------------------
.checkRegel <- function(tni, regel){
  #---check regel
  tp <- sapply(colnames(tni@rowAnnotation), function(i) {
    sum(regel%in%tni@rowAnnotation[, i])
  })
  colid <- names(tp[which.max(tp)])
  idx <- which(tni@rowAnnotation[, colid]%in%regel)
  if(length(idx) < length(regel)) {
    warning("Not all names in 'regulatoryElements' are available in the 'TNI' rowAnnotation!",
            call.=FALSE)
  }
  tp <- tni@rowAnnotation[idx,]
  idx <- tni@regulatoryElements %in% rownames(tp)
  if(sum(idx)==0){
    tp <- paste("NOTE: no names in 'regulatoryElements' has been used to call ",
                "regulons in the provided 'TNI'!", sep="")
    stop(tp, call.=FALSE)
  } else if(sum(idx) < length(regel)){
    tp <- paste("Not all names in 'regulatoryElements' have been used to call ",
                "regulons in the provided 'TNI'!", sep="")
    warning(tp, call.=FALSE)
  }
  regel <- tni@regulatoryElements[idx]
  if(length(regel)<2){
    tp <- paste("NOTE: at least two valid names in 'regulatoryElements'", 
                " are required to call dual regulons!", sep="")
    stop(tp, call.=FALSE)
  }
  return (regel)
}

##------------------------------------------------------------------------------
.survstats <- function(regulonActivity, survData, reg, excludeMid){
  #-- get data
  tumours <- rev(sort(regulonActivity$dif[, reg], decreasing = TRUE))
  regstatus <- regulonActivity$regstatus[names(tumours), reg]
  nclass <- length(unique(regstatus))
  #--- third panel plot (Kaplan-Meier)
  if (excludeMid && nclass%%2 != 0 && nclass > 1){
    rmc <- (nclass + 1)/2
    cols <- cols[-rmc]
    idx <- regstatus != rmc
    regstatus <- regstatus[idx]
    tumours <- tumours[idx]
    nclass <- nclass - 1
  }
  sections <- sort(unique(regstatus))
  ddt <- survData[names(regstatus), ]
  ddt$class <- regstatus
  
  #---log-rank test
  survtb <- c(ChiSquare=NA, Pvalue=NA)
  survdf <- NA
  survft <- NA
  if(nclass > 1){
    survft <- survfit(Surv(time, event) ~ class, data = ddt)
    survdf <- survdiff(Surv(time, event) ~ class, data = ddt)
    pval <- 1 - pchisq(survdf$chisq, length(survdf$n) - 1)
    survtb[] <- c(survdf$chisq,pval)
  }
  res <- list(kmTable=survtb, survdiff=survdf, survfit=survft)
  return(res)
}

##------------------------------------------------------------------------------
.survplot <- function(regulonActivity, kmtb, survdf, survft, reg, endpoint, 
                      xlab, ylab, colorPalette, panelWidths, 
                      excludeMid, attribs, groups){
  #-- get data
  tumours <- rev(sort(regulonActivity$dif[, reg], decreasing = TRUE))
  regstatus <- regulonActivity$regstatus[names(tumours), reg]
  nclass <- length(unique(regstatus))
  
  #-- get colors
  if (is.singleString(colorPalette)){
    if (colorPalette == "red"){
      cols <- pal1(nclass)
    } else if (colorPalette == "blue"){
      cols <- pal2(nclass)
    } else if (colorPalette %in% c("redblue","bluered")){
      cols <- pal3(nclass)
    }
    if(colorPalette=="redblue")
      cols <- rev(cols)
  } else {
    cols <- colorPalette
  }
  if (nclass%%2 == 0){
    rmc <- (nclass/2) + 1
    cols <- cols[-rmc]
  }
  
  #--- adjust graphical parameters
  op <- par(no.readonly = TRUE)
  np <- length(tumours)
  nms <- pretty(c(1, np), eps.correct = 1) + 1
  dp <- nms[2] - nms[1]
  pp <- length(nms)
  if ((abs(nms[pp] - np)/dp) > 0.6) nms <- nms[-pp]
  nms[length(nms)] <- np
  if(!is.null(attribs)){
    mt <- matrix(c(1,2,3), 1, 3)
    mar = c(6.5, 5, 3, 0.7)
  } else {
    mt <- matrix(c(1,1,2), 1, 3)
    mar = c(6.5, 15, 3, 0)
  }
  layout(mt, widths = panelWidths)
  par(mgp = c(2.5, 0.4, 0), mar = mar)
  xlim <- range(tumours) + c(-0.5, 0.5)
  
  #--- first plot (stratification)
  barplot(tumours, space = 0, xlim = c(-2, 2), axes = FALSE, cex.lab = 1.2, 
          col = cols[as.factor(regstatus)], horiz = TRUE, border = NA, 
          axisnames = FALSE, ylab = "Samples", xlab = "", 
          beside = TRUE, lwd = 1)
  mtext("Regulon activity (dES)", 1, adj = 0.5, line = 2, cex = 0.8)
  mtext(reg, 3, adj = 0.1, line = -0.5, cex = 0.8)
  axis(2, at = nms, labels = nms, tcl = -0.2, las = 2, lwd = 1.8, cex.axis = 1.2)
  axis(1, tcl = -0.2, lwd = 1.8, cex.axis = 1.2)
  
  #--- second plot (attribs, optional)
  if(!is.null(attribs)){
    attribs <- attribs[names(regstatus), ]
    par(mar = c(6.5, 0, 3, 0))
    image(t(attribs), col = c("grey95", "black"), axes = FALSE, ylim=c(-0.04,1.04))
    labs <- colnames(attribs)
    axis(1, at = seq(0, 1, length.out = length(labs)), labels = labs, tcl = -0.2, 
         las = 2, lwd = 1.8, cex.axis = 0.8)
    if (!is.null(groups)){
      lanes <- cumsum(groups)[-length(groups)]
      pos <- seq(0, 1, length.out = length(labs))
      pos <- (pos[-length(pos)] + (pos[2:length(pos)]))/2
      par(xpd = TRUE)
      for (i in lanes){
        lines(x = c(pos[i], pos[i]), y = c(-0.2, 1), col = "grey40", lwd = 1.2, 
              lty = "11", lend = 1)
      }
      par(xpd = FALSE)
    }
  }
  #--- set excludeMid
  if (excludeMid && nclass%%2 != 0 && nclass > 1){
    rmc <- (nclass + 1)/2
    cols <- cols[-rmc]
    idx <- regstatus != rmc
    regstatus <- regstatus[idx]
    tumours <- tumours[idx]
    nclass <- nclass - 1
  }
  sections <- sort(unique(regstatus))
  if (length(sections) < length(cols)) cols <- cols[-((length(cols) + 1)/2)]
  
  #-- km plot
  par(mar = c(6.5, 5, 3, 1))
  plot(survft, col = cols, lwd = 1.8, axes = FALSE, cex.lab = 1.2, cex = 0.5, 
       mark.time = TRUE, ylab = ylab, xlab = "")
  mtext(xlab, 1, adj = 0.5, line = 2, cex = 0.8)
  labs <- as.integer(seq(0, endpoint, length.out = 4))
  if (!endpoint %in% labs) labs <- pretty(c(0, endpoint))
  axis(1, at = labs, labels = labs, tcl = -0.2, las = 1, lwd = 1.8, cex.axis = 1.2)
  axis(2, tcl = -0.2, las = 2, lwd = 1.8, cex.axis = 1.2)
  
  #---legends
  if (nclass > 1){
    if (nclass == 2){
      legs <- paste(c("Positive dES", "Negative dES")[1:length(sections)], ": ", 
                    survdf$n, " (", survdf$obs, ")", sep = "")
    } else if (nclass == 3) {
      legs <- paste(c("Positive dES", "undetermined", "Negative dES")[1:length(sections)], ": ", 
                    survdf$n, " (", survdf$obs, ")", sep = "")
    } else {
      legs <- paste("Section ", 1:length(sections), ": ", survdf$n, 
                    "(", survdf$obs,")", sep = "")
    }
    legend("bottomleft", legend = legs, col = cols, bty = "n", pch = 15, 
           cex = 1, pt.cex = 1.5)
    par(xpd=TRUE)
    pval <- kmtb$Adjusted.Pvalue
    pval <- paste("Logrank P: ", format(pval, digits = 3, scientific = TRUE))
    legend("topright", cex = 1, legend = pval, bty = "n", inset = c(0,-0.05))
  }
  par(op)
}

##------------------------------------------------------------------------------
pal1 <- function(nclass){
  pt <- rev(colorRampPalette(brewer.pal(9, "Reds"))(11))
  if (nclass == 1){
    cols <- "grey80"
  } else if(nclass <= 3)
  {
    cols <- c(pt[c(1)], "grey80", pt[c(5)])
  } else if(nclass <= 5){
    cols <- c(pt[c(1, 4)], "grey80", pt[c(7, 9)])
  } else if(nclass <= 7){
    cols <- c(pt[c(1, 3, 5)], "grey80", pt[c(7, 9, 10)])
  } else {
    warning("NOTE: please, provide up to 3 sections for stratification!")
    cols <- "grey80"
  }
  cols
}

##------------------------------------------------------------------------------
pal2 <- function(nclass){
  pt <- rev(colorRampPalette(brewer.pal(9, "Blues"))(11))
  if (nclass == 1){
    cols <- "grey80"
  } else if (nclass <= 3){
    cols <- c(pt[c(1)], "grey80", pt[c(5)])
  } else if (nclass <= 5){
    cols <- c(pt[c(1, 4)], "grey80", pt[c(7, 9)])
  } else if (nclass <= 7){
    cols <- c(pt[c(1, 3, 5)], "grey80", pt[c(7, 9, 10)])
  } else {
    warning("NOTE: please, provide up to 3 sections for stratification!")
    cols <- "grey80"
  }
  cols
}

##------------------------------------------------------------------------------
pal3 <- function(nclass){
  ptreds <- rev(colorRampPalette(brewer.pal(9, "Reds"))(11))
  ptblues <- rev(colorRampPalette(brewer.pal(9, "Blues"))(11))
  if (nclass == 1){
    cols <- "grey80"
  } else if (nclass <= 3){
    cols <- c(ptreds[c(4)], "grey80", rev(ptblues[c(5)]))
  } else if (nclass <= 5){
    cols <- c(ptreds[c(4, 7)], "grey80", rev(ptblues[c(4, 7)]))
  } else if (nclass <= 7){
    cols <- c(ptreds[c(2, 5, 8)], "grey80", rev(ptblues[c(2, 5, 8)]))
  } else {
    warning("NOTE: please, provide up to 3 sections for stratification!")
    cols <- "grey80"
  }
  cols
}

##------------------------------------------------------------------------------
.plotCox = function(coxtb, regs, keycovar, fpath, fname, width, height, 
                    xlim, xlab, ylab, plotpdf){
  #--- assign colors to line representation based on significance
  pal <- c("black", "#008080ff", "grey60", "#d45500ff")
  cols <- rep(NA, nrow(coxtb))
  names(cols) <- rownames(coxtb)
  cols[keycovar] <- 1
  idx1 <- coxtb[regs, "HR"] <= 1
  idx2 <- coxtb[regs, "Upper95"] < 1 | coxtb[regs, "Lower95"] > 1
  cols[regs][idx1 & idx2] <- 2
  cols[regs][!idx2] <- 3
  cols[regs][!idx1 & idx2] <- 4
  cols[] <- pal[cols]
  
  #--- get name labels for the graph
  labs <- rownames(coxtb)
  mxchar <- max(nchar(labs))
  len <- max(mxchar/10, 1)
  if(is.null(keycovar)){
    pal <- c("#008080ff", "grey60", "#d45500ff")
    leg1 <- c("associated, HR<1", "associated, HR>1","not associated")
    if(ylab=="Regulons and other covariates") ylab <- "Regulons"
  } else {
    pal <- c("black", "#008080ff", "grey60", "#d45500ff")
    leg1 <- c("other covariates", "associated, HR<1", "not associated","associated, HR>1")
  }
  
  #---
  xlabs <- .prettylog(xlim)
  xlim <- range(xlim)
  
  #--- plot graph
  nIN = 2
  nOUT = nrow(coxtb) + 1
  ylim <- c(0, nOUT + 1)
  if(plotpdf){
    fname = paste(fname, ".pdf", sep = "")
    pdf(file = paste(fpath, "/", fname, sep = ""), 
        width = width, height = height)
  }
  op <- par(no.readonly = TRUE)
  par(mai = c(0.4, 1.4 * len, 0.8, 0.7), mgp = c(3, 0.4, 0), yaxs = "i", xaxs = "i")
  plot(NA, log = "x", xlim = xlim, ylim = ylim, axes = FALSE, ylab = "", xlab = "")
  segments(xlim[1], nIN:nOUT, coxtb[, "Lower95"], col = "grey85", lwd = 1.5, lty = "21", lend = 2)
  lines(x = c(1, 1), y = c(1, nOUT + 1), lwd = 1.5, col = "grey60", lty = "21", lend = 2)
  segments(coxtb[, "Lower95"], nIN:nOUT, coxtb[, "Upper95"], col = cols, lwd = 1.5)
  points(y = nIN:nOUT, x = coxtb[, "HR"], pch = 18, cex = 1.2, lwd = 1, col = cols)
  axis(3, lwd = 2, cex.axis = 1.2, tck = -0.02, labels = xlabs$labs, at = xlabs$at)
  mtext(xlab, side = 3, line = 1.5, cex = 1.2)
  mtext(ylab, side = 2, line = 1 + 3 * len, cex = 1.2, adj = 0.6)
  labs <- rownames(coxtb)
  mtext(text = labs, side = 2, at = nIN:nOUT, las = 2, cex = 0.8)
  par(xpd = TRUE)
  if(is.null(keycovar)){
    legend(y = 0.5, x = 1, legend = leg1, col = pal[c(1,3,2)], pch = 18, lwd = 1.2, 
           xjust = 0.5, horiz=T, yjust = 0.5, bty = "n", cex = 0.65, pt.cex = 0.9, 
           x.intersp=0.5)
  } else {
    legend(y = 0, x = 1, legend = leg1, ncol = 2, col = pal, pch = 18, lwd = 1.2, 
           xjust = 0.5, yjust = 0.5, bty = "n", cex = 0.65, pt.cex = 0.9)
  }
  
  if (plotpdf){
    tp1 <- c("NOTE: file '",fname,"' should be available either in the working directory or\n")
    tp2 <- c("in a user's custom directory!\n")
    message(tp1,tp2)
    dev.off()
  }
  par(op)
}

##------------------------------------------------------------------------------
.prettylog <- function(xlim, n=10){
  # urd <- function(d, x){
  #   lxd <- log10(x/d)
  #   rlxd <- unique(c(floor(lxd), ceiling(lxd)))
  #   d * 10^rlxd
  # }
  xlim <- sort(xlim)
  if (length(xlim) == 2){
    tp <- seq(xlim[1], xlim[2], by = xlim[1])
    tp <- sort(unique(c(1, xlim, tp)))
    # xlabs <- sort(urd(1, tp))
    xlabs <- sort(c(1, xlim))
    #---
    tp <- unlist(sapply(1:(length(xlabs) - 1), function(i)
    {
      pretty(xlabs[i:(i + 1)], n = n)
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
  } else {
    xat <- xlabs <- xlim
  }
  idx <- xat >= 10
  tp <- c(format(xat[!idx], scientific = FALSE, digits=1), format(xat[idx], nsmall = 0))
  tp[!xat %in% xlabs] <- NA
  tp[xat == 1] <- "1.0"
  list(labs = tp, at = xat)
}
##------------------------------------------------------------------------------
tns.set <- function(object, para = NULL, what){
  if (what == "status"){
    object@status["Preprocess"] <- "[x]"
    object@status["Activity"] <- "[ ]"
    object@status["KM"] <- "[ ]"
    object@status["KmInt"] <- "[ ]"
    object@status["Cox"] <- "[ ]"
    object@status["CoxInt"] <- "[ ]"
  } else if (what == "regulonActivity"){
    object@results$regulonActivity <- para
    object@status["Activity"] <- "[x]"
  } else if (what == "KM"){
    object@results$KM <- para
    object@status["KM"] <- "[x]"
  } else if (what == "Cox"){
    object@results$Cox <- para
    object@status["Cox"] <- "[x]"
  } else if (what == "CoxInt"){
    object@results$CoxInt <- para
    object@status["CoxInt"] <- "[x]"
  } else if (what == "KmInt"){
    object@results$KmInt <- para
    object@status["KmInt"] <- "[x]"
  } else if (what == "survivalData"){
    object@survivalData <- para
  } else if (what == "para"){
    object@para <- para
  }
  return(object)
}

##------------------------------------------------------------------------
##returns rejection threshold for methods in 'p.adjust'
p.threshold <- function (pvals, alpha=0.05, method="BH"){
  pvals <- sort(pvals)
  padj <- p.adjust(pvals, method = method)
  thset <- which(padj <= alpha)
  if(length(thset)>0){
    mx1 <- mx2 <- which.max(thset)
    if(mx2<length(padj)) mx2 <- mx2 + 1
    th <- (pvals[mx1] + min(pvals[mx2],alpha) ) / 2
  } else {
    th <- min(c(alpha,pvals))
  }
  return(th)
}



