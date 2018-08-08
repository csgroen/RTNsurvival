
##------------------------------------------------------------------------------
.survstatsDuals <- function(regulonActivity, survData, regs, excludeMid)
{
  
  #--- survData
  survData <- survData[,c("time","event")]
  survData <- survData[complete.cases(survData),]
  
  #--- tabstatus
  tabstatus <- regulonActivity$regstatus[rownames(survData), regs]
  tabstatus <- data.frame(tabstatus)
  idx <- rowSums(tabstatus==0)>0
  if(excludeMid){
    tabstatus <- tabstatus[!idx,]
    sections <- 1:4
    names(sections) <- c("-|-","-|+","+|-","+|+")
  } else {
    tabstatus[idx,] <- 0
    sections <- 1:5
    names(sections) <- c("0|0","-|-","-|+","+|-","+|+")
  }
  survData <- survData[rownames(tabstatus),]
  
  #--- regstatusChar
  tp1 <- as.character(tabstatus[[regs[1]]])
  tp1[tp1=="1"]  <- "+"
  tp1[tp1=="-1"] <- "-"
  tp2 <- as.character(tabstatus[[regs[2]]])
  tp2[tp2=="1"] <- "+"
  tp2[tp2=="-1"] <- "-"
  regstatusChar <- paste(tp1,tp2, sep = "|")
  names(regstatusChar) <- rownames(tabstatus)
  
  #--- regstatusNum
  regstatusNum <- sections[regstatusChar]
  names(regstatusNum) <- names(regstatusChar)
  nclass <- length(unique(regstatusNum))
  ddt <- survData[names(regstatusNum), ]
  ddt$class <- regstatusNum
  
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
.survplotDuals <- function(regulonActivity, survData, regs, endpoint,
                           excludeMid, ylab, xlab, colorPalette)
{
  
  #--- survData
  survData <- survData[,c("time","event")]
  survData <- survData[complete.cases(survData),]
  
  #--- tabdiff
  tabdiff <- regulonActivity$dif[rownames(survData), regs]
  survData <- survData[rownames(tabdiff),]
  
  #--- tabstatus
  tabstatus <- regulonActivity$regstatus[rownames(survData), regs]
  tabstatus <- data.frame(tabstatus)
  idx <- rowSums(tabstatus==0)>0
  if(excludeMid){
    tabstatus <- tabstatus[!idx,]
    sections <- 1:4
    names(sections) <- c("-|-","-|+","+|-","+|+")
  } else {
    tabstatus[idx,] <- 0
    sections <- 1:5
    names(sections) <- c("0|0","-|-","-|+","+|-","+|+")
  }
  survData <- survData[rownames(tabstatus),]
  
  #--- regstatusChar
  tp1 <- as.character(tabstatus[[regs[1]]])
  tp1[tp1=="1"]  <- "+"
  tp1[tp1=="-1"] <- "-"
  tp2 <- as.character(tabstatus[[regs[2]]])
  tp2[tp2=="1"] <- "+"
  tp2[tp2=="-1"] <- "-"
  regstatusChar <- paste(tp1,tp2, sep = "|")
  names(regstatusChar) <- rownames(tabstatus)
  
  #--- regstatusNum
  regstatusNum <- sections[regstatusChar]
  names(regstatusNum) <- names(regstatusChar)
  
  #-- get colors
  if (is.singleString(colorPalette)){
    if (colorPalette == "red"){
      cols <- pal1(4)
    } else if (colorPalette == "blue"){
      cols <- pal2(4)
    } else if (colorPalette %in% c("redblue","bluered")){
      cols <- pal3(4)
    }
    if(colorPalette!="redblue") 
      cols <- rev(cols)
  } else {
    cols <- colorPalette
  }
  if(excludeMid){
    cols <- cols[-3]
  } else {
    cols <- cols[c(3,1,2,4,5)]
  }
  
  #--- adjusting graphical parameters
  op <- par(no.readonly=TRUE)
  par(mar = c(4, 4, 2, 2), mgp = c(2, 0.4, 0), cex=1)
  #-- survival analysis
  ddt <- survData[names(regstatusNum), ]
  ddt$class <- regstatusNum
  res1 <- survfit(Surv(time, event) ~ class, data = ddt)
  plot(res1, col = cols, lwd = 1.8, axes = FALSE, cex = 0.5, mark.time = TRUE, ylab = "", xlab = "")
  title(ylab = ylab, adj = 0.5, cex.lab = 1.2, mgp = c(2.2, 0.4, 0))
  title(xlab = xlab, adj = 0.5, cex.lab = 1.2, mgp = c(1.6, 0.4, 0))
  labs <- as.integer(seq(0, endpoint, length.out = 4))
  if (!endpoint %in% labs) labs <- pretty(c(0, endpoint))
  axis(1, at = labs, labels = labs, tcl = -0.2, las = 1, lwd = 1.8, cex.axis = 1.2)
  axis(2, tcl = -0.2, las = 2, lwd = 1.8, cex.axis = 1.2)
  #---log-rank test
  lrstats <- c(chisq=NA, p=NA)
  res2 <- survdiff(Surv(time, event) ~ class, data = ddt)
  pval <- 1 - pchisq(res2$chisq, length(res2$n) - 1)
  lrstats[] <- c(res2$chisq,pval)
  #---legends
  par(xpd=TRUE)
  legs <- names(sections); legs[1] <- "undetermined"
  legs <- paste(legs, " : ", res2$n, "(", res2$obs,")", sep = "")
  legend("bottomleft", legend = rev(legs), col = rev(cols), bty = "n", pch = 15, title.adj = 0, adj = 0,
         title = paste(paste(regs, collapse = " | "), "\nregulon status",sep=""), inset = c(0.01,0),
         cex = 0.8, pt.cex = 1)
  pval <- paste("Logrank P: ", format(pval, digits = 3, scientific = TRUE))
  legend("topright", cex = 0.8, legend = pval, bty = "n", inset = c(0,-0.05))
  par(op)
  return(lrstats)
}
.namesCorrect <- function(regs) {
  xregs <- gsub("-|\\+|\\.|:|\\*|,|;", "_", regs)
  xregs <- gsub("\\s", "", xregs)
  xregs
}

##------------------------------------------------------------------------------
.getSurvplotCols <- function(regulonActivity, regs, excludeMid, colorPalette){
  
  #--- tabdiff
  tabdiff <- regulonActivity$dif[, regs]
  
  #--- tabstatus
  tabstatus <- regulonActivity$regstatus[, regs]
  tabstatus <- data.frame(tabstatus)
  idx <- rowSums(tabstatus==0)>0
  if(excludeMid){
    tabstatus <- tabstatus[!idx,]
    sections <- 1:4
    names(sections) <- c("-|-","-|+","+|-","+|+")
  } else {
    tabstatus[idx,] <- 0
    sections <- 1:5
    names(sections) <- c("0|0","-|-","-|+","+|-","+|+")
  }
  
  #--- regstatusChar
  tp1 <- as.character(tabstatus[[regs[1]]])
  tp1[tp1=="1"]  <- "+"
  tp1[tp1=="-1"] <- "-"
  tp2 <- as.character(tabstatus[[regs[2]]])
  tp2[tp2=="1"] <- "+"
  tp2[tp2=="-1"] <- "-"
  regstatusChar <- paste(tp1,tp2, sep = "|")
  names(regstatusChar) <- rownames(tabstatus)
  
  #--- regstatusNum
  regstatusNum <- sections[regstatusChar]
  names(regstatusNum) <- names(regstatusChar)
  
  #-- get colors
  if (is.singleString(colorPalette)){
    if (colorPalette == "red"){
      cols <- pal1(4)
    } else if (colorPalette == "blue"){
      cols <- pal2(4)
    } else if (colorPalette %in% c("redblue","bluered")){
      cols <- pal3(4)
    }
    if(colorPalette!="redblue") 
      cols <- rev(cols)
  } else {
    cols <- colorPalette
  }
  if(excludeMid){
    cols <- cols[-3]
  } else {
    cols <- cols[c(3,1,2,4,5)]
  }
  regstatusCol <- cols[regstatusNum]
  names(regstatusCol) <- names(regstatusNum)
  res <- list(numb=regstatusNum, char=regstatusChar, 
              cols=regstatusCol)
  return(res)
}

