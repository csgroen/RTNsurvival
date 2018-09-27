#' Plot regulon activity and categorical covariates
#' 
#' This method plots regulon activity for a given regulon in all samples and 
#' adds covariate tracks to evaluate the regulon activity distribution. The
#' samples are order by regulon activity for that particular regulon.
#' 
#' Automatic dummy encoding is available for categorical variables.
#' 
#' @param tns A A \linkS4class{TNS} object.
#' @param regs An optional string vector specifying regulons to plot.
#' @param attribs A character vector of column names from the survivalData. All
#' attribs should be either binary encoded or categorical variables for plotting.
#' Available attribs can be checked by running colnames(tnsGet(tns, "survivalData"))
#' @param fname A string. The name of the file in which the plot will be saved
#' @param fpath A string. The path to the directory where the plot will be saved
#' @param plotpdf A logical value. If TRUE, a pdf file is created instead of 
#' plotting to the graphics device.
#' @param plotbatch A logical value. If TRUE, plots for all regulons are saved 
#' in the same file. If FALSE, each plot for each regulon is saved in a different file.
#' @param panelHeights A numeric vector of length 2 specifying the relative heights
#' of the panels (regulon activity plot, and covariate tracks)
#' @param width A numeric value. Represents the width of the plot.
#' @param height A numeric value. Represents the height of the plot.
#' @param dummyEncode A logical value. If TRUE, all categorical variables are
#' dummy encoded. If FALSE, categorical variables are represented as one track 
#' and a legend is added to the plot.
#' @param divs A numeric vector of division positions in the covariate tracks.
#' 
#' @return A plot of regulon activity and covariate tracks.
#' @examples 
#'  # load survival data
#'  data(survival.data)
#'  # load TNI-object
#'  data(stni, package = "RTN")
#'  
#'  # create TNS object
#'  stns <- tni2tnsPreprocess(stni, survivalData = survival.data,
#'  keycovar = c('Grade','Age'), time = 1, event = 2)
#'  stns <- tnsGSEA2(stns)
#'  
#'  # plot only binary covariates
#'  tnsPlotCovariates(stns, "MYB", 
#'  attribs = c("ER+", "ER-", "PR+", "PR-", "LumA", "LumB", "Basal", "Her2",
#'  "Normal"),
#'  divs = c(2, 4))
#'  
#'  # also dummy encode categorical variables (LN and Grade)
#'  tnsPlotCovariates(stns, "MYB", 
#'  attribs = c("ER+", "ER-", "PR+", "PR-", "LumA", "LumB", "Basal", "Her2",
#'  "Normal", "LN", "Grade"),
#'  divs = c(2, 4, 9, 12))
#'  
#'  # don't dummy encode categorical variables
#'  tnsPlotCovariates(stns, "MYB", attribs = c("ER+", "ER-", "PR+", "PR-",
#'  "LumA", "LumB", "Basal", "Her2", "Normal", "Grade"), divs = c(2, 4, 9),
#'  dummyEncode = FALSE)
#'  
#' @importFrom RColorBrewer brewer.pal
#' @importFrom data.table melt
#' @importFrom stats na.exclude na.pass model.matrix
#' @import ggplot2
#' @importFrom egg ggarrange
#' @docType methods
#' @rdname tnsPlotCovariates-methods
#' @aliases tnsPlotCovariates
#' @export
setMethod("tnsPlotCovariates", "TNS", 
          function(tns, regs = NULL, attribs = NULL, fname = "covarplot", 
                   fpath = ".", plotpdf = FALSE, plotbatch = FALSE,
                   panelHeights = c(1,1), width = 5.3, height = 4,
                   dummyEncode = TRUE, divs = NULL) {
    #-- Parameter checks
    .tns.checks(tns, type = "Activity")
    .tns.checks(regs, type = "regs")
    .tns.checks(attribs, tns@survivalData, type = "attribs2")
    .tns.checks(fname, type = "fname")
    .tns.checks(fpath, type = "fpath")
    .tns.checks(plotpdf, type = "plotpdf")
    .tns.checks(plotbatch, type = "plotbatch")
    .tns.checks(width, type = "width")
    .tns.checks(height, type = "height")
    .tns.checks(panelHeights, type = "panelHeights")
    .tns.checks(dummyEncode, type = "dummyEncode")
    .tns.checks(divs, attribs, type = "divs")
    tnstatus <- tnsGet(tns, what = "status")
    if(tnstatus["Activity"] != "[x]")
        stop("NOTE: TNS object needs to be evaluated by 'tnsGSEA2' or 'tnsAREA3'!",
             call. = FALSE)

    #-- Get data
    regact <- tnsGet(tns, "regulonActivity")$dif
    regstatus <- tnsGet(tns, "regulonActivity")$regstatus
    regstatus <- apply(regstatus, 1:2, as.character)
    survData <- tnsGet(tns, "survivalData")
    
    #-- Get regs
    if (is.null(regs)) {
        regs <- colnames(regact)
    }
    
    #-- Check dummyEncode
    if(!is.logical(dummyEncode) && !(dummyEncode %in% colnames(survData))) {
        stop("`dummyEncode` must be either a logical value or a character vector of names of columns to dummy encode.")
    }
    
    #-- Create plotData
    plotData <- data.frame(rownames(regact), regact[,regs], regstatus[,regs])
    colnames(plotData) <- c("Sample_name", regs, paste0(regs, "_status"))
    
    #-- attribs preprocess
    if (is.null(attribs)) {
        all_attribs <- plotData[,0]
    } else {
        covars <- survData[,attribs]
        if (isFALSE(dummyEncode)){
            all_attribs <- covars
            
            #-- Check which covars are not binary
            idx <- apply(covars, 2, function(covar) {
                !all(covar %in% c(1, 0))
            })
            
            #-- Treat exceptions (class = 0 or 1)
            non_dummy_attribs <- as.data.frame(all_attribs[,idx])
            non_dummy_attrib_names <- colnames(all_attribs)[idx]
            fixed_ndattribs <- sapply(1:sum(idx), function(i) {
                col <- non_dummy_attribs[,i]
                attrib_name <- non_dummy_attrib_names[i]
                if (any(as.character(col) %in% c("0", "1"))) {
                    col <- paste(attrib_name, col, sep = "_")
                }
                return(col)
            })
            all_attribs <- cbind(all_attribs[,!idx], fixed_ndattribs)
            colnames(all_attribs) <- c(colnames(all_attribs)[!idx], non_dummy_attrib_names)
            
            } else {
            #-- Check which covars are binary
            idx <- apply(covars, 2, function(covar) {
                !all(covar %in% c(1, 0))
            })
            bincovars <- covars[,!idx]
            
            #-- Get divs
            if(is.null(divs) && all(idx)) {
                nlevels <- apply(covars, 2, function(col) { length(unique(col)) })
                divs <- cumsum(nlevels)[1:length(nlevels)-1]
            }
            
            if (isTRUE(dummyEncode)) {
                #-- Make dummy variables for non-binary covars
                encoded_covars_ls <- lapply(names(idx)[idx], dummyEncodeCovar, covars)
                encoded_covars <- do.call(cbind, encoded_covars_ls)
            } else if (is.character(dummyEncode)) {
                nonDE <- names(idx)[!(names(idx)[idx] %in% dummyEncode)]
                encoded_covars_ls <- lapply(dummyEncode, dummyEncodeCovar, covars)
                encoded_covars <- do.call(cbind, encoded_covars_ls)
                encoded_covars <- cbind(encoded_covars, covars[,nonDE])
                colnames(encoded_covars)[colnames(encoded_covars) == "covars[, nonDE]"] <- nonDE
            }
            
            #-- Add binary and encoded non-binary covars to plotData
            if(is.null(encoded_covars)) {
                all_attribs <- bincovars[rownames(plotData),]
            } else {
                all_attribs <- cbind(encoded_covars[rownames(plotData),], bincovars[rownames(plotData),])
            }
        }
    }
    #-- Fix table order
    og_order <- unlist(sapply(attribs, grep, colnames(all_attribs), fixed = TRUE))
    all_attribs <- all_attribs[,og_order]
    
    #-- Add to plot data
    plotData  <- cbind(plotData, all_attribs)
    attrib_names <- colnames(all_attribs)
    
    #-- Plot covars
    allPlots <- lapply(regs, ggPlotCovariates, plotData, 
                                        attrib_names, panelHeights, dummyEncode, divs)
    
    if (plotpdf) {
        #-- Treat fname
        if (fname == "covarplot") {
            if (length(regs) == 1) {
                fname <- paste(fname, regs, sep = "_")
            } else if (plotbatch) {
                fname <- paste(fname, "regs", sep = "_")
            }
        } else {
            fname <- gsub(".pdf", '',fname, ignore.case = TRUE)
        }
        
        #-- Plot
        if (length(regs) == 1) { #-- Plot pdf one reg
            ggsave(filename = paste0(fpath, "/", fname, ".pdf"),
                   plot = allPlots[[1]]$grid_plot, height = height, width = width)
        } else if (plotbatch) { #-- Plot batch multiple regs
            pdf(file = paste0(fpath, "/", fname, ".pdf"), 
                width = width, height = height)
            lapply(allPlots, function(plots) {
                print(plots[["grid_plot"]])
            })
            dev.off()
        } else { #-- Plot each reg in a pdf
            for (i in 1:length(regs)) {
                new_fname <-  paste(fname, regs[i], sep = "_")
                ggsave(filename = paste0(fpath, "/", new_fname, ".pdf"),
                       plot = allPlots[[i]]$grid_plot, height = height, width = width)
            }
        }
        msg <- c("NOTE: file '",fname,"' should be available either in the working directory or in a user's custom directory!\n")
        message(msg)
        
    } else { #-- Plot to the graphics device
        for (i in 1:length(allPlots)) {
            print(allPlots[[i]]$grid_plot)
        }
    }
    
    #-- Return
    plot_list <- lapply(allPlots, "[[", 2)
    names(plot_list) <- regs
    return(invisible(plot_list))
})


ggPlotCovariates <- function(reg, plotData, attrib_names, panelHeights, dummyEncode,
                             divs) {
    #-- Copy data
    plotData_reg <- plotData
    
    #-- Change `reg` name (for plotting function)
    colnames(plotData_reg)[colnames(plotData_reg) == reg] <- "reg"
    colnames(plotData_reg)[colnames(plotData_reg) == paste0(reg, "_status")] <- "reg_status"
    
    #-- Reorder samples
    plotData_reg <- plotData_reg[order(plotData_reg$reg),]
    plotData_reg$Samples <- 1:nrow(plotData_reg)
    
    #-- Get colors
    n <- length(unique(plotData_reg$reg_status))
    pal <- pal3(n)
    
    #-- First plot, (regulon activity + stratification)
    if(length(attrib_names) == 0) {
        p1 <- ggDesPlot(plotData_reg, pal, reg, xaxis = "bottom")
        plot <- list(grid_plot = p1, ggplots = p1)
        return(plot)
    }
    p1 <- ggDesPlot(plotData_reg, pal, reg)
    
    #-- Melt data for second plot
    attribData <- as.data.frame(apply(plotData_reg[,attrib_names], 1:2, "as.character"))
    attribData$Samples <- plotData_reg$Samples
    plotData_melt <- suppressWarnings(melt(attribData, id.vars = "Samples",
                                      measure.vars = attrib_names,
                                      variable.name = "Covariates"))
    
    #-- Second plot (covariate tracks)
    p2 <- ggCovariateTracks(plotData_melt, dummyEncode, divs, attrib_names)
    
    #-- Align
    grid_plot <- ggarrange(p1, p2, nrow = 2, 
                           heights = panelHeights, 
                           draw = FALSE,
                           padding = unit(3, "line"))
    ggplots <- list(p1, p2)
    
    return(list(grid_plot = grid_plot, ggplots = ggplots))
}

ggDesPlot <- function(plotData_reg, pal, reg, xaxis = "top", flipPlot = FALSE) {
    p <- ggplot(plotData_reg, aes_string("Samples", "reg")) +
        geom_bar(aes_string(fill = "reg_status"), stat = "identity", width = 1) +
        annotate("text", x = 0, y = 1.7, label = reg, hjust = -0.2) +
        scale_fill_manual(values = pal) +
        scale_y_continuous(name = "Regulon activity (dES)", 
                           limits = c(-2, 2), expand = c(0,0)) +
        guides(fill = FALSE) +
        theme_classic() +
        theme(plot.margin = unit(c(2,4,2,4), "mm"))
    if (xaxis == "top") {
        p <- p + scale_x_continuous(limits = c(0, nrow(plotData_reg)), 
                           expand = c(0,0), position = "top")
    } else {
        p <- p + scale_x_continuous(limits = c(0, nrow(plotData_reg)), 
                                    expand = c(0,0))
    }
    if (flipPlot) {
        p <- p + coord_flip()
    }
    return(p)
    
} 

ggCovariateTracks <- function(plotData_melt, dummyEncode, divs, attrib_names) {
    #-- Get colors
    allattr <- na.exclude(unique(plotData_melt$value))
    if (isFALSE(dummyEncode) && !any(grepl("0", allattr))) {
        gbpal <- NULL
        nattr <- length(allattr)
    } else {
        gbpal <- c("0" = "grey95", "1" = "black")
        nattr <- length(allattr) - 2
    }
    attrnames <- allattr[!(allattr %in% c("0", "1"))]
    
    if (nattr <= 12 && nattr > 0) {
        addcols <- brewer.pal(nattr, "Set3")
        names(addcols) <- attrnames
    } else if (nattr > 0) {
        addcols <- colorRampPalette(brewer.pal(12, "Set3"))(nattr)
        names(addcols) <- attrnames
    } else {
        addcols <- NULL
    }
    fullpal <- c(gbpal, addcols)
    plotData_melt$Covariates <- factor(plotData_melt$Covariates,
                                       levels = rev(levels(plotData_melt$Covariates)))
    
    p2 <- ggplot(plotData_melt) +
        geom_tile(aes_string("Samples", "Covariates", fill = "value")) +
        scale_fill_manual(values = fullpal, breaks = attrnames) +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        theme_classic() +
        theme(axis.title = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.y = element_blank(),
              axis.line.x = element_blank(),
              legend.title = element_blank())
    
    if (!is.null(divs)) {
        p2 <- p2 + 
            geom_hline(yintercept = length(attrib_names)-divs+0.5, 
                       color = "white", size = 3/length(attrib_names)*5)
    }
    return (p2)
}

dummyEncodeCovar <- function (covar_name, covars_tb) {
    covar <- as.factor(covars_tb[,covar_name])
    
    #-- treat exception: NA in categorical data
    og.op <- options()
    options(na.action=na.pass)
    
    #-- dummy encode
    encoded_covar <- as.data.frame(model.matrix(~ 0 + covar),
                                   row.names = rownames(covars_tb))
    options(og.op)
    
    #-- fix names
    level_names <- gsub("^covar", "", colnames(encoded_covar))
    colnames(encoded_covar) <- paste(covar_name, level_names, sep = "_")
    
    return(encoded_covar)
}

pal3 <- function(nclass){
    ptreds <- rev(colorRampPalette(brewer.pal(9, "Reds"))(11))
    ptblues <- rev(colorRampPalette(brewer.pal(9, "Blues"))(11))
    if (nclass == 1){
        cols <- "grey80"
    } else if (nclass <= 3){
        cols <- c(ptreds[c(4)], "grey80", rev(ptblues[c(5)]))
        names(cols) <- as.character(seq(1, -1, -1))
        cols
    } else if (nclass <= 5){
        cols <- c(ptreds[c(4, 7)], "grey80", rev(ptblues[c(4, 7)]))
        names(cols) <- as.character(seq(2, -2, -1))
        cols
    } else if (nclass <= 7){
        cols <- c(ptreds[c(2, 5, 8)], "grey80", rev(ptblues[c(2, 5, 8)]))
        names(cols) <- as.character(seq(3, -3, -1))
        cols
    } else {
        warning("NOTE: please, provide up to 3 sections for stratification!")
        cols <- "grey80"
    }
    cols
}
