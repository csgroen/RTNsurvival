#' Subgroup Regulon Enrichment for TNS-class objects
#' 
#' This method evaluates which regulons are enriched in sample groups, given a
#' grouping variable. It performs Fisher's Exact Test whether a regulon is
#' positively or negatively enriched in a subgroup using regulon activity.
#' 
#' @param tns A A \linkS4class{TNS} object.
#' @param subgroup a character vector. It must be the name of a column in the
#' survivalData featuring the grouping information as a categorical variable.
#' @param regs An optional string vector specifying regulons to use for the analysis.
#' @param pValueCutoff a single numeric value specifying the cutoff for 
#' p-values considered significant.
#' @param pAdjustMethod a single character value specifying the p-value 
#' adjustment method to be used (see 'p.adjust' for details).
#' 
#' @return A TNS-class object with the results of the subgroup regulon enrichment
#' added to the results slot. To recover the results, use tnsGet(tns, "regulonEnrichment")
#' 
#' @examples 
#' # load survival data
#' data(survival.data)
#' # load TNI-object
#' data(stni, package = "RTN")
#' 
#' # create TNS object
#' stns <- tni2tnsPreprocess(stni, survivalData = survival.data,
#'                           keycovar = c('Grade','Age'), time = 1, event = 2)
#' stns <- tnsGSEA2(stns)
#' 
#' # run subgroup regulon enrichment analysis
#' stns <- tnsSRE(stns, "ER+")
#' 
#' # plot the result
#' tnsPlotSRE(stns)
#' @importFrom data.table rbindlist dcast
#' @importFrom stats fisher.test
#' @docType methods
#' @rdname tnsSRE-methods
#' @aliases tnsSRE
#' @export

setMethod("tnsSRE", "TNS", 
          function(tns, subgroup, regs = NULL, pValueCutoff = 0.05, 
                             pAdjustMethod = "BH") {
    #-- Basic check + get survData
    .tns.checks(tns, type = "Activity")
    survData <- tnsGet(tns, "survivalData")
    
    #-- Checks
    .tns.checks(subgroup, survData, type = "subgroup")
    .tns.checks(regs, type = "regs")
    .tns.checks(pValueCutoff, type = "pValueCutoff")
    .tns.checks(pAdjustMethod, type = "pAdjustMethod")
    
    #-- Get data
    regact <- tnsGet(tns, "regulonActivity")$dif
    all_regels <- tni.get(tnsGet(tns, "TNI"), "regulatoryElements")
    
    #-- regulatoryElements check
    if (!is.null(regs)) {
        if(all(regs %in% names(all_regels))) {
            regact <- regact[,regs] 
        } else {
            stop("All `regs` must be in the the regulatoryElements of the TNI-object.")
        }
    }
    
    grouping <- split(rownames(survData), survData[,subgroup])
    
    #-- Create Results table
    res_list <- lapply(names(grouping), function(group) {
        group_res <- rbindlist(
            lapply(colnames(regact), .regulonGroupFET, regact, group, grouping, pValueCutoff))
    })
    res_tb <- rbindlist(res_list)
    res_tb$FET_pAdjusted <- p.adjust(res_tb$FET_pValue, method = pAdjustMethod)
    
    #-- Use pAdjusted for Enrichment Test:
    res_tb[res_tb$FET_pAdjusted > pValueCutoff, "Enrichment_mode"] <- "0"
    
    #-- Reorder
    res_tb <- res_tb[order(res_tb$FET_pAdjusted),]
    
    #-- Add to tns
    tns@results$subgroupEnrichment[[subgroup]] <- res_tb
    tns@para$srFET$pValueCutoff[[subgroup]] <- pValueCutoff
    
    return(tns)
})

#' Plot Subgroup Regulon Enrichment for TNS-class objects
#' 
#' This method plots the results of the subgroup regulon enrichment analysis in
#' a heatmap. The rows of the heatmap represent enriched regulons, while the
#' columns show the subgroups. The plotted values correspond to average regulon
#' activity for a regulon in a subgroup. Enriched values can be marked.
#' 
#' @param tns A A \linkS4class{TNS} object.
#' @param subgroup a character vector. It must be the name of a column in the
#' survivalData featuring the grouping information as a categorical variable.
#' @param by one of 'nGroups' or 'groupTop'. If by = 'nGroups', the nGroupsEnriched
#' value will be used to select regulons. If by = 'groupTop', 'nTopEnriched' will
#' be used to select regulons for plotting.
#' @param nGroupsEnriched a single integer. It represents in how many subgroups
#' a regulon has to be enriched for it to appear in the rows of the heatmap.
#' @param nTopEnriched a single integer. If by = 'groupTop', this represents how
#' regulons will be shown for each group (duplicates are removed. The top regulons
#' are chosen by significance.
#' @param breaks a numerical vector of breaks for the heatmap. 
#' @param markEnriched a single logical value. If TRUE, asterisks are added to
#' cells of heatmap that were found to be significant by tnsSRE.
#' @param ... parameters passed to pheatmap::pheatmap for customization.
#' 
#' @return A heatmap of the subgroup regulon enrichment results.
#' 
#' @examples 
#' # load survival data
#' data(survival.data)
#' # load TNI-object
#' data(stni, package = "RTN")
#' 
#' # create TNS object
#' stns <- tni2tnsPreprocess(stni, survivalData = survival.data,
#'                           keycovar = c('Grade','Age'), time = 1, event = 2)
#' stns <- tnsGSEA2(stns)
#' 
#' # run subgroup regulon enrichment analysis
#' stns <- tnsSRE(stns, "ER+")
#' tnsPlotSRE(stns)
#' 
#' @importFrom data.table rbindlist dcast
#' @importFrom pheatmap pheatmap
#' @docType methods
#' @rdname tnsPlotSRE-methods
#' @aliases tnsPlotSRE
#' @export
setMethod("tnsPlotSRE", "TNS", 
          function(tns, subgroup = NULL, by = "nGroups",
                       nGroupsEnriched = 1, nTopEnriched = 10, 
                       breaks = seq(-1.5, 1.5, 0.1),
                       markEnriched = FALSE, ...) {
    
    #-- Basic check + get survData
    .tns.checks(tns, type = "Activity")
    survData <- tnsGet(tns, "survivalData")
    
    #-- Checks
    .tns.checks(nGroupsEnriched, type = "nGroupsEnriched")
    .tns.checks(nTopEnriched, type = "nTopEnriched")
    .tns.checks(breaks, type = "breaks")
    .tns.checks(markEnriched, type = "markEnriched")
    
    
    #-- Subgroup check
    if (is.null(subgroup)) {
        if (length(tns@results$subgroupEnrichment) == 1) {
            fet_res <- tns@results$subgroupEnrichment[[1]]
        } else if (length(tns@results$subgroupEnrichment) > 1) {
            stop("`subgroup` must be specified for plotting")
        } else {
            stop("`tnsSRE` must be run before `tnsPlotSRE`")
        }
    } else {
        if (subgroup %in% names(tns@results$subgroupEnrichment)) {
            fet_res <- tns@results$subgroupEnrichment[[subgroup]]
        } else {
            stop("`tnsSRE` must be run with `subgroup` before plotting.")
        }
    }
    
    fet_ls <- split(fet_res, fet_res$Regulon)
    
    #-- by check
    if(by == "nGroups") {
        #-- Filter to keep regulons enriched in nGroupsEnriched
        idx <- sapply(fet_ls, function(regtb) { 
            sum(regtb$Enrichment_mode != 0) >= nGroupsEnriched
        })
    } else if (by == "groupTop") {
        fet_Gls <- split(fet_res, fet_res$Group)
        
        #-- Filter to keep top n regulons in each group
        idx <- sapply(fet_Gls, function(grouptb) {
            unlist((grouptb[order(grouptb$FET_pValue), "Regulon"][1:nTopEnriched]))
        })
        idx <- unique(as.vector(idx))
    } else {
        stop("`by` is an invalid option")
    }
    
    filt_fet_res <- rbindlist(fet_ls[idx])
    
    #-- Spread data to a dES_average matrix
    cast_res <- as.data.frame(
        dcast(filt_fet_res, Regulon ~ Group, value.var = "dESaverage"))
    rownames(cast_res) <- cast_res$Regulon
    cast_res <- cast_res[,-1]
    
    colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(length(breaks) - 1)
    
    
    if(markEnriched) {
        #-- Spread p-values for "*"s
        cast_ast <- as.data.frame(dcast(filt_fet_res, 
                                        Regulon ~ Group, 
                                        value.var = "FET_pValue"))
        rownames(cast_ast) <- cast_res$Regulon
        cast_ast <- cast_ast[,-1]
        
        #-- Get analysis pval
        if (is.null(subgroup))
            pval <- tnsGet(tns, "para")$srFET$pValueCutoff
        else
            pval <- tnsGet(tns, "para")$srFET[subgroup]
        
        #-- Transform p < pval into "*", otherwise ""
        idx <- cast_ast < pval
        cast_ast[idx] <- "*"
        cast_ast[!idx] <- ""
        
        #-- Plot
        pheatmap(cast_res,
                 color = colors,
                 breaks = breaks,
                 border_color = NA, 
                 display_numbers = cast_ast,
                 ...)
    } else {
        colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(length(breaks) - 1)
        pheatmap(cast_res,
                 color = colors,
                 breaks = breaks,
                 border_color = NA,
                 ...)
    }
    
    invisible(cast_res)
})

.regulonGroupFET <- function(reg, regact, group, grouping, pValueCutoff) {
    #-- Getting group position in grouping list
    grn <- which(names(grouping) == group)
    
    #-- Initializing result
    res <- list(
        Regulon = reg,
        Group = group,
        Enrichment_mode = 0,
        FET_pValue = NA,
        dESaverage = 0,
        lowerCI = NA
        )
    
    #-- For one regulon
    reg_dES <- regact[,reg]
    idx <- reg_dES >= 0
    act <- names(reg_dES[idx])
    rep <- names(reg_dES[!idx])
    
    #-- Positively enriched test
    length(act)
    ct <- matrix(c(
        length(intersect(act, grouping[[grn]])),
        length(intersect(rep, grouping[[grn]])),
        length(intersect(act, unlist(grouping[-grn]))),
        length(intersect(rep, unlist(grouping[-grn])))
    ), nrow = 2, dimnames = list(Regulon = c("Active", "Repressed"),
                                 Group = c("Belong", "Don't belong")))
    
    ft_act <- fisher.test(ct, alternative = "greater")
    
    #-- Negatively enriched test
    ft_rep <- fisher.test(ct[2:1,], alternative = "greater")
    
    #-- Add results
    ft <- list(ft_act, ft_rep)
    idx <- which.min(c(ft_act$p.value, ft_rep$p.value))
    res$FET_pValue <- ft[[idx]]$p.value
    res$lowerCI <- ft[[idx]]$conf.int[1]
    res$Enrichment_mode <- ifelse(idx == 1, "1", "-1")
    res$Enrichment_mode <- ifelse(res$FET_pValue > pValueCutoff, 0, res$Enrichment_mode)
    res$dESaverage <- mean(reg_dES[grouping[[grn]]])
    
    return(res)
}
