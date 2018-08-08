## This function is used for argument checking
.tns.checks <- function(object1, object2 = NULL, type)
{
  if (type == "survivalData"){
    if (!is.data.frame(object1)) 
      stop("'survivalData' must be a data frame.", call. = FALSE)
    if (is.null(colnames(object1))) 
      stop("columns in 'survivalData' must be named.", call. = FALSE)
  } 
  else if (type == "time"){
    if (!is.singleInteger(object1) && !is.singleString(object1)) 
      stop("'time' should be either a character or integer value !\n", 
           call. = FALSE)
    if (is.character(object1)){
      if (!object1 %in% colnames(object2)) 
        stop("'time' parameter must be present in 'survivalData' 
             colnames.", 
             call. = FALSE)
      object1 <- which(colnames(object2) == object1)
    }
    vals <- object2[, object1]
    if (!is.numeric(vals)){
      stop("'time' data must be numeric.", call. = FALSE)
    } else if (any(vals < 0, na.rm = TRUE)){
      stop("values in 'time' data must be >= 0", call. = FALSE)
    }
    return(object1)
  } 
  else if (type == "event"){
    if (!is.singleInteger(object1) && !is.singleString(object1)) 
      stop("'event' should be either a character or integer 
           value !\n", 
           call. = FALSE)
    if (is.character(object1)){
      if (!object1 %in% colnames(object2)) 
        stop("'event' parameter must be present in 'survivalData' 
             colnames.", 
             call. = FALSE)
      object1 <- which(colnames(object2) == object1)
    }
    vals <- object2[, object1]
    if (!is.numeric(vals)){
      stop("'event' data must be numeric.", call. = FALSE)
    } else if (!all.binaryValues(vals)){
      stop("'event' data must be either binary or logical.", call. = FALSE)
    }
    return(object1)
  } 
  else if (type == "Keycovars"){
    if(!is.null(object1)) {
      if (!is.character(object1)) 
        stop("'keycovar' must be a character vector.", call. = FALSE)
      if (!all(object1 %in% colnames(object2))) 
        stop("All strings in 'keycovar' must be colnames in 'survivalData'", 
             call. = FALSE)
      for (col in object1){
        if (!is.numeric(object2[, col])) 
          stop("All values in 'keycovar' columns must be numeric.", 
               call. = FALSE)
      }
    }
  } 
  else if (type == "samples"){
    if (is.null(object1)){
      object1 <- rownames(object2)
    } else {
      if (!is.character(object1)) 
        stop("'samples' must be a character vector.", call. = FALSE)
      if (!all(object1 %in% rownames(object2))) 
        stop("All strings in 'samples' must appear in the rownames of 'survivalData'", 
             call. = FALSE)
    }
    return(object1)
  } 
  else if (type == "survival_cox"){
    if (nrow(object1) < 50)
      warning("If the number of samples in 'survivalData' is too small, 
              coxph function may not converge.")
  } 
  else if (type == "fpath"){
    if (!is.singleString(object1)) 
      stop("'fpath' must be a single character.", call. = FALSE)
    if (!dir.exists(object1)) 
      stop("'fpath' does not lead to an existing directory.", call. = FALSE)
  } 
  else if (type == "fname"){
    if (!is.singleString(object1)) 
      stop("'fname' must be a single character.", call. = FALSE)
    #---check name
    validname <- gsub("[^0-9A-Za-z\\.]", '_',object1)
    if(validname!=object1){
      stop("NOTE: please provide 'fname' without special charaters or path information!",
           call. = FALSE)
    }
  } 
  else if (type == "ylab"){
    if (!is.singleString(object1)) 
      stop("'ylab' must be a single character.", call. = FALSE)
  } 
  else if (type == "xlab"){
    if (!is.singleString(object1)) 
      stop("'xlab' must be a single character.", call. = FALSE)
  } 
  else if (type == "regs"){
    if (!is.null(object1)){
      if (!all.characterValues(object1)) 
        stop("'regs' must be a character vector.", call. = FALSE)
    }
  } 
  else if(type == "regulatoryElements"){
    if(!all.characterValues(object1) || any(duplicated(object1)) ){
      stop("NOTE: 'regulatoryElements' should be unique character values !", 
           call. = FALSE)
    }
  } 
  else if (type == "attribs"){
    if (!is.null(object1)){
      if (is.list(object1)) 
        object1 <- unlist(object1)
      if (is.character(object1)){
        if (!all(object1 %in% colnames(object2))) 
          stop("all 'attribs' must be listed in the 'survivalData' colnames, at the 'tns' object!", 
               call. = FALSE)
      } else if (!all.integerValues(object1)){
        stop("'attribs' must be either a vector or a list of vectors, with either integer or character values!", 
             call. = FALSE)
      }
    }
  } 
  else if (type == "colorPalette"){
    len <- (object2 * 2) + 1
    tp1 <- "'colorPalette' must be 'red', 'blue', 'redblue' or 'bluered'"
    message <- paste(tp1,", or a vector with ", len," valid colors.", sep="")
    if(is.singleString(object1)){
      if (!object1 %in% c("red", "blue", "redblue","bluered"))
        stop(message, call. = FALSE)
    } else if(!is.color(object1) || length(object1)!=len){
      stop(message, call. = FALSE)
    }
  }
  else if (type == "excludeMid"){
    if (!is.singleLogical(object1)) 
      stop("'excludeMid' must be a logical value.", call. = FALSE)
  }
  else if (type == "plotpdf"){
    if (!is.singleLogical(object1)) 
      stop("'plotpdf' must be a logical value.", call. = FALSE)
  } 
  else if (type == "showdata"){
    if (!is.singleLogical(object1)) 
      stop("'showdata' must be a logical value.", call. = FALSE)
  } 
  else if (type == "plotbatch"){
    if (!is.singleLogical(object1)) 
      stop("'plotbatch' must be a logical value.", call. = FALSE)
  } 
  else if (type == "qqkeycovar"){
    if (!is.singleLogical(object1)) 
      stop("'qqkeycovar' must be a logical value..", call. = FALSE)
  } 
  else if (type == "sortregs"){
    if (!is.singleLogical(object1)) 
      stop("'sortregs' must be logical value.", call. = FALSE)
  } 
  else if (type == "checklog"){
    if (!is.singleLogical(object1)) 
      stop("'checklog' must be logical value.", call. = FALSE)
  } 
  else if(type=="pValueCutoff"){
    if(!is.singleNumber(object1) || object1 > 1 || object1 < 0) 
      stop("'pValueCutoff' should be an integer or numeric value >=0 and <=1 !", 
           call. = FALSE)
  } 
  else if(type=="pAdjustMethod"){
    tp <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
    if(!is.singleString(object1) || !(object1 %in% tp)) 
      stop("'pAdjustMethod' should be any one of: ", 
           paste(tp, collapse = ", "), call. = FALSE)
  } 
  else if (type == "verbose"){
    if (!is.singleLogical(object1)) 
      stop("'verbose' must be logical value.", call. = FALSE)
  } 
  else if (type == "stepFilter"){
    if (!is.singleLogical(object1)) 
      stop("'stepFilter' must be logical value.", call. = FALSE)
  } 
  else if (type == "ntop"){
    if (!is.null(object1)){
      if (!is.singleInteger(object1) || object1 <= 0) 
        stop("'ntop' should be an integer value > 0.", call. = FALSE)
    }
  }
  else if (type == "width"){
    if (!is.singleNumber(object1)) 
      stop("'width' must be a single numeric values.", call. = FALSE)
  }
  else if (type == "height"){
    if (!is.singleNumber(object1)) 
      stop("'height' must be a single numeric values.", call. = FALSE)
  }
  else if (type == "endpoint"){
    if (!is.null(object1) && !is.singleNumber(object1)) 
      stop("'endpoint' must be a numeric value.", call. = FALSE)
  } 
  else if (type == "aSample"){
    if (!object1 %in% rownames(object2)) 
      stop("'aSample' must be present inside the 'TNS' class object, 
           slots 'tni' and 'survivalData'", 
           call. = FALSE)
  } 
  else if (type == "refsamp"){
    if (!is.null(object1)){
      if (!object1 %in% colnames(object2)) 
        stop("'refsamp' must be samples present in 'gexp' inside 
             TNS object1.", call. = FALSE)
    }
    } 
  else if (type == "plotype"){
    if(is.singleString(object1)) 
      tp <- c("2D","3D")
    if(!is.singleString(object1)){
      stop("'plotype' must be a single string.", call. = FALSE) 
    } else if (!(object1 %in% tp)){
      stop("'plotype' must be one of '",paste(tp, collapse ="', '"),"'", call. = FALSE)
    }
  } 
  else if (type == "xlim"){
    if(!is.numeric(object1) || length(object1) != 2) 
      stop("'xlim' must be a numeric vector of length 2.", call. = FALSE)
  } 
  else if (type == "xlim_log"){
    if(!is.numeric(object1) || length(object1) != 2) 
      stop("'xlim' must be a numeric vector of length 2.", call. = FALSE)
    if(any(object1<=0))
      stop("'xlim' must be > 0 in log space.", call. = FALSE)
  }
  else if (type == "ylim"){
    if(!is.numeric(object1) || length(object1) != 2) 
      stop("'ylim' must be a numeric vector of length 2.", call. = FALSE)
  } 
  else if (type == "ylim_log"){
    if(!is.numeric(object1) || length(object1) != 2) 
      stop("'ylim' must be a numeric vector of length 2.", call. = FALSE)
    if(any(object1<=0))
      stop("'ylim' must be > 0 in log space.", call. = FALSE)
  } 
  else if (type == "zlim"){
    if(!is.numeric(object1) || length(object1) != 2) 
      stop("'zlim' must be a numeric vector of length 2.", call. = FALSE)
  }
  else if (type == "zlim_log"){
    if(!is.numeric(object1) || length(object1) != 2) 
      stop("'zlim' must be a numeric vector of length 2.", call. = FALSE)
    if(any(object1<=0))
      stop("'zlim' must be > 0 in log space.", call. = FALSE)
  }
  else if (type == "xlim_reg"){
    if(!is.numeric(object1) || length(object1) != 2) 
      stop("'xlim' must be a numeric vector of length 2.", call. = FALSE)
    if(min(object1)<(-2) || max(object1)>2)
      stop("'xlim' must be in the range [-2,2].", call. = FALSE)
  } 
  else if (type == "ylim_reg"){
    if(!is.numeric(object1) || length(object1) != 2) 
      stop("'ylim' must be a numeric vector of length 2.", call. = FALSE)
    if(min(object1)<(-2) || max(object1)>2)
      stop("'ylim' must be in the range [-2,2].", call. = FALSE)
  }
  else if (type == "hlim_log"){
    if(!is.numeric(object1) || length(object1) != 2) 
      stop("'hlim' must be a numeric vector of length 2.", call. = FALSE)
    if(any(object1<=0))
      stop("'hlim' must be > 0 in log space.", call. = FALSE)
  }
  else if (type == "dualreg"){
    if(!is.singleString(object1)) 
      stop("'dualreg' must be a single string.", call. = FALSE)
    tp <- unlist(strsplit(object1, split = "~", fixed=TRUE))
    if(length(tp)!=2)
      stop("'dualreg' does not follow the expected syntax, e.g., 'reg1~reg2'", 
           call. = FALSE)
  } 
  else if (type == "panelWidths"){
    if (!is.numeric(object1) || length(object1) != 3) 
      stop("'panelWidths' must be a numeric vector of length 3.", 
           call. = FALSE)
    if (object1[1] == 0 || object1[3] == 0) 
      stop("The width of the first and third panels cannot be 0.", 
           call. = FALSE)
  } 
  else if (type == "TNI"){
    if(class(object1)!='TNI')
      stop("NOTE: 'tni' must be an object of class 'TNI'!", call. = FALSE)
    if (object1@status["Preprocess"] != "[x]") 
      stop("NOTE: TNI object requires preprocessing in the RTN package!")
    if (object1@status["Permutation"] != "[x]") 
      stop("NOTE: TNI object needs to be evaluated by 'tni.permutation' in the RTN package!")
    if (object1@status["DPI.filter"] != "[x]") 
      stop("NOTE: TNI object needs to be evaluated by 'tni.dpi.filter' in the RTN package!")
  }
  else if (type == "TNSpreprocess"){
    if(class(object1)!='TNS')
      stop("NOTE: 'tns' must be an object of class 'TNS'!", call. = FALSE)
    if(object1@status["Preprocess"] != "[x]") 
      stop("NOTE: TNS object requires preprocessing!", call. = FALSE)
  }
  else if (type == "TNSgsea2"){
    if(class(object1)!='TNS')
      stop("NOTE: 'tns' must be an object of class 'TNS'!", call. = FALSE)
    if(object1@status["Preprocess"] != "[x]") 
      stop("NOTE: TNS object requires preprocessing!", call. = FALSE)
    if(object1@status["GSEA2"] != "[x]") 
      stop("NOTE: TNS object needs to be evaluated by 'tnsGSEA2'!", 
           call. = FALSE)
  } 
  else if(type == "MBR"){
    if(class(object1) != "MBR")
      stop("'mbr' must be an object of class 'MBR'", call. = FALSE)
    status <- mbrGet(object1, "status")
    if(status["Association"] != "[x]") {
      stop("'mbr' must be evaluated by 'mbrAssociation'", call. = FALSE)
    }
  } 
  else if(type == "dualregs"){
    if(class(object1) != "character" & !is.null(object1)) {
      stop("'dualregs' must be a character vector or NULL.", call. = FALSE)
    }
    if(!is.null(object1)){
      duals <- mbrGet(object2, "dualRegulons")
      if(!all(object1 %in% duals)){
        stop("All elements of 'dualregs' must be present in the 'mbr' object",
             call. = FALSE)
      }
    }
  } 
  else if(type == "CBpal"){
    valid.pals <- c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", 
                    "Spectral", "Accent", "Dark2", "Paired", "Pastel1", "Pastel2", 
                    "Set1", "Set2", "Set3", "Blues", "BuGn", "BuPu", "GnBu", "Greens", 
                    "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", 
                    "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd")
    if (!(object1 %in% valid.pals)) {
      stop("'pal' must be a valid palette of RColorBrewer.", call. = FALSE)
    }
  } 
  else if(type == "nSections"){
    if (!(object1 %in% 1:3)) {
      stop("'nSections' must be a number between 1 and 3.", call. = FALSE)
    }
  } 
  else if(type == "setMid"){
    if (!is.singleLogical(object1)) 
      stop("'setMid' must be logical value.", call. = FALSE)
  }
  else if(type == "cols") {
    if(!is.color(object1))
      stop("NOTE: 'cols' should be a vector with valid colors!", 
           call.=FALSE)
  }
  else if(type == "hcols") {
    if(!is.color(object1) || length(object1)!=2)
      stop("NOTE: 'hcols' should be a vector (length = 2) with valid colors!", 
           call.=FALSE)
  }
}

is.singleNumber <- function(para){
    (is.integer(para) || is.numeric(para)) && length(para) == 1L && !is.na(para)
}
is.singleInteger <- function(para){
    lg <- (is.integer(para) || is.numeric(para)) && length(para) == 1L && !is.na(para)
    if (lg) lg <- (para/ceiling(para)) == 1
    return(lg)
}
is.singleString <- function(para){
    is.character(para) && length(para) == 1L && !is.na(para)
}
is.singleLogical <- function(para){
    is.logical(para) && length(para) == 1L && !is.na(para)
}
all.binaryValues <- function(para){
    all(para %in% c(0, 1, NA))
}
all.integerValues <- function(para){
    lg <- (all(is.integer(para)) || all(is.numeric(para))) && !any(is.na(para))
    if (lg) lg <- all((para/ceiling(para)) == 1)
    return(lg)
}
all.characterValues <- function(para){
    all(is.character(para)) && !any(is.na(para))
}
is.color <- function(x){
  res <- try(col2rgb(x),silent=TRUE)
  return(!"try-error"%in%class(res))
}
