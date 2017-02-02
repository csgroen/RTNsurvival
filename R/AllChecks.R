##This function is used for argument checking
.tns.checks <- function(object1, object2 = NULL, type)
{
    #---
    if (type == "survivalData") {
        if (!is.data.frame(object1))
            stop("'survivalData' must be a data frame.",call.=FALSE)
        if(is.null(colnames(object1)))
            stop("columns in 'survivalData' must be named.",call.=FALSE)
    }
    #---
    else if (type == "Time")
    {
        if( !is.singleInteger(object1) && !is.singleString(object1))
            stop("'time' should be either a character or integer value !\n",
                 call.=FALSE)
        if (is.character(object1)){
            if(!object1 %in% colnames(object2))
                stop("'time' parameter must be present in 'survivalData' 
                     colnames.",call.=FALSE)
            object1 <- which(colnames(object2) == object1)
        }
        vals <- object2[,object1]
        if (!is.numeric(vals)){
            stop("'time' data must be numeric.",call.=FALSE) 
        } else if(any(vals < 0, na.rm=TRUE)){
            stop("values in 'time' data must be >= 0",call.=FALSE)
        }
        return(object1)
    }
    
    else if (type == "Event")
    {
        if( !is.singleInteger(object1) && !is.singleString(object1) )
            stop("'event' should be either a character or integer 
                 value !\n",call.=FALSE)
        if (is.character(object1)){
            if(!object1 %in% colnames(object2))
                stop("'event' parameter must be present in 'survivalData' 
                     colnames.",call.=FALSE)
            object1 <- which(colnames(object2) == object1)
        }
        vals <- object2[,object1]
        if (!is.numeric(vals)){
            stop("'event' data must be numeric.",call.=FALSE) 
        } else if(!all.binaryValues(vals)){
            stop("'event' data must be either binary or logical.",call.=FALSE)
        }
        return(object1)
    }
    
    #---
    else if (type == "Keycovars") {
        if (!is.character(object1))
            stop ("'keycovar' must be a character vector.",call.=FALSE)
        if (!all(object1 %in% colnames(object2)))
            stop("All strings in 'keycovar' must be colnames in 'survivalData'",
                 call.=FALSE)
        for (col in object1) {
            if (!is.numeric(object2[,col]))
                stop("All values in 'keycovar' columns must be numeric.",
                     call.=FALSE)
        }
    }
    
    #---
    else if (type == "Samples") {
        if(is.null(object1)){
            object1 <- rownames(object2)
        } else {
            if (!is.character(object1))
                stop("'samples' must be a character vector.",call.=FALSE)
            
            if(!all(object1 %in% rownames(object2)))
                stop("All strings in 'samples' must appear in the rownames of 
             'survivalData'",call.=FALSE)
        }
        return(object1)
    }
    #---
    else if (type == "survival_cox")
    {
        if (nrow(object1) < 75)
            warning("If the number of samples in survivalData is too small, 
              coxph function may not converge.")
    }
    #---
    else if (type == "Path")
    {
        if (!is.singleString(object1))
            stop("'fpath' must be a single character.",call.=FALSE)
        
        if (!dir.exists(object1))
            stop("fpath does not lead to an existing directory.",call.=FALSE)
        
    }
    #---
    else if (type == "Fname")
    {
        if (!is.singleString(object1))
            stop("'fname' must be a single character.",call.=FALSE)
    }
    #---
    else if (type == "Ylab")
    {
        if (!is.singleString(object1))
            stop("'ylab' must be a single character.",call.=FALSE)
    }
    
    #---
    else if (type == "Xlab")
    {
        if (!is.singleString(object1))
            stop("'xlab' must be a single character.",call.=FALSE)
    }
    
    #---
    else if (type == "Regs")
    {
        if (!is.null(object1)){
            if (!all.characterValues(object1))
                stop("'regs' must be a character vector.",call.=FALSE)
        }
    }
    
    #---
    else if (type == "Attribs")
    {
        if (!is.null(object1)){
            if(is.list(object1))
                object1<-unlist(object1)
            if(is.character(object1)){
                if(!all(object1%in%colnames(object2)))
                    stop("all 'attribs' must be listed in the 'survivalData' 
colnames, at the 'tns' object!",call.=FALSE)
            } else if( !all.integerValues(object1) ){
                stop("'attribs' must be either a vector or a list of vectors, 
             with either integer or character values!",call.=FALSE)
            }
        }
    }
    
    #---
    else if (type == "Pal")
    {
        if ( !object1 %in% c("red","blue","redblue")){
            len <- (object2*2)+1
            message <- paste("'pal' must be 'red', 'blue' or 'redblue', 
                       or a vector of colors of length ", len)
            if (length(object1) < len)
                stop (message,call.=FALSE)
        }
    }
    
    #---
    else if (type == "ExcludeMid")
    {
        if (!is.singleLogical(object1))
            stop ("'excludeMid' must be a logical value.",call.=FALSE)
    }
    
    #---
    else if (type == "FlipCols")
    {
        if (!is.singleLogical(object1))
            stop ("'flipcols' must be a logical value.",call.=FALSE)
    }
    
    #---
    else if (type == "PlotPDF")
    {
        if (!is.singleLogical(object1))
            stop ("'plotpdf' must be a logical value.",call.=FALSE)
    }
    
    #---
    else if (type == "PlotBatch")
    {
        if (!is.singleLogical(object1))
            stop ("'plotbatch' must be a logical value.",call.=FALSE)
    }
    
    #---
    else if (type == "QQCovar")
    {
        if (!is.singleLogical(object1))
            stop ("'qqkeycovar' must be a logical value..",call.=FALSE)
    }
    
    #---
    else if (type == "SortRegs")
    {
        if (!is.singleLogical(object1))
            stop ("'sortregs' must be logical value.",call.=FALSE)
    }
    
    #---
    else if (type == "Log")
    {
        if (!is.singleLogical(object1))
            stop ("'log' must be logical value.",call.=FALSE)
    }
    
    #---
    else if(type=="pValueCutoff") {
        if(!is.singleNumber(object1) || object1>1 || object1<0)
            stop("'pValueCutoff' should be an integer or numeric value >=0 
                 and <=1 !",call.=FALSE)
    }
    
    #---
    else if(type=="pAdjustMethod") {
        tp <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", 
                "none")
        if(!is.singleString(object1) || !(object1 %in% tp))
            stop("'pAdjustMethod' should be any one of: ",paste(tp, 
                                                                collapse = ", "),call.=FALSE)
    }
    
    #---
    else if (type == "Verbose")
    {
        if (!is.singleLogical(object1))
            stop ("'verbose' must be logical value.",call.=FALSE)
    }
    
    #---
    else if(type=="Ntop") {
        if(!is.null(object1)){
            if (!is.singleInteger(object1) || object1<=0)
                stop("'ntop' should be an integer value > 0.",call.=FALSE) 
        }
    }
    
    #---
    else if (type == "WidthHeight")
    {
        if (!is.singleNumber(object1) || !is.singleNumber(object2))
            stop("'width' and 'height' must be single numeric values.",
                 call.=FALSE)
    }
    
    #---
    else if (type == "EndPoint")
    {
        if (!is.singleNumber(object1))
            stop("'endPoint' must be a numeric value.",call.=FALSE)
    }
    
    #---
    else if (type == "aSample")
    {
        if (!object1 %in% rownames(object2))
            stop ("'aSample' must be present inside the 'TNS' class object, 
                  slots 'tni' and 'survivalData'",call.=FALSE)
    }
    
    #--
    else if (type == "Refsamp")
    {
        if (!is.null(object1))
        {
            if (!object1 %in% colnames(object2))
                stop ("'refsamp' must be samples present in 'gexp' inside 
                      TNS object1.",call.=FALSE)
        }
    }
    
    #---
    else if (type == "Xlim")
    {
        if (!is.numeric(object1) || length(object1) != 2)
            stop ("'xlim' must be a numeric vector of length 2.")
        
        if (object1[1] > object1[2])
            stop("The first element of 'xlim' vector must be smaller than 
                 the second.",call.=FALSE)
    }
    
    #---
    else if (type == "panelWidths")
    {
        if (!is.numeric(object1) || length(object1) != 3)
            stop ("'panelWidths' must be a numeric vector of length 2.")
    }
    
    #---
    else if (type == "status")
    {
        if(object1@status["Preprocess"]!="[x]")
            stop("NOTE: TNS object requires preprocessing!")
        if(object1@status["GSEA2"]!="[x]")
            stop("NOTE: TNS object needs to be evaluated by 'tnsGSEA2'!")
    }
}


is.singleNumber<-function(para){
    (is.integer(para) || is.numeric(para)) && length(para)==1L && !is.na(para)
}
is.singleInteger<-function(para){
    lg <- (is.integer(para) || is.numeric(para)) && length(para)==1L && !is.na(para)
    if(lg) lg <- (para / ceiling(para)) == 1
    return(lg)
}
is.singleString<-function(para){
    is.character(para) && length(para) == 1L && !is.na(para)
}
is.singleLogical<-function(para){
    is.logical(para) && length(para) == 1L && !is.na(para)
}
all.binaryValues<-function(para){
    all( para %in% c(0, 1, NA) )
}
all.integerValues<-function(para){
    lg <- ( all(is.integer(para)) || all(is.numeric(para)) ) && !any(is.na(para))
    if(lg) lg <- all ( (para / ceiling(para)) == 1 )
    return(lg)
}
all.characterValues<-function(para){
    all(is.character(para)) && !any(is.na(para))
}

