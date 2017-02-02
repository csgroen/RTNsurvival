# Unit tests fot TNS-class methods
test_tns <- function(){
    #-- loading data
  data(survival.data)
  data(dt4rtn, package="RTN")
  tfs4test <- dt4rtn$tfs[c("PTTG1", "FOXM1")]
  geneIDs <- dt4rtn$gexpIDs
  samples <- colnames(dt4rtn$gexp)[1:50]
  gexp <- dt4rtn$gexp[, samples]
  survival.data <- survival.data[samples, c(1,2,3)]
  
   #-- tni creation
  rtni <- new("TNI", gexp = gexp, transcriptionFactors = tfs4test)
  rtni <- tni.preprocess(rtni, gexpIDs = geneIDs, verbose = F)
  rtni <- tni.permutation(rtni, nPermutations = 1, verbose = F)
  rtni <- tni.bootstrap(rtni, nBootstraps = 1, verbose = F)
  rtni <- tni.dpi.filter(rtni, verbose = F)
  
   #-- tns creation
  
  rtns <- tns.preprocess(rtni, survival.data, keycovar = c("Grade"))
  checkTrue(class(rtns) == "TNS")
  
  #-- GSEA2
  rtns <- tns.gsea2(rtns)
  EScores <- tns.get(rtns, what = "EScores")
  checkTrue(is.list (EScores) && length(EScores) == 3)
  
}
