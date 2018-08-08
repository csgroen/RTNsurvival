# Unit tests fot TNS-class methods
test_tns <- function(){
  #--- dataset for demonstration
  data("stni", package = "RTN")
  data(survival.data)
  # new TNS object
  stns <- tni2tnsPreprocess(stni, survivalData = survival.data, 
                            keycovar = c('Grade','Age'), 
                            time = 1, event = 2)
  checkTrue(class(stns) == "TNS")
  #-- GSEA2
  stns <- tnsGSEA2(stns, verbose=FALSE)
  regulonActivity <- tnsGet(stns, what = "regulonActivity")
  checkTrue(is.list(regulonActivity) && length(regulonActivity) >= 3)
}
