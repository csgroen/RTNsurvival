\name{RTNsurvival-package}
\alias{RTNsurvival-package}
\alias{RTNsurvival}
\docType{package}

\title{
Performs survival analysis using transcriptional networks inferred by the RTN package.
}

\description{
This package provides classes and methods to perform survival analysis using transcriptional networks inferred by the RTN package, including Kaplan-Meier and multivariate survival analysis using Cox's regression model.
}

\details{
\tabular{ll}{
Package: \tab RTNsurvival\cr
Type: \tab Package\cr
Depends: \tab R (>= 3.3.0), methods\cr
Imports: \tab RTN, survival, RColorBrewer\cr
Suggests: \tab Fletcher2013b, knitr, rmarkdown, BiocStyle, RUnit, BiocGenerics\cr
License: \tab Artistic-2.0\cr
biocViews: \tab NetworkInference, NetworkEnrichment, GeneRegulation, GeneExpression, GraphAndNetwork, GeneSetEnrichment\cr
}
}

\section{Index}{

\tabular{ll}{
\link{TNS-class}: \tab an S4 class for survival survival analysis using RTN transcriptional networks.\cr
\link{tni2tnsPreprocess}: \tab a preprocessing method for objects of class TNS.\cr
\link{tnsGSEA2}: \tab two-tailed gene set enrichment analysis (GSEA) over a list of regulons.\cr
\link{tnsPlotGSEA2}: \tab plot results from the two-tailed GSEA.\cr
\link{tnsKM}: \tab Kaplan-Meier analysis for TNS class objects.\cr
\link{tnsPlotKM}: \tab Kaplan-Meier plots for TNS class objects.\cr
\link{tnsCox}: \tab Cox regression analysis for TNS class objects.\cr
\link{tnsPlotCox}: \tab Cox plots for TNS class objects.\cr
\link{tnsGet}: \tab Get information from slots in a TNS object.\cr
\link{tnsInteraction}: \tab A generic call to 'tnsCoxInteraction' and 'tnsKmInteraction'.\cr
\link{tnsKmInteraction}: \tab Kaplan-Meier analysis for dual regulons.\cr
\link{tnsPlotKmInteraction}: \tab Plot results from Kaplan-Meier analysis for dual regulons.\cr
\link{tnsCoxInteraction}: \tab Cox regression analysis for dual regulons.\cr
\link{tnsPlotCoxInteraction}: \tab Plot results from Cox regression analysis for dual regulons.\cr
\link{tnsPlotGSEA2}: \tab Plot 2-tailed GSEA for a sample from a TNS.\cr
}

Further information is available in the vignettes by typing \code{vignette("RTNsurvival")}. Documented 
topics are also available in HTML by typing \code{help.start()} and selecting the RTNsurvival package 
from the menu.

}

\author{
Clarice S. Groeneveld, Vinicius S. Chagas, Gordon Robertson, ..., Kerstin Meyer, Mauro A. A. Castro
}

\references{
Fletcher M.N.C. et al., \emph{Master regulators of FGFR2 signalling and breast cancer risk.} Nature Communications, 4:2464, 2013.

Castro M.A.A. et al., \emph{Regulators of genetic risk of breast cancer identified by integrative network analysis.} Nature Genetics, 48:12-21, 2016.
}

\keyword{ package }
