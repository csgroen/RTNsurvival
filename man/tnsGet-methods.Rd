% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllMethods.R
\docType{methods}
\name{tnsGet,TNS-method}
\alias{tnsGet,TNS-method}
\alias{tnsGet}
\title{Get information from slots in a TNS object}
\usage{
\S4method{tnsGet}{TNS}(tns, what)
}
\arguments{
\item{tns}{A \linkS4class{TNS} object.}

\item{what}{A string specifying what should be retrieved from the object.
Options: 'status','survivalData', 'regulonActivity', 'TNI', 'para', 'kmTable', 'kmFit', 
'coxTable', 'coxFit', 'kmInteractionTable', 'kmInteractionFit', 
'coxInteractionTable', 'coxInteractionFit', and 'regulatoryElements'.}
}
\value{
Content from slots in the \linkS4class{TNS} object.
}
\description{
Get information from individual slots in a TNS object and any
available results from a previous analysis.
}
\examples{
# load survival data
data(survival.data)

# load TNI-object
data(stni, package = "RTN")

stns <- tni2tnsPreprocess(stni, survivalData = survival.data, 
keycovar = c('Grade','Age'), time = 1, event = 2)
stns <- tnsGSEA2(stns)
regulonActivity <- tnsGet(stns, 'regulonActivity')

}
