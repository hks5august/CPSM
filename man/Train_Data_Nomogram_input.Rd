\name{Train_Data_Nomogram_input}
\alias{Train_Data_Nomogram_input}
\docType{data}
\title{
Train Data Nomogram input
}
\description{
This is Train Data Nomogram input.
}
\usage{data("Train_Data_Nomogram_input")}
\format{
  A data frame with 158 observations on the following 34 variables.
  \describe{
    \item{\code{Age}}{a numeric vector}
    \item{\code{subtype}}{a character vector}
    \item{\code{gender}}{a character vector}
    \item{\code{race}}{a character vector}
    \item{\code{ajcc_pathologic_tumor_stage}}{a logical vector}
    \item{\code{histological_type}}{a character vector}
    \item{\code{histological_grade}}{a character vector}
    \item{\code{treatment_outcome_first_course}}{a character vector}
    \item{\code{radiation_treatment_adjuvant}}{a character vector}
    \item{\code{sample_type}}{a character vector}
    \item{\code{type}}{a character vector}
    \item{\code{OS}}{a numeric vector}
    \item{\code{OS.time}}{a numeric vector}
    \item{\code{DSS}}{a numeric vector}
    \item{\code{DSS.time}}{a numeric vector}
    \item{\code{DFI}}{a numeric vector}
    \item{\code{DFI.time}}{a numeric vector}
    \item{\code{PFI}}{a numeric vector}
    \item{\code{PFI.time}}{a numeric vector}
    \item{\code{OS_month}}{a numeric vector}
    \item{\code{OS_month}}{a numeric vector}
    \item{\code{ALG6}}{a numeric vector}
    \item{\code{ARHGAP11A}}{a numeric vector}
    \item{\code{DESI1}}{a numeric vector}
    \item{\code{GALNT7}}{a numeric vector}
    \item{\code{GJD3}}{a numeric vector}
    \item{\code{GPC1}}{a numeric vector}
    \item{\code{H2BC5}}{a numeric vector}
    \item{\code{HOXD12}}{a numeric vector}
    \item{\code{RNF185}}{a numeric vector}
    \item{\code{TANGO2}}{a numeric vector}
    \item{\code{UNG}}{a numeric vector}
    \item{\code{ZNF648}}{a numeric vector}
    \item{\code{PI}}{a numeric vector}
  }
}
\details{
This is Train Data Nomogram input.
}
\value{
  A data frame with columns X and Y.
}
\source{
https://portal.gdc.cancer.gov/
}

\examples{
data(Train_Data_Nomogram_input)
## maybe str(Train_Data_Nomogram_input)
## plot(Train_Data_Nomogram_input)
}
\keyword{datasets}
