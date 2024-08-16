\name{Train_Clin}
\alias{Train_Clin}
\docType{data}
\title{
Train Clin
}
\description{
This is Train Clin.
}
\usage{data("Train_Clin")}
\format{
  A data frame with 158 observations on the following 20 variables.
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
  }
}
\details{
This is Train Clin.
}
\value{
  A data frame with columns X and Y.
}
\source{
https://portal.gdc.cancer.gov/
}

\examples{
data(Train_Clin)
## maybe str(Train_Clin)
## plot(Train_Clin)
}
\keyword{datasets}
