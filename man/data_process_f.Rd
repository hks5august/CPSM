% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_process_f.R
\name{data_process_f}
\alias{data_process_f}
\title{Data Processing function}
\arguments{
\item{data}{A data frame containing clinical and expression data.}

\item{col_num}{An integer indicating the column number where expression data
starts.}

\item{surv_time}{A character string specifying the name of the survival time
column in the data.}
}
\value{
A data frame combining processed clinical and expression data.
Includes a new column \code{OS_month} for survival time in months, and retains
only rows with positive \code{OS_month}.
}
\description{
This function converting OS time (in days) into months and
removing samples where OS time information is missing. Note: In the Example
data OS time labelled is "OS.time" and OS event labelled as "OS"
}
\examples{
library(SummarizedExperiment)
data(Example_TCGA_LGG_FPKM_data, package = "CPSM")
combined_df <- cbind(as.data.frame(colData(Example_TCGA_LGG_FPKM_data))
[, -ncol(colData(Example_TCGA_LGG_FPKM_data))],
t(as.data.frame(assay(Example_TCGA_LGG_FPKM_data, "expression"))))
data_process_f(data = combined_df,
  col_num = 20, surv_time = "OS.time"
)

}
