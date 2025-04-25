\name{Train_results}
\docType{data}
\alias{Train_results}
\title{Training Results for Survival Prediction}
\description{
A data frame containing survival prediction results for training samples.
The data includes actual risk group labels, predicted risk groups, probabilities for each class, and overall survival data.
}
\format{
A data frame with 158 observations and 8 variables:
\describe{
  \item{Sample_ID}{Character vector with the unique sample identifier (e.g., "TCGA-E1-5318-01").}
  \item{Actual}{Factor indicating the actual risk group: "High_Risk" or "Low_Risk".}
  \item{Predicted_Risk_Group}{Factor indicating the predicted risk group: "High_Risk" or "Low_Risk".}
  \item{High_Risk_Prob}{Numeric value representing the probability of the sample being classified as High Risk.}
  \item{Low_Risk_Prob}{Numeric value representing the probability of the sample being classified as Low Risk.}
  \item{Prediction_Prob}{Numeric value representing the maximum class probability (confidence in prediction).}
  \item{OS_month}{Integer representing the overall survival time in months.}
  \item{OS_event}{Integer indicating the survival event: 1 = death, 0 = censored.}
}
}
\usage{data(Train_results)}
\source{
The data was derived from TCGA (The Cancer Genome Atlas) using survival prediction models.
}

