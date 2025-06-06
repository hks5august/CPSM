# CPSM: Cancer Patient Survival model
# Introduction
CPSM is an R package that provides a comprehensive computational pipeline for predicting survival probabilities and risk groups in cancer patients. It includes dedicated modules to perform key steps such as data preprocessing, training/test splitting, and normalization. CPSM enables feature selection through univariate cox-regression survival analysis, feature selection though LASSO method,  and calculates a LASSO-based Prognostic Index (PI) score. It supports the development of predictive models using different feature sets and offers a suite of visualization tools, including survival curves based on predicted probabilities, barplots of predicted mean and median survival times, Kaplan-Meier (KM) plots overlaid with individual survival predictions, and nomograms for estimating 1-, 3-, 5-, and 10-year survival probabilities. Together, these functionalities make CPSM a powerful and versatile tool for survival analysis in cancer research.

# Installation from the Bioconductor

## Install development version (with latest updates) from the Bioconductor

To install this package, start R (version "4.5") and enter:
```{r, warning=FALSE, message=FALSE, eval=FALSE}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("CPSM")

```

## Install last released version from the Bioconductor

To install this package, start R (version "4.4") and enter the code provided:
```{r, warning=FALSE, message=FALSE, eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("CPSM")
```


# Installation from Github:
```{r, warning=FALSE, message=FALSE, eval=FALSE}
#Step1: First Install remote package
install.packages("remotes")
#load remotes package
library("remotes")
#Step2: install CPSM package
remotes::install_github("hks5august/CPSM", local = TRUE, , dependencies=TRUE)
```

# Input Data
The example input data object, **`Example_TCGA_LGG_FPKM_data`**, contains data for **184 LGG cancer samples** as rows and various features as columns. Gene expression data is represented in **FPKM values**. The dataset includes **11 clinical and demographic features**, **4 types of survival data** (with both time and event information), and **19,978 protein-coding genes**. The clinical and demographic features in the dataset include `Age`, `subtype`, `gender`, `race`, `ajcc_pathologic_tumor_stage`, `histological_type`, `histological_grade`, `treatment_outcome_first_course`, `radiation_treatment_adjuvant`, `sample_type`, and `type`. The four types of survival data included are **Overall Survival (OS)**, **Progression-Free Survival (PFS)**, **Disease-Specific Survival (DSS)**, and **Disease-Free Survival (DFS)**. In the dataset, the columns labeled **OS**, **PFS**, **DSS**, and **DFS** represent event occurrences, while the columns **OS.time**, **PFS.time**, **DSS.time**, and **DFS.time** provide survival times (in days).

```{r , warning=FALSE, message=FALSE}
library(CPSM)
library(SummarizedExperiment)
set.seed(7) # set seed
data(Example_TCGA_LGG_FPKM_data, package = "CPSM")
Example_TCGA_LGG_FPKM_data
```

# Step 1- Data Processing 

## Description 
The **data_process_f** function converts OS time (in days) into months and removes samples where OS/OS.time information is missing. 
## Required inputs
To use this function, the input data should be provided in TSV format. Additionally, you need to define `col_num` (the column number at which clinical, demographic, and survival information ends, e.g., 20), `surv_time` (the name of the column that contains survival time information, e.g., `OS.time`), and `output` (the desired name for the output, e.g., "New_data"). 

## Example Code
```{r }
data(Example_TCGA_LGG_FPKM_data, package = "CPSM")
combined_df <- cbind(
  as.data.frame(colData(Example_TCGA_LGG_FPKM_data))
  [, -ncol(colData(Example_TCGA_LGG_FPKM_data))],
  t(as.data.frame(assay(
    Example_TCGA_LGG_FPKM_data,
    "expression"
  )))
)
New_data <- data_process_f(combined_df, col_num = 20, surv_time = "OS.time")
str(New_data[1:10])
```
## Outputs
After data processing, the output object **`New_data`** is generated, which contains 176 samples. This indicates that the function has removed 8 samples where OS/OS.time information was missing. Moreover, a new 21st column, **`OS_month`**, is added to the data, containing OS time values in months.

# Step 2 - Split Data into Training and Test Subset
## Description 
Before proceeding further, we need to split the data into training and test subsets for feature selection and model development. 
## Required inputs
The output from the previous step, **`New_data`**, serves as the input for this process. Next, you need to define the fraction (e.g., 0.9) by which to split the data into training and test sets. For example, setting `fraction = 0.9` will divide the data into 90% for training and 10% for testing. Additionally, you should specify names for the training and test outputs (e.g., `train_FPKM` and `test_FPKM`).

## Example Code
```{r}
data(New_data, package = "CPSM")
# Call the function
result <- tr_test_f(data = New_data, fraction = 0.9)
# Access the train and test data
train_FPKM <- result$train_data
str(train_FPKM[1:10])
test_FPKM <- result$test_data
str(test_FPKM[1:10])
```
## Outputs
After the train-test split, two new output objects are generated: **`train_FPKM`** and **`test_FPKM`**. The **`train_FPKM`** object contains 158 samples, while **`test_FPKM`** contains 18 samples. This indicates that the **`tr_test_f`** function splits the data in a 90:10 ratio.

# Step 3 - Data Normalization
## Description 
In order to select features and develop ML models, the data must be normalized. Since the expression data is available in terms of FPKM values, the **`train_test_normalization_f`** function will first convert the FPKM values into a log scale using the formula [log2(FPKM+1)], followed by quantile normalization. The training data will be used as the target matrix for the quantile normalization process. 
## Required inputs
For this function, you need to provide the training and test datasets obtained from the previous step (Train/Test Split). Additionally, you must specify the column number where clinical information ends (e.g., 21) in the input datasets. Finally, you need to define output names for the resulting datasets: **`train_clin_data`** (which contains only clinical information from the training data), **`test_clin_data`** (which contains only clinical information from the test data), **`train_Normalized_data_clin_data`** (which contains both clinical information and normalized gene expression values for the training samples), and **`test_Normalized_data_clin_data`** (which contains both clinical information and normalized gene expression values for the test samples).

## Example Code
```{r }
# Step 3 - Data Normalization
# Normalize the training and test data sets
data(train_FPKM, package = "CPSM")
data(test_FPKM, package = "CPSM")
Result_N_data <- train_test_normalization_f(
  train_data = train_FPKM,
  test_data = test_FPKM,
  col_num = 21
)
# Access the Normalized train and test data
Train_Clin <- Result_N_data$Train_Clin
Test_Clin <- Result_N_data$Test_Clin
Train_Norm_data <- Result_N_data$Train_Norm_data
Test_Norm_data <- Result_N_data$Test_Norm_data
str(Train_Clin[1:10])
str(Train_Norm_data[1:10])
```
## Outputs
After running the function, four outputs objects are generated: **`Train_Clin`** (which contains only clinical features from the training data), **`Test_Clin`** (which contains only clinical features from the test data), **`Train_Norm_data`** (which includes clinical features and normalized gene expression values for the training samples), and **`Test_Norm_data`** (which includes clinical features and normalized gene expression values for the test samples).


# Step 4a - Prognostic Index (PI)  Score Calculation
## Description 
To create a survival model, the next step is to calculate the Prognostic Index (PI) score. The PI score is based on the expression levels of features selected by the LASSO regression model and their corresponding beta coefficients. For example, suppose five features (**G1**, **G2**, **G3**, **G4**, **G5**) are selected by the LASSO method, and their associated coefficients are **B1**, **B2**, **B3**, **B4**, and **B5**, respectively. The PI score is then computed using the following formula:


**PI score = G1 * B1 + G2 * B2 + G3 * B3 + G4 * B4 + G5 * B5**

## Required inputs
For this function, you need to provide the normalized training data object  (**Train_Norm_data**) and test data object (**Test_Norm_data**) obtained from the previous step (**train_test_normalization_f**). Additionally, you must specify the column number (`col_num`) where clinical features end (e.g., 21), the number of folds (`nfolds`) for the LASSO regression method (e.g., 5), and the survival time (`surv_time`) and survival event (`surv_event`) columns in the data (e.g., `OS_month` and `OS`, respectively). The LASSO regression is implemented using the **`glmnet`** package. Finally, you need to define names of output object to store the results, which will include the selected LASSO features and their corresponding PI values. 

## Example Code
```{r, warning=FALSE, message=FALSE, fig.width=7, fig.height=4 }
# Step 4 - Lasso PI Score
data(Train_Norm_data, package = "CPSM")
data(Test_Norm_data, package = "CPSM")
Result_PI <- Lasso_PI_scores_f(
  train_data = Train_Norm_data,
  test_data = Test_Norm_data,
  nfolds = 5,
  col_num = 21,
  surv_time = "OS_month",
  surv_event = "OS"
)
Train_Lasso_key_variables <- Result_PI$Train_Lasso_key_variables
Train_PI_data <- Result_PI$Train_PI_data
Test_PI_data <- Result_PI$Test_PI_data
str(Train_PI_data[1:10])
str(Test_PI_data[1:10])
plot(Result_PI$cvfit)
```
## Outputs
The **`Lasso_PI_scores_f`** function generates the following outputs objects:
1. **`Train_Lasso_key_variables`**: A list of features selected by LASSO along with their beta coefficient values.
2. **`Train_Cox_Lasso_Regression_lambda_plot`**: The Lasso regression lambda plot.
3. **`Train_PI_data`**: This dataset contains the expression values of genes selected by LASSO along with the PI score in the last column for the training samples.
4. **`Test_PI_data`**: This dataset contains the expression values of genes selected by LASSO along with the PI score in the last column for the test samples.


# Step 4b - Univariate  Survival Significant Feature Selection
## Description 
In addition to the Prognostic Index (PI) score, the **`Univariate_sig_features_f`** function in the CPSM package allows for the selection of significant features based on univariate cox-regression survival analysis. This function identifies features with a p-value less than 0.05, which are able to stratify high-risk and low-risk survival groups. The stratification is done by using the median expression value of each feature as a cutoff.
## Required inputs
To use this function, you need to provide the normalized training (**Train_Norm_data**) and test (**Test_Norm_data**) dataset objects, which were obtained from the previous step (**train_test_normalization_f**). Additionally, you must specify the column number (`col_num`) where the clinical features end (e.g., 21), as well as the names of the columns containing survival time (`surv_time`, e.g., `OS_month`) and survival event information (`surv_event`, e.g., `OS`). Furthermore, you need to define output names for the resulting datasets that will contain the expression values of the selected genes. These outputs will be used to store the significant genes identified through univariate survival analysis.

## Example Code
```{r, warning=FALSE, message=FALSE }
# Step 4b - Univariate  Survival Significant Feature Selection.
data(Train_Norm_data, package = "CPSM")
data(Test_Norm_data, package = "CPSM")
Result_Uni <- Univariate_sig_features_f(
  train_data = Train_Norm_data,
  test_data = Test_Norm_data,
  col_num = 21,
  surv_time = "OS_month",
  surv_event = "OS"
)
Univariate_Suv_Sig_G_L <- Result_Uni$Univariate_Survival_Significant_genes_List
Train_Uni_sig_data <- Result_Uni$Train_Uni_sig_data
Test_Uni_sig_data <- Result_Uni$Test_Uni_sig_data
Uni_Sur_Sig_clin_List <- Result_Uni$Univariate_Survival_Significant_clin_List
Train_Uni_sig_clin_data <- Result_Uni$Train_Uni_sig_clin_data
Test_Uni_sig_clin_data <- Result_Uni$Test_Uni_sig_clin_data
str(Univariate_Suv_Sig_G_L[1:10])
```
## Outputs
The **`Univariate_sig_features_f`** function generates the following output objects:
1. **`Univariate_Surv_Sig_G_L`**: A table of univariate significant genes, along with their corresponding coefficient values, hazard ratio (HR) values, p-values, and C-Index values.
2. **`Train_Uni_sig_data`**: This dataset contains the expression values of the significant genes selected by univariate survival analysis for the training samples.
3. **`Test_Uni_sig_data`**: This dataset contains the expression values of the significant genes selected by univariate survival analysis for the test samples.


# Step 5 - Prediction model development for survival probability of patients
## Description 
After selecting significant features using LASSO or univariate survival analysis, the next step is to develop a machine learning (ML) prediction model to estimate the survival probability of patients. The **`MTLR_pred_model_f`** function in the CPSM package provides several options for building prediction models based on different feature sets. These options include:
- **Model_type = 1**: Model based on only clinical features
- **Model_type = 2**: Model based on PI score
- **Model_type = 3**: Model based on PI score + clinical features
- **Model_type = 4**: Model based on significant univariate features
- **Model_type = 5**: Model based on significant univariate features + clinical features


For this analysis, we are interested in developing a model based on the PI score (i.e., **Model_type = 2**). 
## Required inputs
To use this function, the following inputs are required:
1. **Training data with only clinical features**
2. **Test data with only clinical features**
3. **Model type** (e.g., **2** for a model based on PI score)
4. **Training data with PI score**
5. **Test data with PI score**
6. **`Clin_Feature_List`** (e.g., **Key_PI_list**), a list of features to be used for building the model
7. **`surv_time`**: The name of the column containing survival time in months (e.g., `OS_month`)
8. **`surv_event`**: The name of the column containing survival event information (e.g., `OS`)

These inputs will allow the **`MTLR_pred_model_f`** function to generate a prediction model for the survival probability of patients based on the provided data.

## Model for only Clinical features
## Example  Code
```{r, warning=FALSE, message=FALSE, error = TRUE }
data(Train_Clin, package = "CPSM")
data(Test_Clin, package = "CPSM")
data(Key_Clin_feature_list, package = "CPSM")
Result_Model_Type1 <- MTLR_pred_model_f(
  train_clin_data = Train_Clin,
  test_clin_data = Test_Clin,
  Model_type = 1,
  train_features_data = Train_Clin,
  test_features_data = Test_Clin,
  Clin_Feature_List = Key_Clin_feature_list,
  surv_time = "OS_month",
  surv_event = "OS"
)
survCurves_data <- Result_Model_Type1$survCurves_data
mean_median_survival_tim_d <- Result_Model_Type1$mean_median_survival_time_data
survival_result_bas_on_MTLR <- Result_Model_Type1$survival_result_based_on_MTLR
Error_mat_for_Model <- Result_Model_Type1$Error_mat_for_Model
```


## Model for PI
## Example  Code
```{r, warning=FALSE, message=FALSE, error = TRUE}
data(Train_Clin, package = "CPSM")
data(Test_Clin, package = "CPSM")
data(Train_PI_data, package = "CPSM")
data(Test_PI_data, package = "CPSM")
data(Key_PI_list, package = "CPSM")
Result_Model_Type2 <- MTLR_pred_model_f(
  train_clin_data = Train_Clin,
  test_clin_data = Test_Clin,
  Model_type = 2,
  train_features_data = Train_PI_data,
  test_features_data = Test_PI_data,
  Clin_Feature_List = Key_PI_list,
  surv_time = "OS_month",
  surv_event = "OS"
)
survCurves_data <- Result_Model_Type2$survCurves_data
mean_median_surviv_tim_da <- Result_Model_Type2$mean_median_survival_time_data
survival_result_b_on_MTLR <- Result_Model_Type2$survival_result_based_on_MTLR
Error_mat_for_Model <- Result_Model_Type2$Error_mat_for_Model
```


## Model for Clinical features + PI
## Example  Code
```{r, warning=FALSE, message=FALSE, error = TRUE}
data(Train_Clin, package = "CPSM")
data(Test_Clin, package = "CPSM")
data(Train_PI_data, package = "CPSM")
data(Test_PI_data, package = "CPSM")
data(Key_Clin_features_with_PI_list, package = "CPSM")
Result_Model_Type3 <- MTLR_pred_model_f(
  train_clin_data = Train_Clin,
  test_clin_data = Test_Clin,
  Model_type = 3,
  train_features_data = Train_PI_data,
  test_features_data = Test_PI_data,
  Clin_Feature_List = Key_Clin_features_with_PI_list,
  surv_time = "OS_month",
  surv_event = "OS"
)
survCurves_data <- Result_Model_Type3$survCurves_data
mean_median_surv_tim_da <- Result_Model_Type3$mean_median_survival_time_data
survival_result_b_on_MTLR <- Result_Model_Type3$survival_result_based_on_MTLR
Error_mat_for_Model <- Result_Model_Type3$Error_mat_for_Model
```


## Model for Univariate + Clinical features
## Example  Code
```{r, warning=FALSE, message=FALSE, error = TRUE}
data(Train_Clin, package = "CPSM")
data(Test_Clin, package = "CPSM")
data(Train_Uni_sig_data, package = "CPSM")
data(Test_Uni_sig_data, package = "CPSM")
data(Key_univariate_features_with_Clin_list, package = "CPSM")
Result_Model_Type5 <- MTLR_pred_model_f(
  train_clin_data = Train_Clin,
  test_clin_data = Test_Clin,
  Model_type = 4,
  train_features_data = Train_Uni_sig_data,
  test_features_data = Test_Uni_sig_data,
  Clin_Feature_List = Key_univariate_features_with_Clin_list,
  surv_time = "OS_month",
  surv_event = "OS"
)
survCurves_data <- Result_Model_Type5$survCurves_data
mean_median_surv_tim_da <- Result_Model_Type5$mean_median_survival_time_data
survival_result_b_on_MTLR <- Result_Model_Type5$survival_result_based_on_MTLR
Error_mat_for_Model <- Result_Model_Type5$Error_mat_for_Model
```

## Outputs
After implementing the **`MTLR_pred_model_f`** function, the following outputs are generated:


1. **Model_with_PI.RData**: This object contains the trained model based on the input data.
2. **survCurves_data**: This object contains the predicted survival probabilities for each patient at various time points. This data can be used to plot survival curves for patients.
3. **mean_median_survival_time_data**: Object containing the predicted mean and median survival times for each patient in the test data. This data can be used to generate bar plots illustrating the predicted survival times.
4. **Error_mat_for_Model**: Object containing performance metrics of the model based on the training and test data. It includes the following key performance scores:
   - **C-Index** = 0.81
These outputs allow you to evaluate the model's performance and visualize survival probabilities and survival times for the training test data.


# Step 6 - Survival curves/plots for individual patient
## Description
To visualize the survival of patients, we use the **`surv_curve_plots_f`** function, which generates survival curve plots based on the **`survCurves_data`** obtained from the previous step (after running the **`MTLR_pred_model_f`** function). This function also provides the option to highlight the survival curve of a specific patient. 

## Required Inputs
The function requires two inputs:
1. **Surv_curve_data**: The data object containing predicted survival probabilities for all patients.
2. **Sample ID**: The ID of the specific patient (e.g., `TCGA-TQ-A8XE-01`) whose survival curve you want to highlight.

## Example Code
```{r, warning=FALSE, message=FALSE, error = TRUE , fig.width=7, fig.height=4}
# Create Survival curves/plots for individual patients
data(survCurves_data, package = "CPSM")
plots <- surv_curve_plots_f(
  Surv_curve_data = survCurves_data,
  selected_sample = "TCGA-TQ-A7RQ-01"
)
# Print the plots
print(plots$all_patients_plot)
print(plots$highlighted_patient_plot)
```
## Outputs
After running the function, two output plots are generated:
1. **Survival curves for all patients** in the test data, displayed with different colors for each patient.
2. **Survival curves for all patients (in black)** with the selected patient highlighted in **red**.

These plots allow for easy visualization of individual patient survival in the context of the overall test data.


# Step 7 - Predicted mean and median survival time of individual patients 
## Description
To visualize the predicted survival times for patients, we use the **`mean_median_surv_barplot_f`** function, which generates bar plots for the mean and median survival times based on the data obtained from **Step 5** after running the **`MTLR_pred_model_f`** function. This function also provides the option to highlight a specific patient on the bar plot.

## Required Inputs
This function requires two inputs:
1. **surv_mean_med_data**: The data containing the predicted mean and median survival times for all patients.
2. **Sample ID**: The ID of the specific patient (e.g., `TCGA-TQ-A7RQ-01`) whose bar plot should be highlighted.

## Example Code
```{r, warning=FALSE, message=FALSE, error = TRUE , fig.width=7, fig.height=4}
data(mean_median_survival_time_data, package = "CPSM")
plots_2 <- mean_median_surv_barplot_f(
  surv_mean_med_data =
    mean_median_survival_time_data,
  selected_sample = "TCGA-TQ-A7RQ-01"
)
# Print the plots
print(plots_2$mean_med_all_pat)
print(plots_2$highlighted_selected_pat)
```

## Outputs

After running the function, two output bar plots are generated:
1. **Bar plot for all patients** in the test data, where the red-colored bars represent the mean survival time, and the cyan/green-colored bars represent the median survival time.
2. **Bar plot for all patients with a highlighted patient** (indicated by a dashed black outline). This plot shows that the highlighted patient has predicted mean and median survival times of 81.58 and 75.50 months, respectively.

These plots provide a clear comparison of the predicted survival times for all patients and the highlighted individual patient.

# Step 8 – Risk-Group Prediction of Test Samples Based on Selected Features
## Description
To predict the survival-based risk group of test samples (i.e., **high-risk** with shorter survival or **low-risk** with longer survival), we use the `predict_survival_risk_group_f()` function provided in the **CPSM** package. This function implements a **randomForestSRC**-based prediction approach for survival risk classification. Thhis function first defines actual risk groups in the training data using the **median overall survival time**:

- **Low-Risk**: Survival time > median  
- **High-Risk**: Survival time ≤ median

Multiple Random Survival Forest (RSF) models are then trained using different values for `ntree`:  
**10, 20, 50, 100, 250, 500, 750, 1000**.

The model with the best performance (e.g., highest accuracy) is selected automatically. This best-performing model is used to predict the risk group of test samples, along with prediction probabilities.

## Required Inputs

- `selected_train_data`: A data frame with normalized expression values for selected features and survival information (`OS_month`, `OS_event`) for the training set.
- `selected_test_data`: A data frame with normalized expression values for the same features for the test set.
- `Feature_List`: A character vector containing the names of selected features to be used in the model.

## Example Code

```{r, warning=FALSE, message=FALSE, error=TRUE}
# Load example data from CPSM package
data(Train_PI_data, package = "CPSM")
data(Test_PI_data, package = "CPSM")
data(Key_PI_list, package = "CPSM")

# Predict survival-based risk groups for test samples
Results_Risk_group_Prediction <- predict_survival_risk_group_f(
  selected_train_data = Train_PI_data,
  selected_test_data = Test_PI_data,
  Feature_List = Key_PI_list
)

#Performance of the best model on Training and Test data
Best_model_Prediction_results<- Results_Risk_group_Prediction$misclassification_results 
print(head(Best_model_Prediction_results))

#Prediction results of the best model on Training set
Test_results <-  Results_Risk_group_Prediction$Test_results #Prediction resulst on Test data
print(head(Test_results))

```

## Outputs

The output is a list that includes:
1. Best prediction model
2. Performance metrics (accuracy, sensitivity, specificity, Error rate, etc.) for training and test data
3. Predicted risk groups  with prediction probability values for training samples
4. Predicted risk groups  with prediction probability values for test samples

User can use these results for further validation and visualization, such as overlaying test sample survival curves on the training KM plot (see next step).

# Step 9 – Visual Overlay of Predicted Test Sample on Kaplan-Meier Curve

## Description
To visually evaluate how a specific test sample compares to survival risk groups defined in the training dataset, we use the `km_overlay_plot_f()` function.
This function overlays the predicted **survival curve of a selected test sample** onto the Kaplan-Meier (KM) survival plot derived from the training data. This visual comparison helps determine how closely the test sample aligns with population-level survival trends.
## Required Inputs
It requres requires following inputs
- `Train_results`:  
  A data frame containing predicted risk groups, survival times (`OS_month`), event status (`OS_event`), and additional training data.  
  **Row names** must correspond to sample IDs.

- `Test_results`:  
  A data frame with predicted risk groups and prediction probabilities for the test dataset.  
  **Row names** must correspond to sample IDs.

- `survcurve_te_data`:  
  A data frame with predicted survival probabilities over multiple time points for test samples (that we obtained from Step 5).  

- `selected_sample`:  
  The **sample ID** (matching a row in `Test_results`) for which the test survival curve should be plotted.

## Example Code

```{r, fig.height=6, fig.width=8, warning=FALSE, message=FALSE}
# Load example results
data(Train_results, package = "CPSM")
data(Test_results, package = "CPSM")
data(survCurves_data, package = "CPSM")

# Select a test sample to visualize
sample_id <- "TCGA-TQ-A7RQ-01"  

# Generate KM overlay plot
KM_plot <- km_overlay_plot_f(
  Train_results = Train_results,
  Test_results = Test_results,
  survcurve_te_data = survCurves_data,
  selected_sample = sample_id
)

# Display plot
KM_plot
```

## Outputs

1. Constructs Kaplan-Meier curves for training samples, stratified by **Low-Risk** and **High-Risk** groups based on median survival time.
2. Overlays the test sample’s predicted survival probability curve as a **dashed green line**.
3. Annotates the plot with Sample ID, Predicted risk group, Prediction probability

This visualization is useful for:
- Displaying individual patient patterns in a survival context of training samples
- Verifying predicted risk classifications

# Step 10 - Nomogram based on Key features
## Description
The **`Nomogram_generate_f`** function in the CPSM package allows you to generate a nomogram plot based on user-defined clinical and other relevant features in the data. For example, we will generate a nomogram using six features: Age, Gender, Race, Histological Type, Sample Type, and PI score. 

## Required Inputs
To create the nomogram, we need to provide the following inputs:
1. **Train_Data_Nomogram_input**: A dataset containing all the features, where samples are in the rows and features are in the columns.
2. **feature_list_for_Nomogram**: A list of features (e.g., Age, Gender, etc.) that will be used to generate the nomogram.
3. **surv_time**: The column name containing survival time in months (e.g., `OS_month`).
4. **surv_event**: The column name containing survival event information (e.g., `OS`).

## Example Code
```{r, warning=FALSE, message=FALSE, error = TRUE, fig.width=7, fig.height=6 }
data(Train_Data_Nomogram_input, package = "CPSM")
data(feature_list_for_Nomogram, package = "CPSM")
Result_Nomogram <- Nomogram_generate_f(
  data = Train_Data_Nomogram_input,
  Feature_List = feature_list_for_Nomogram,
  surv_time = "OS_month",
  surv_event = "OS"
)
C_index_mat <- Result_Nomogram$C_index_mat
```
## Outputs 
After running the function, the output is a nomogram that predicts the risk (e.g., Event risk such as death), as well as the 1-year, 3-year, 5-year, and 10-year survival probabilities for patients based on the selected features.The nomogram provides a visual representation to estimate the patient's survival outcomes over multiple time points, helping clinicians make more informed decisions.

# SessionInfo

As last part of this document, we call the function "sessionInfo()", 
which reports the version numbers of R and all the packages used 
in this session. It is good practice to always keep such a record 
as it will help to trace down what has happened in case that an R 
script ceases to work because the functions have been changed in 
a newer version of a package.


```{r}
sessionInfo()
```

# References
1. Kuhn, Max (2008). “Building Predictive Models in R Using the caret Package.” 
Journal of Statistical Software, 28(5), 1–26. doi:10.18637/jss.v028.i05, 
https://www.jstatsoft.org/index.php/jss/article/view/v028i05.
2. Bolstad B (2024). preprocessCore: A collection of pre-processing functions. 
R package version 1.66.0, https://github.com/bmbolstad/preprocessCore.
3. Horikoshi M, Tang Y (2018). ggfortify: Data Visualization Tools for 
Statistical Analysis Results. https://CRAN.R-project.org/package=ggfortify.
4. Therneau T (2024). A Package for Survival Analysis in R. R package version 
3.7-0, https://CRAN.R-project.org/package=survival.
5. Terry M. Therneau, Patricia M. Grambsch (2000). Modeling Survival Data: 
Extending the Cox Model. Springer, New York. ISBN 0-387-98784-3.
6. Kassambara, A., Kosinski, M., Biecek, P., & Scheipl, F. (2021). survminer: 
Drawing survival curves using 'ggplot2' (Version 0.4.9) [R package]. CRAN. 
https://doi.org/10.32614/CRAN.package.survminer 
7. Haider, H. (2019). MTLR: Survival Prediction with Multi-Task Logistic 
Regression (Version 0.2.1) [R package]. CRAN. 
https://doi.org/10.32614/CRAN.package.MTLR
8. Wickham H, François R, Henry L, Müller K, Vaughan D (2023). dplyr: A Grammar 
of Data Manipulation. R package version 1.1.4, 
https://github.com/tidyverse/dplyr, https://dplyr.tidyverse.org.
9. Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. 
Springer-Verlag New York. ISBN 978-3-319-24277-4, 
https://ggplot2.tidyverse.org.
10. Zhou, H., Cheng, X., Wang, S., Zou, Y., & Wang, H. (2022). SurvMetrics: 
Predictive Evaluation Metrics in Survival Analysis (Version 0.5.0) [R package]. 
CRAN. https://doi.org/10.32614/CRAN.package.SurvMetrics
11. Simon N, Friedman J, Tibshirani R, Hastie T (2011). “Regularization Paths 
for Cox's Proportional Hazards Model via Coordinate Descent.” Journal of 
Statistical Software, 39(5), 1–13. doi:10.18637/jss.v039.i05.
12. Gerds TA (2023). pec: Prediction Error Curves for Risk Prediction Models in 
Survival Analysis. R package version 2023.04.12, 
https://CRAN.R-project.org/package=pec.
13. Heagerty, P. J., & Saha-Chaudhuri, P. (2022). survivalROC: Time-Dependent 
ROC Curve Estimation from Censored Survival Data (Version 1.0.3.1) [R package]. 
CRAN. https://doi.org/10.32614/CRAN.package.survivalROC
14. Harrell, F. E. Jr. (2024). rms: Regression Modeling Strategies 
(Version 6.8-1) [R package]. CRAN. https://doi.org/10.32614/CRAN.package.rms
15. Sing T, Sander O, Beerenwinkel N, Lengauer T (2005). “ROCR: visualizing 
classifier performance in R.” Bioinformatics, 21(20), 7881. 
http://rocr.bioinf.mpi-sb.mpg.de.
16. Bates, D., Maechler, M., Jagan, M., Davis, T. A., Karypis, G., Riedy, J., 
Oehlschlägel, J., & R Core Team. (2024). Matrix: Sparse and Dense Matrix 
Classes and Methods (Version 1.7-0) [R package]. CRAN. 
https://doi.org/10.32614/CRAN.package.Matrix
17. Harrell, F. E. Jr., & Dupont, C. (2024). Hmisc: Harrell Miscellaneous 
(Version 5.1-3) [R package]. CRAN. https://doi.org/10.32614/CRAN.package.Hmisc
Wickham H (2007). “Reshaping Data with the reshape Package.” Journal of 
Statistical Software, 21(12), 1–20. http://www.jstatsoft.org/v21/i12/.
18. Das P, Roychowdhury A, Das S, Roychoudhury S, Tripathy S (2020). 
“sigFeature: Novel Significant Feature Selection Method for Classification of 
Gene Expression Data Using Support Vector Machine and t Statistic.” Frontiers 
in genetics, 11, 247. doi:10.3389/fgene.2020.00247, 
https://www.frontiersin.org/article/10.3389/fgene.2020.00247.
