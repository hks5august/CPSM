
 # CPSM: Cancer Patient Survival Model

## Introduction
The CPSM R-package is an advanced computational pipeline designed to predict the survival probability of cancer patients with precision and efficiency. This package automates critical steps in survival analysis, including data preprocessing, normalization, and splitting datasets into training and testing subsets. It identifies significant features using univariate survival analysis and LASSO COX-regression, and also generates LASSO-based prognostic index (PI) scores to optimize predictive performance. CPSM  builds robust survival prediction models using selected clinical, molecular, and integrated feature sets. To support data interpretation and clinical applications, the package includes powerful visualization tools, such as survival curves, bar plots for predicted mean and median survival times, and nomograms. Designed for multi-omics data, CPSM simplifies complex workflows, empowering researchers to uncover novel biomarkers and advance precision oncology.

![CPSM_workflow_figure](https://github.com/user-attachments/assets/1bc0fd51-06e5-4739-881d-9c49f312f5fc)

Figure: The workflow of the CPSM package represents different steps performed by various functions of the CPSM package.


# Follow the Steps to Install the CPSM package on your Local R system:
```r
#Step1: First Install remote package
install.packages("remotes") 
#load remotes package
library("remotes")
#Step2: install CPSM package
remotes::install_github("hks5august/CPSM", local = TRUE, , dependencies=TRUE)
#or use the following command
remotes::install_github("hks5august/CPSM", ref = "v1.0.0", dependencies=TRUE)
# Check if package get installed, load package
library("CPSM") 
```

## Load Packages
```{r, warning=FALSE, message=FALSE }
#Load CPSM packages
library(CPSM)
#Load other required packages
library(preprocessCore)
library(ggfortify)
library(survival)
library(survminer)
library(ggplot2)
library(MASS)
library(MTLR)
library(dplyr)
library(SurvMetrics)
library(pec)
library(glmnet)
library(reshape2)
library(rms)
library(Matrix)
library(Hmisc)
library(survivalROC)
library(ROCR)
```

```{r }
#set seed
set.seed(7)
```


## Input Data
Example Input data: "Example_TCGA_LGG_FPKM_data" is a tab separated file.  It contains Samples (184 LGG Cancer Samples) in the rows and Features in the columns. Gene Expression is available in terms of FPKM values in the data.
Features information: In the data there are 11 clinical + demographic, 4 types survival with time and event information  and 19,978 protein coding genes.
Clinical and demographic features: Clinical demographic features that are present in this example data include Age,  subtype,  gender,  race,  ajcc_pathologic_tumor_stage,  histological_type, histological_grade,  treatment_outcome_first_course, radiation_treatment_adjuvant,  sample_type,  type.
Types of Survival: 4 types of Survival include OS (overall survival), PFS (progression-free survival), DSS (disease-specific survival), DFS (Disease-free survival). In the data, column names OS, PFS, DSS and DFS represent event information, while  OS.time, PFS.time, DSS.time and DFS.time indicate survival time in days.

## Step 1- Data Processing 
This function converts OS time (in days) into months and then removes samples where OS/OS.time information is missing.
Here, we need to provide input data in tsv or txt  format. Further, we needs to define col_num (column number at which clinical/demographic and survival information ends,e.g. 20,  surv_time (name of column which contain survival time (in days) information, e.g. OS.time ) and output file name, e.g.  “New_data.txt”

```{r }
data(Example_TCGA_LGG_FPKM_data, package = "CPSM")
New_data <- data_process_f(Example_TCGA_LGG_FPKM_data, col_num=20, surv_time = "OS.time")
str(New_data[1:10])
```

After data processing, we got a new output file “New_data”, which contains 176 samples. Thus, data_process_f function removes 8 samples where OS/OS time information is missing. Besides, here is a new 21st column in the data with  column name “OS_month”  where OS time is available in months.


## Step 2 - Split Data into Training and Test Subset
Before proceeding further, we need to split our data into training and test subset for the purpose of feature selection and model development. Here, we need output from the previous step as an input ( which was “New_data.txt”). Next we need to define the fraction (e.g. 0.9) by which we want to split data into training and test. Thus, fraction=0.9 will split data into 90% training and 10% as test set. Besides, we also need to provide training and set output names (e.g. train_FPKM.txt,test_FPKM.txt )


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

After the train-test split, we got a two new outputs:  “train_FPKM”, “test_FPKM”, where, train_FPKM contains 158 samples and test_FPKM contains 18 samples. Thus, tr_test_f function splits data into 90:10 ratio. 

## Step 3 - Data Normalization
Next to select features and develop ML models, data must be normalized. Since, expression is available in terms of FPKM values. Thus,  `train_test_normalization_f` function will first convert FPKM value into log scale [log2(FPKM+1) followed by quantile normalization using the “preprocessCore” package. Here, training data will be used as a target matrix for quantile normalization.  Here, we need to provide training and test datasets (that we obtained from the  previous step of Train/Test Split). Further, we need to provide column number where clinical information ends (e.g. 21) in the input datasets. Besides, we also need to provide output files names (train_clin_data (which contains only Clinical information of training data), test_clin_data (which contains only Clinical information of training data), train_Normalized_data_clin_data (which contains Clinical information and normalized values of genes of training samples), test_Normalized_data_clin_data (which contains Clinical information and normalized values of genes of test samples).


## Step 3 - Data Normalization
```{r }

# Normalize the training and test data sets
data(train_FPKM, package = "CPSM")
data(test_FPKM, package = "CPSM")
Result_N_data <- train_test_normalization_f(train_data = train_FPKM,
                            test_data = test_FPKM,
                            col_num = 21)
# Access the Normalized train and test data
Train_Clin <- Result_N_data$Train_Clin
Test_Clin <- Result_N_data$Test_Clin
Train_Norm_data <- Result_N_data$Train_Norm_data
Test_Norm_data <- Result_N_data$Test_Norm_data
str(Train_Clin[1:10])
str(Train_Norm_data[1:10])
```

After, running the function, we obtained 4 outputs: Train_Clin - Contains only Clinical features, Test_Clin - contains only Clinical features of Test samples;  Train_Norm_data - Clinical features with normalized values of genes for training samples; Test_Norm_data - Clinical features with normalized values of genes for test samples.


## Step 4a - Prognostic Index (PI)  Score Calculation
Next to create a survival model, we  will  create a Prognostic Index (PI)  Score. PI score is calculated based on the expression of the features selected by the LASSO regression model and their beta coefficients. For instance, 5 features  (G1, G2, G3, G4, and G5 and their coefficient values are B1, B2, B3, B4, and B5, respectively) selected by the LASSO method. Then PI score will be computed as following:

PI score  = G1*B1 + G2*B2 + G3 * B3 + G4*B4+ G5*B5

Here, we need to provide Normalized training (Train_Norm_data) and test data (Test_Norm_data)as input data that we have obtained from the previous function “train_test_normalization_f”. Further, we need to provide col_num n column number at which clinical features ends (e.g. 21), nfolds  (number of folds  e.g. 5) for the LASSO regression method to select features. We implemented LASSO using the “glmnet” package.  Further, we need to provide surv_time (name of column containing survival time in months, e.g. OS_month) and surv_event (name of column containing survival event information, e.g. OS) information in the data. Besides, we also need to provide names and training and test output file names to store data containing LASSO genes and PI values.


```{r, warning=FALSE, message=FALSE, fig.width=7, fig.height=4 }
# Step 4 - Lasso PI Score
data(Train_Norm_data, package = "CPSM")
data(Test_Norm_data, package = "CPSM")
Result_PI <- Lasso_PI_scores_f(train_data = Train_Norm_data,
                  test_data = Test_Norm_data, 
                  nfolds=5, 
                  col_num=21, 
                  surv_time = "OS_month", 
                  surv_event = "OS")
Train_Lasso_key_variables <- Result_PI$Train_Lasso_key_variables
Train_PI_data <- Result_PI$Train_PI_data
Test_PI_data <- Result_PI$Test_PI_data
str(Train_PI_data[1:10])
str(Test_PI_data[1:10])
```

Thus, Lasso_PI_scores_f gave us following outputs:
1. Train_Lasso_key_variables: List of features selected by LASSO and their beta coefficient values
2. Train_Cox_Lasso_Regression_lamda_plot: Lasso Regression Lambda plot.
3. Train_PI_data: It contains expression of genes selected by LASSO and PI score in the last column for training samples.
4. Test_PI_data: It contains expression of genes selected by LASSO and PI score in the last column for test samples.

![lasso_plot](https://github.com/user-attachments/assets/51abf244-672c-46eb-b65e-36b054f81f8e)
Figure: Lasso Regression Lambda plot


## Step 4b - Univariate  Survival Significant Feature Selection
Besides PI score, with the “Univariate_sig_features_f” function of CPSM package, we can select significant (p-value <0.05) features based on univariate survival analysis. These features are selected based on their capability to stratify high-risk and low-risk survival groups using the cut off value of their median expression.  
Here, we need to provide Normalized training (Train_Norm_data.txt) and test data (Test_Norm_data.txt)as input data that we have obtained from the previous function “train_test_normalization_f”. Further, we need to provide a “col_num” (e.g 21)column number at which clinical features ends. Further, we need to provide surv_time (name of column containing survival time in months, e.g. OS_month) and surv_event (name of column containing survival event information, e.g. OS) information in the data. Besides, we also need to provide names and training and test output file names to store data containing expression of selected genes.


```{r, warning=FALSE, message=FALSE }
#Step 4b - Univariate  Survival Significant Feature Selection.
data(Train_Norm_data, package = "CPSM")
data(Test_Norm_data, package = "CPSM")
Result_Uni <- Univariate_sig_features_f(train_data = Train_Norm_data, 
                          test_data = Test_Norm_data, 
                          col_num=21, 
                          surv_time = "OS_month" , 
                          surv_event = "OS")
Univariate_Survival_Significant_genes_List <- Result_Uni$Univariate_Survival_Significant_genes_List
Train_Uni_sig_data <- Result_Uni$Train_Uni_sig_data
Test_Uni_sig_data <- Result_Uni$Test_Uni_sig_data
str(Univariate_Survival_Significant_genes_List[1:10])
```

Thus, Univariate_sig_features_f  gave us following outputs:
Univariate_Survival_Significant_genes_List: a table of univariate significant genes along with their corresponding coefficient values, HR value, P-values, C-Index values.
Train_Uni_sig_data: It contains expression of significant genes selected by univariate survival analysis for training samples.
Test_Uni_sig_data: It contains expression of significant genes selected by univariate survival analysis for test samples.

## Step 5 - Prediction model development for survival probability of patients
After selecting significant or key features using LASSO or Univariate survival analysis, next we want to develop an ML prediction model to predict survival probability of patients. MTLR_pred_model_f function of CPSM give us multiple options to develop models including Only Clinical features (Model_type=1), PI score (Model_type=2), PI Score + Clinical features (Model_type=3), Significant Univariate features (Model_type=4), Significant Univariate features Clinical features (Model_type=5) using MTLR package. Further, here, we were interested in developing a model based on PI score. Thus, we need to provide following inputs: (1) Training data with only clinical features, (2) Test data with only clinical features, (3) Model type (e.g. 2, since we want to develop model based on PI score), (4) Training data with PI score , (5) Test data with PI score, (6) Clin_Feature_List (e.g. Key_PI_list.txt), a list of features which will be  used to build model . Furthermore, we also need to provide surv_time (name of column containing survival time in months, e.g. OS_month) and surv_event (name of column containing survival event information, e.g. OS) information in the clinical data

### Model1 for only Clinical features
```{r, warning=FALSE, message=FALSE, error = TRUE }
data(Train_Clin, package = "CPSM")
data(Test_Clin, package = "CPSM")
data(Key_Clin_feature_list, package = "CPSM")
Result_Model_Type1 <- MTLR_pred_model_f(train_clin_data = Train_Clin, 
                      test_clin_data = Test_Clin, 
                      Model_type = 1, 
                      train_features_data = Train_Clin, 
                      test_features_data = Test_Clin, 
                      Clin_Feature_List = Key_Clin_feature_list, 
                      surv_time = "OS_month", 
                      surv_event = "OS")
survCurves_data <- Result_Model_Type1$survCurves_data
mean_median_survival_time_data <- Result_Model_Type1$mean_median_survival_time_data
survival_result_based_on_MTLR <- Result_Model_Type1$survival_result_based_on_MTLR
Error_mat_for_Model <- Result_Model_Type1$Error_mat_for_Model
```




### Model2 based on PI score
```{r, warning=FALSE, message=FALSE, error = TRUE}
data(Train_Clin, package = "CPSM")
data(Test_Clin, package = "CPSM")
data(Train_PI_data, package = "CPSM")
data(Test_PI_data, package = "CPSM")
data(Key_PI_list, package = "CPSM")
Result_Model_Type2 <- MTLR_pred_model_f(train_clin_data = Train_Clin, 
                      test_clin_data = Test_Clin, 
                      Model_type = 2, 
                      train_features_data = Train_PI_data , 
                      test_features_data = Test_PI_data , 
                      Clin_Feature_List = Key_PI_list, 
                      surv_time = "OS_month", 
                      surv_event = "OS")
survCurves_data <- Result_Model_Type2$survCurves_data
mean_median_survival_time_data <- Result_Model_Type2$mean_median_survival_time_data
survival_result_based_on_MTLR <- Result_Model_Type2$survival_result_based_on_MTLR
Error_mat_for_Model <- Result_Model_Type2$Error_mat_for_Model
```


### Model3  based  Clinical features + PI score
```{r, warning=FALSE, message=FALSE, error = TRUE}
data(Train_Clin, package = "CPSM")
data(Test_Clin, package = "CPSM")
data(Train_PI_data, package = "CPSM")
data(Test_PI_data, package = "CPSM")
data(Key_Clin_features_with_PI_list, package = "CPSM")
Result_Model_Type3 <- MTLR_pred_model_f(train_clin_data = Train_Clin, 
                      test_clin_data = Test_Clin, 
                      Model_type = 3, 
                      train_features_data = Train_PI_data,
                      test_features_data = Test_PI_data,
                      Clin_Feature_List = Key_Clin_features_with_PI_list, 
                      surv_time = "OS_month",  
                      surv_event = "OS")
survCurves_data <- Result_Model_Type3$survCurves_data
mean_median_survival_time_data <- Result_Model_Type3$mean_median_survival_time_data
survival_result_based_on_MTLR <- Result_Model_Type3$survival_result_based_on_MTLR
Error_mat_for_Model <- Result_Model_Type3$Error_mat_for_Model
```

### Model4 based on Univariate features (Genes + Clinical features)
```{r, warning=FALSE, message=FALSE, error = TRUE}
data(Train_Clin, package = "CPSM")
data(Test_Clin, package = "CPSM")
data(Train_Uni_sig_data, package = "CPSM")
data(Test_Uni_sig_data, package = "CPSM")
data(Key_univariate_features_with_Clin_list, package = "CPSM")
Result_Model_Type5 <- MTLR_pred_model_f(train_clin_data = Train_Clin, 
                      test_clin_data = Test_Clin, 
                      Model_type = 4, 
                      train_features_data = Train_Uni_sig_data, 
                      test_features_data = Test_Uni_sig_data,
                      Clin_Feature_List =Key_univariate_features_with_Clin_list, 
                      surv_time = "OS_month", 
                      surv_event = "OS")
survCurves_data <- Result_Model_Type5$survCurves_data
mean_median_survival_time_data <- Result_Model_Type5$mean_median_survival_time_data
survival_result_based_on_MTLR <- Result_Model_Type5$survival_result_based_on_MTLR
Error_mat_for_Model <- Result_Model_Type5$Error_mat_for_Model
```

After, implementing MTLR_pred_model_f function , we got following outputs:
1. Model_with_PI.RData : Model on training data
2. survCurves_data : Table containing predicted survival probability of each patient at different time points. This data can be further used to plot the survival curve of patients.
3. mean_median_survival_time_data : Table containing predicted mean and median survival time  of each patient in the test data. This data can be further used for bar plots.
4. Error_mat_for_Model : Table containing performance parameters obtained on test data based on prediction model. It contains IBS score (Integrated Brier Score) =0.192, C-Index =0.81.


## Step 6 - Survival curves/plots for individual patient
Next to visualize survival of patients, we will plot survival curve plots using the surv_curve_plots_f function  based on the data “survCurves_data ” that we obtained from the previous step after running the  MTLR_pred_model_f function. Further,  the surv_curve_plots_f function also allows highlighting a specific patient on the curve. Thus the function  needs only two inputs: 1) Surv_curve_data, (2) Sample ID of a specific patient (e.g. TCGA-DB-A4XF-01) that needs to be highlighted.


```{r, warning=FALSE, message=FALSE, error = TRUE , fig.width=3, fig.height=3}
#Create Survival curves/plots for individual patients
data(survCurves_data, package = "CPSM")
surv_curve_plots_f(Surv_curve_data = survCurves_data,
                   selected_sample = "TCGA-DB-A4XF-01")
```

Here, we obtained two output plots:
1. Survival curves for all patients in the test data with different colors 
2. Survival curves for all patients (in black) and highlighted patient (red) in the test data 

![survival_curve1](https://github.com/user-attachments/assets/8547ed51-b5f4-4061-98e7-d58417c52f59)
Figure: Survival curves for all patients in the test data.

![survival_curve2](https://github.com/user-attachments/assets/5d11ca81-3125-43b7-94f5-1de98826f31a)
Figure: Survival curves for all patients (in black/grey) and highlighted patient (red).

## Step 7 - Bar Plot for predicted  mean and median survival time of individual patients
Next, to visualize the predicted survival time of patients, we will plot the barplot for mean/median using “mean_median_surv_barplot_f” function based on the data that we obtained from step 5 after running the  MTLR_pred_model_f function. Further,  the mean_median_surv_barplot_f function also allows highlighting a specific patient on the curve. Thus the function  needs only two inputs: 1) surv_mean_med_data, (2) Sample ID of a specific patient (e.g. TCGA-DB-A4XF-01) that needs to be highlighted.


```{r, warning=FALSE, message=FALSE, error = TRUE }
data(mean_median_survival_time_data, package = "CPSM")
mean_median_surv_barplot_f(surv_mean_med_data = mean_median_survival_time_data, 
                           selected_sample = "TCGA-DB-A4XF-01")
```

Here, we obtained two output plots:
1. Barplot for all patients in the test data, where the red color bar represents mean survival and cyan/green color bar represents median survival time. 
2. Barplot for all patients with a highlighted patient (dashed black outline) in the test data. It shows  this patient has a predicted mean and median survival is 81.58 and 75.50 months.

![barplot1](https://github.com/user-attachments/assets/c44244f6-72df-4a4f-b237-4cbc7630b141)
Figure: Barplot for mean and median survival time of all patients in the data, where the red color bar represents mean survival and cyan/green color bar represents median survival time.

![barplot2](https://github.com/user-attachments/assets/87ce14f9-2f27-4410-8a97-a7918e81434f)
Figure: Barplot for mean and median survival time of all patients (grey color) with a highlighted patient (colored) in the test data.

## Development of Nomogram based on Key features
Next, the Nomogram_generate_f function of CPSM will allow users to generate a nomogram plot for their data (training data containing all samples) based on user-defined clinical and other features in their data. For instance, we will generate a nomogram based on 6 features (Age, gender, race, histological_type, sample_type, PI). Here, we will provide data containing all the features (Samples in rows and features in columns) (e.g. Train_Data_Nomogram_input) and a list of features (feature_list_for_Nomogram) based on which we want to generate a nomogram.  Further, we also need to provide surv_time (name of column containing survival time in months, e.g. OS_month) and surv_event (name of column containing survival event information, e.g. OS) information in the data.


```{r, warning=FALSE, message=FALSE, error = TRUE, fig.width=7, fig.height=6 }
data(Train_Data_Nomogram_input, package = "CPSM")
data(feature_list_for_Nomogram, package = "CPSM")
Result_Nomogram <- Nomogram_generate_f(data = Train_Data_Nomogram_input,  
                    Feature_List = feature_list_for_Nomogram, 
                    surv_time = "OS_month", 
                    surv_event = "OS")
C_index_mat <- Result_Nomogram$C_index_mat
```

Here, we will get a Nomogram based on features that we provide. This nomogram can predict Risk (Event risk, eg, Death), 1-year, 3-year, 5-year and 10 years survival of patients.

![Nomogram](https://github.com/user-attachments/assets/206f91cc-178d-4d14-a128-2a695ad95e0e)
## SessionInfo

As last part of this document, we call the function "sessionInfo()", 
which reports the version numbers of R and all the packages used 
in this session. It is good practice to always keep such a record 
as it will help to trace down what has happened in case that an R 
script ceases to work because the functions have been changed in 
a newer version of a package.


```{r}
sessionInfo()
```

## Citation Info
CPSM: R-package of an Automated Machine Learning Pipeline for Predicting the Survival Probability of Single Cancer Patient. 
(2024) Harpreet Kaur, Pijush Das, Kevin Camphausen, Uma T Shankavaram. bioRxiv 2024.11.14.623597; doi: https://doi.org/10.1101/2024.11.14.623597

## References


1. Kuhn, Max (2008). “Building Predictive Models in R Using the caret Package.” Journal of Statistical Software, 28(5), 1–26. doi:10.18637/jss.v028.i05, https://www.jstatsoft.org/index.php/jss/article/view/v028i05.
2. Bolstad B (2024). preprocessCore: A collection of pre-processing functions. R package version 1.66.0, https://github.com/bmbolstad/preprocessCore.
3. Horikoshi M, Tang Y (2018). ggfortify: Data Visualization Tools for Statistical Analysis Results. https://CRAN.R-project.org/package=ggfortify.
4. Therneau T (2024). A Package for Survival Analysis in R. R package version 3.7-0, https://CRAN.R-project.org/package=survival.
5. Terry M. Therneau, Patricia M. Grambsch (2000). Modeling Survival Data: Extending the Cox Model. Springer, New York. ISBN 0-387-98784-3.
6. Kassambara, A., Kosinski, M., Biecek, P., & Scheipl, F. (2021). survminer: Drawing survival curves using 'ggplot2' (Version 0.4.9) [R package]. CRAN. https://doi.org/10.32614/CRAN.package.survminer 
7. Haider, H. (2019). MTLR: Survival Prediction with Multi-Task Logistic Regression (Version 0.2.1) [R package]. CRAN. https://doi.org/10.32614/CRAN.package.MTLR
8. Wickham H, François R, Henry L, Müller K, Vaughan D (2023). dplyr: A Grammar of Data Manipulation. R package version 1.1.4, https://github.com/tidyverse/dplyr, https://dplyr.tidyverse.org.
9. Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org.
10. Zhou, H., Cheng, X., Wang, S., Zou, Y., & Wang, H. (2022). SurvMetrics: Predictive Evaluation Metrics in Survival Analysis (Version 0.5.0) [R package]. CRAN. https://doi.org/10.32614/CRAN.package.SurvMetrics
11. Simon N, Friedman J, Tibshirani R, Hastie T (2011). “Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent.” Journal of Statistical Software, 39(5), 1–13. doi:10.18637/jss.v039.i05.
12. Gerds TA (2023). pec: Prediction Error Curves for Risk Prediction Models in Survival Analysis. R package version 2023.04.12, https://CRAN.R-project.org/package=pec.
13. Heagerty, P. J., & Saha-Chaudhuri, P. (2022). survivalROC: Time-Dependent ROC Curve Estimation from Censored Survival Data (Version 1.0.3.1) [R package]. CRAN. https://doi.org/10.32614/CRAN.package.survivalROC
14. Harrell, F. E. Jr. (2024). rms: Regression Modeling Strategies (Version 6.8-1) [R package]. CRAN. https://doi.org/10.32614/CRAN.package.rms
15. Sing T, Sander O, Beerenwinkel N, Lengauer T (2005). “ROCR: visualizing classifier performance in R.” Bioinformatics, 21(20), 7881. http://rocr.bioinf.mpi-sb.mpg.de.
16. Bates, D., Maechler, M., Jagan, M., Davis, T. A., Karypis, G., Riedy, J., Oehlschlägel, J., & R Core Team. (2024). Matrix: Sparse and Dense Matrix Classes and Methods (Version 1.7-0) [R package]. CRAN. https://doi.org/10.32614/CRAN.package.Matrix
17. Harrell, F. E. Jr., & Dupont, C. (2024). Hmisc: Harrell Miscellaneous (Version 5.1-3) [R package]. CRAN. https://doi.org/10.32614/CRAN.package.Hmisc
Wickham H (2007). “Reshaping Data with the reshape Package.” Journal of Statistical Software, 21(12), 1–20. http://www.jstatsoft.org/v21/i12/.
