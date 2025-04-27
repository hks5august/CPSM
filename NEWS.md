# CPSM 0.99.2 (Unreleased)
- Fixed indentation issues in R source files.
- Removed redundant `Maintainer` field from the DESCRIPTION file.
# CPSM 0.99.3
## Updates
- Updated `DESCRIPTION` file to include `file LICENSE` in the `License` field.
- Removed indentation inconsistencies in all `R` scripts for better readability.
- Bumped the version from 0.99.2 to 0.99.3 as per Bioconductor guidelines.

## Fixes
- Addressed `BiocCheck` notes on coding practices.
# CPSM 0.99.4 (Unreleased)
- Removed indentation inconsistencies in all `R` scripts for better readability.
- Bumped the version from 0.99.3 to 0.99.4 as per Bioconductor guidelines.

# CPSM 0.99.5 (Unreleased)
- Updated the function previously dependent on the SurvMetrics package. The computation is now implemented using the survival and pec packages, removing the dependency on SurvMetrics.
- Resolved issues with the NAMESPACE file:
- Removed exportPattern("^[[:alpha:]]+") from the NAMESPACE file.
- Bumped the version from 0.99.4 to 0.99.5 as per Bioconductor guidelines

# CPSM 0.99.6 (Unreleased)
- Updated vigenette file "CPSM.Rmd" included YAML header of the vignette
- corrected input in MTLR_pred_model_f.R file
- updated .gitignore file
- Bumped the version from 0.99.5 to 0.99.6 as per Bioconductor guidelines

# CPSM 0.99.7 (Unreleased)
- Updated vigenette file "CPSM.Rmd", updated installation instruction as per Bioconductor guidelines
- updated data_process_f.R
- Updated DESCRIPTION  file to add author and version update
- Bumped the version from 0.99.6 to 0.99.7 as per Bioconductor guidelines

# CPSM 1.0.0 (Released)
## First Official release of the **CPSM** package on Bioconductor.
- Includes core functionality for:
- Data preprocessing, and data normalization
- Feature selection for survival modeling
- Model prediction and evaluation
- Visualization using survival cureves plots, mean-median barplot, nomograms.
## Documentation
- Initial vignette `CPSM.Rmd` added, detailing step-by-step usage of the package.
- Includes installation instructions, use-case workflows, and examples.
- Ensured compliance with Bioconductor submission guidelines.
- Bumped version from **0.99.7** (development) to **1.0.0** (release).

# CPSM 1.1.0 
## New Features
- **`predict_survival_risk_group_f()`**:
  - Added a new function to build Random Forest-based survival risk group classifiers using selected molecular features.
  - Automatically selects the best-performing model based on prediction accuracy across various `ntree` values (10 to 1000).
  - Outputs include predicted risk groups, prediction probabilities, misclassification summary, best model, and tree-wise performance matrix.
- **`km_overlay_plot_f()`**:
  - Introduced a new function to generate Kaplan-Meier survival plots overlayed with individual test sample survival curves.
  - Useful for visualizing predicted test sample risk in the context of population-level survival groups.
  - Includes annotated output with sample ID, predicted group, and probability.
## Updated data files, added new files
## Improvement
- Updated MTLR_pred_model_f function to compute MAE and improve the structure of function. 
## Documentation
- Updated Description file
- Updated vignette file `CPSM.Rmd`:
  - Included usage details and examples for new functions.
  - Enhanced step-by-step explanation of prediction and visualization workflows.
## Version update
- Bumped version from **1.0.0** to **1.1.0**

# CPSM 1.1.1
##Improvement
- Updated MTLR_pred_model_f function to make it simple
## Version update
- Bumped version from **1.1.0** to **1.1.1**.

