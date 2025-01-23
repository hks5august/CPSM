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
