# Environmmental drivers of early-onset colorectal cancer in the United States

This folder includes code for the statistical analysis in Chen J, Terry MB, Piero D, Chin H, Hu J, Yang W. Environmental drivers of the rising incidence of early-onset colorectal cancer in the United States. 2023. (Under review by The International Journal of Cancer)

- mvm_compare.R (Testing models without interaction effects; ie, Step 1 of the analysis)
- mvm_interact.R (Testing models with interaction effects; ie, Step 2 of the analysis)
- mvm_validate.R (5-fold cross-validation)
- mvm_counterfactual.R (Counterfactual modeling)

Of note, all functions needed to run the above analysis are in utils.R, crc_mvm.R, and crc_mvm_other.R.
