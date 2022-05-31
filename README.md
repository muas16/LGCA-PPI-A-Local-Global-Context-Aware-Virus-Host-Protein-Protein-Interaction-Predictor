# LGCA-PPI-A-Local-Global-Context-Aware-Virus-Host-Protein-Protein-Interaction-Predictor

## Requirements
* Python 3
* Scikit-learn
* Pandas
* Numpy 1.14.3+

### Training

For performing experiment on Denovo or Yong_et_al Datasets, Run:

    python main.py


For performing experiment on Berman Dataset, Run:

    python main_kfold.py


### Parameter Settings:
* Set the **encoding_type** parameter such as ["Global", "Local","Local_Global"] (***NOTE***: you can also set one of them)
* set parameter "dir" (Path of folder where all results and encoding will be saved)
* set parameter "datasets" (Path of datasets files)
*In case of Berman Dataset, also set kfold parameter

## Datssets
Download the Datasets from :
https://sds_genetic_analysis.opendfki.de/HVI/Download/
