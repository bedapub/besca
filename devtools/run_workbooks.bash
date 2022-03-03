#!/bin/bash

## the script should run without any error before major release

## pre-conditions:
## 1. you are in the root path of besca, and 
## 2. the conda/venv environment for besca development is activated
## 3. besca_dev kernel has been installed (see install_besca_anew_local.bash)

for notebook in "minimal_notebook.ipynb" "standard_workflow_besca2.ipynb" "celltype_annotation_besca.ipynb" "Signature_exports.ipynb" "Testing_Notebook.ipynb"; do
    echo Running notebook "$notebook"
    jupyter nbconvert --to notebook --execute workbooks/"$notebook"
    jupyter nbconvert --clear-output workbooks/"$notebook"
done
