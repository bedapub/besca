#!/bin/bash

kernel=besca_dev

## pre-conditions:
## 1. you are in the root path of besca, and 
## 2. the conda/venv environment for besca development is activated
pip uninstall besca -y
pip install .
python -m ipykernel install --user --name "$kernel" --display-name "$kernel"

