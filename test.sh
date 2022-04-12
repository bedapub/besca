#!/bin/bash

python -V
conda init bash
source ~/.bashrc
source activate besca_dev
ls
pytest --doctest-modules -W ignore::PendingDeprecationWarning