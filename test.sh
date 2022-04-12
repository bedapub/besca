#!/bin/bash

python -V
conda init bash
source ~/.bashrc
conda activate besca_dev
ls
pytest --doctest-modules -W ignore::PendingDeprecationWarning