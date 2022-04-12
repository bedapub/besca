#!/bin/bash

whereis conda
python -V
conda init bash
source ~/.bashrc
source ~/anaconda3/etc/profile.d/conda.sh
conda activate besca_dev
ls
pytest --doctest-modules -W ignore::PendingDeprecationWarning