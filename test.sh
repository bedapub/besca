#!/bin/bash

whereis conda
python -V
conda init bash
export PATH=/opt/conda/bin:$PATH
in ~/.bashrc
source ~/.bashrc
source ~/anaconda3/etc/profile.d/conda.sh
conda activate besca_dev
ls
pytest --doctest-modules -W ignore::PendingDeprecationWarning