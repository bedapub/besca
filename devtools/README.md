Scripts and tools used by BESCA developers
===
Jitao David Zhang

## Test besca package in a virtual environment

The following commands can be used to test `besca` in a virtual environment.

```
## load python 3.8
ml purge
ml load Python/3.8.2-GCCcore-9.3.0
## create and activate a virtual environment
python -m venv test_besca
source test_besca/bin/activate
test_besca/bin/python -m pip install --upgrade pip
## clone besca and install it into the virtual environment
git clone --branch scanpy1.8.2 git@github.com:Accio/besca.git
cd test_besca/bin
./python -m pip install ../../besca
## run notebooks
cd ../../besca
bash devtools/run_workbooks.bash
## capture the version that we use
../../test_besca/bin/pip freeze > devtools/requirements.txt
```

You can use `pip install -r devtools/requirements.txt` to get an identical environment in which the test was done.
