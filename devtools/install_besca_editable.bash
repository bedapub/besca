#!/bin/bash

kernel=besca_dev

## pre-conditions:
## 1. you are in the root path of besca, and 
## 2. uv is installed (https://docs.astral.sh/uv/)
uv sync
uv run python -m ipykernel install --user --name "$kernel" --display-name "$kernel"

