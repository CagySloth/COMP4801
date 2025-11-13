#!/bin/bash

# Wrapper for dataset/simulate.py
# Usage: ./scripts/simulate.sh --ploidy 2 --num-variants 1000 --num-reads 5000 ...

python dataset/simulate.py "$@"
