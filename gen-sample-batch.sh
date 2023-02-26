#!/bin/bash

# Create arguments
Rscript gen-sample-args.R
n=$(sed -n "$=" args/arg-gen-sample.txt)

# Perform batch of simulation rounds
for ((i = 1; i <= $n; i++))
do
  echo "Round number:" $i"/"$n
  arg=$(sed -n "${i} p" args/arg-gen-sample.txt)
  Rscript ${arg}
done

