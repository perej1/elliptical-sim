#!/bin/bash

# Remove old figures
rm figures/*

# Create arguments
Rscript gen-args.R
n=$(sed -n "$=" arg.txt)

# Perform batch of simulation rounds
for ((i = 1; i <= $n; i++))
do
  echo "Plot number:" $i"/"$n
  arg=$(sed -n "${i} p" arg.txt)
  Rscript ${arg}
done
rm arg.txt

