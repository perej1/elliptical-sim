#!/bin/bash

# Create arguments
Rscript gen-args.R

n=$(sed -n "$=" args/arg-gen-sample.txt)
echo "Create samples"
for ((i = 1; i <= $n; i++))
do
  echo "Scenario number:" $i"/"$n
  arg=$(sed -n "${i} p" args/arg-gen-sample.txt)
  Rscript ${arg}
done

n=$(sed -n "$=" args/arg-estimate.txt)
echo "Compute extreme quantile region estimates"
for ((i = 1; i <= $n; i++))
do
  echo "Scenario number:" $i"/"$n
  arg=$(sed -n "${i} p" args/arg-estimate.txt)
  Rscript ${arg}
done

n=$(sed -n "$=" args/arg-plot.txt)
echo "Plotting"
for ((i = 1; i <= $n; i++))
do
  echo "Scenario number:" $i"/"$n
  arg=$(sed -n "${i} p" args/arg-plot.txt)
  Rscript ${arg}
done
