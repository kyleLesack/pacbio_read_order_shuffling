#!/usr/bin/env sh
input=$(realpath $1)
outdir=$(realpath $2)
#output=$(realpath $3)
mkdir -p ${outdir}
tail -q -n 1 ${input} 
