#!/bin/bash

rm kernel_params.h

maxRange=${1-7}

for range in `seq 0 $maxRange`
do
	paramString=$(cat profile_optimum/range"$range"Stats.txt | tail -1 | awk -F " " '{print $3}')

	echo $paramString | sed 's/[^0-9^_]//g' | awk -v var="$range" -F "_" '{print "#define UNROLLS_" var" "$1"\n#define SNUMREG_" var " "$2"\n#define SDIVINT_" var " "$3"\n#define SDIVINDM_" var " "$4"\n#define SFDIVINDM_" var " "$4".0f"}' >> kernel_params.h
done

cp kernel_params.h ../lib/AstroAccelerate/kernel_params.h
