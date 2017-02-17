#!/bin/bash

rm rangeParams.h

for range in {0,1,2,3,4,5,6,7}
do
	paramString=$(cat range"$range"Stats.txt | tail -1 | awk -F " " '{print $3}')

	echo $paramString | sed 's/[^0-9^_]//g' | awk -v var="$range" -F "_" '{print "#define UNROLLS_" var" "$1"\n#define SNUMREG_" var " "$2"\n#define SDIVINT_" var " "$3"\n#define SDIVINDM_" var " "$4"\n#define SFDIVINDM_" var " "$4".0f"}' >> rangeParams.h
done
