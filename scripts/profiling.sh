#!/bin/bash

maxThreads=1024

rm -rf profile_results

mkdir profile_results

for acc in {8,10,12}
do
	for divint in {8,10,12}
	do
		for divindm in {72,75}
		do
			threads=$(echo "$divint * $divindm" | bc)
			if [ $threads -le $maxThreads ]; then
				cd ../lib
				pwd
				rm params.h
				echo "#define NUMREG $acc" > ./params.txt
				echo "#define DIVINT $divint" >> ./params.txt
				echo "#define DIVINDM $divindm" >> ./params.txt
				echo "#define FDIVINDM $divindm.0f" >> ./params.txt
				echo "#define SNUMREG $acc" >> ./params.txt
				echo "#define SDIVINT $divint" >> ./params.txt
				echo "#define SDIVINDM $divindm" >> ./params.txt
				echo "#define SFDIVINDM $divindm.0f" >> ./params.txt
				echo "#define CARD 0" >> ./params.txt
				echo "#define NOPSSHIFT 5" >> ./params.txt
				echo "#define NOPSLOOP 3" >> ./params.txt
				echo "#define NDATAPERLOOP 1" >> ./params.txt
				echo "#define BINDIVINT 6" >> ./params.txt
				echo "#define BINDIVINF 32" >> ./params.txt
				echo "#define CT 256" >> ./params.txt
				echo "#define CF 2" >> ./params.txt
				echo "#define NOPS 4.0" >> ./params.txt

				mv params.txt AstroAccelerate/params.h

				make clean
			
				regcount=$(make -j 8 gpu=sm_35 2>&1 | grep -A2 shared_dedisperse_kernel | tail -1 | awk -F" " '{print $5}')

#				mv dedisperse-gpu ..

				cd ../scripts/

				./astro-accelerate.sh ../input_files/ska_bin_8.txt > profile_results/a"$acc"_t"$divint"_dm"$divindm"_r"$regcount".dat
			
				echo "acc: $acc    divint: $divint    divindm: $divindm    reg: $regcount"
			
			fi
		done
	done
done
