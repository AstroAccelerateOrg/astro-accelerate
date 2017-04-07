#!/bin/bash

maxThreads=1024

rm -rf profile_results

mkdir profile_results

for unroll in {16,32,64}
do
	for acc in {6,8,10}
	do
		for divint in {6,8,10}
		do
			for divindm in {32,40}
			do
				threads=$(echo "$divint * $divindm" | bc)
				if [ $threads -le $maxThreads ]; then
					cd ../lib
					pwd
					rm AstroAccelerate/params.h
					echo "#define ACCMAX 350" > ./params.txt
					echo "#define ACCSTEP 11" >> ./params.txt
					echo "#define UNROLLS $unroll" >> ./params.txt
					echo "#define SNUMREG $acc" >> ./params.txt
					echo "#define SDIVINT $divint" >> ./params.txt
					echo "#define SDIVINDM $divindm" >> ./params.txt
					echo "#define SFDIVINDM $divindm.0f" >> ./params.txt
					echo "#define CARD 1" >> ./params.txt
					echo "#define NOPSSHIFT 5" >> ./params.txt
					echo "#define NOPSLOOP 3" >> ./params.txt
					echo "#define NDATAPERLOOP 1" >> ./params.txt
					echo "#define BINDIVINT 6" >> ./params.txt
					echo "#define BINDIVINF 32" >> ./params.txt
					echo "#define CT 256" >> ./params.txt
					echo "#define CF 2" >> ./params.txt
					echo "#define NOPS 4.0" >> ./params.txt
					echo "#define STATST 128" >> ./params.txt
					echo "#define STATSLOOP 8" >> ./params.txt

					mv params.txt AstroAccelerate/params.h
	
					make clean
			
					regcount=$(make -j 16 2>&1 | grep -A2 shared_dedisperse_kernel | tail -1 | awk -F" " '{print $5}')

					cd ../scripts/

					./astro-accelerate.sh ../input_files/one_bit.txt > profile_results/u"$unroll"_a"$acc"_t"$divint"_dm"$divindm"_r"$regcount".dat
					cp ../lib/AstroAccelerate/params.h profile_results/u"$unroll"_a"$acc"_t"$divint"_dm"$divindm"_r"$regcount".h
			
					echo "unrolls: $unroll	acc: $acc    divint: $divint    divindm: $divindm    reg: $regcount"
			
				fi
			done
		done
	done
done

optimum=$(grep "Real" * | awk -F" " '{print $4" "$1}' | sort -n | tail -1 | awk -F" " '{print $2}' | awk -F"." '{print $1".h"}')

cp profile_results/$optimum ../lib/AstroAccelerate/params.h
cd ../lib
pwd
make clean
regcount=$(make -j 16 2>&1 | grep -A2 shared_dedisperse_kernel | tail -1 | awk -F" " '{print $5}')
cd ../scripts/


