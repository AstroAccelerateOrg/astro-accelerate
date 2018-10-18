################################################################################
# profiling.sh
#
# Description: A brute force optimiser for astro-accelerate.
# Usage:       Please run the script from the top-level directory.
# Notice:      Please do not commit the optimiser output to the repository.
################################################################################

#!/usr/bin/env bash
inFile="$1"

rm -rf profile_results

mkdir profile_results

for unroll in {4,8,16,32}
do
    for acc in {6,8,10,12,14,16}
    do
	for divint in {8,10,12,14,16,18}
	do
	    for divindm in {20,25,32,40,50,60}
	    do		
		pwd
		cat include/params.hpp > ./params.txt
		rm include/params.hpp
		echo "#define UNROLLS $unroll" >> ./params.txt
		echo "#define SNUMREG $acc" >> ./params.txt
		echo "#define SDIVINT $divint" >> ./params.txt
		echo "#define SDIVINDM $divindm" >> ./params.txt
		echo "#define SFDIVINDM $divindm.0f" >> ./params.txt
		
		mv params.txt include/params.hpp
		
		make clean
		
		regcount=$(make -j 16 2>&1 | grep -A2 shared_dedisperse_kernel | tail -1 | awk -F" " '{print $5}')
		
		cp include/params.hpp profile_results/u"$unroll"_a"$acc"_t"$divint"_dm"$divindm"_r"$regcount".h
		
		./astro-accelerate.sh $inFile > profile_results/u"$unroll"_a"$acc"_t"$divint"_dm"$divindm"_r"$regcount".dat
		
		echo "unrolls: $unroll	acc: $acc    divint: $divint    divindm: $divindm    reg: $regcount"
	    done
	done
    done
done

optimum=$(grep "Real" profile_results/* | awk -F" " '{print $4" "$1}' | sort -n | tail -1 | awk -F" " '{print $2}' | awk -F"." '{print $1".h"}')

cp $optimum include/params.hpp
pwd
make clean
regcount=$(make -j 16 2>&1 | grep -A2 shared_dedisperse_kernel | tail -1 | awk -F" " '{print $5}')
cd scripts/

echo "FINISED OPTIMISATION"
