#!/bin/bash

maxThreads=1024

rm -rf profile_optimum
mkdir profile_optimum

for range in {0,1,2,3,4,5,6,7}
#for range in {0,1}
#for range in {6,7}
	do
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
						echo "#define CARD 0" >> ./params.txt
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
						echo "#define FILTER_OUT_RANGES 1" >> ./params.txt
						echo "#define RANGE_TO_KEEP $range" >> ./params.txt

						echo "" >> ./params.txt
						echo "//Added by Karel Adamek" >> ./params.txt
						echo "#define WARP 32" >> ./params.txt
						echo "#define HALF_WARP 16" >> ./params.txt
						echo "#define MSD_ELEM_PER_THREAD 8" >> ./params.txt
						echo "#define MSD_WARPS_PER_BLOCK 16" >> ./params.txt
						echo "#define THR_ELEM_PER_THREAD 4" >> ./params.txt
						echo "#define THR_WARPS_PER_BLOCK 4" >> ./params.txt
						echo "#define PD_NTHREADS 512" >> ./params.txt
						echo "#define PD_NWINDOWS 2" >> ./params.txt
						echo "#define PD_MAXTAPS 16" >> ./params.txt
						echo "#define PD_SMEM_SIZE 1280" >> ./params.txt
						echo "#define PD_FIR_ACTIVE_WARPS 2" >> ./params.txt
						echo "#define PD_FIR_NWINDOWS 2" >> ./params.txt


						mv params.txt AstroAccelerate/params.h
		
						make clean
				
						regcount=$(make -j 16 2>&1 | grep -A2 shared_dedisperse_kernel | tail -1 | awk -F" " '{print $5}')

						cd ../scripts/

						./astro-accelerate.sh ../input_files/ska_ania.txt > profile_results/u"$unroll"_a"$acc"_t"$divint"_dm"$divindm"_r"$regcount".dat
						cp ../lib/AstroAccelerate/params.h profile_results/u"$unroll"_a"$acc"_t"$divint"_dm"$divindm"_r"$regcount".h
				
						echo "unrolls: $unroll	acc: $acc    divint: $divint    divindm: $divindm    reg: $regcount"
				
					fi
				done
			done
		done
	done

	#optimum=$(grep "Real" ./profile_results/* | awk -F" " '{print $4" "$1}' | sort -n | tail -1 | awk -F" " '{print $2}' | awk -F":" '{print $1".h"}')
	optimum=$(grep "Real" ./profile_results/* | awk -F" " '{print $4" "$1}' | sort -n | tail -1 | awk -F" " '{print $2}' | awk -F":" '{print $1}' | awk -F"." '{print "."$2}')

	echo "optimum: "$optimum

	cp $optimum.h ../lib/AstroAccelerate/params.h
	cp $optimum.h profile_optimum/params_$range.h
	cp $optimum.dat profile_optimum/params_$range.dat

done

#cd ../lib
#pwd
#make clean
#regcount=$(make -j 16 2>&1 | grep -A2 shared_dedisperse_kernel | tail -1 | awk -F" " '{print $5}')
#cd ../scripts/


