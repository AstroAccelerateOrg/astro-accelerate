#!/bin/bash

maxThreads=1024
numTrials=5

rm -rf profile_optimum_averaged
mkdir profile_optimum_averaged

rm -rf profile_results
mkdir profile_results

rm profile_optimum_averaged/stats.txt

unroll=16
for unroll in {16,32,64}
do
	for acc in {6,8,10,12,14}
	do
		for divint in {6,8,10,12,14}
		do
			for divindm in {32,40,50,60}
			do
				threads=$(echo "$divint * $divindm" | bc)
				if [ $threads -le $maxThreads ]; then
					totalSpeedup=0

					for trial in $(seq "$numTrials")
					do
						trialString=trial"$trial"
						if [ $trial -eq 1 ]; then
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
							echo "#define FILTER_OUT_RANGES 0" >> ./params.txt
							echo "#define RANGE_TO_KEEP 0" >> ./params.txt

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

							paramString=u"$unroll"_a"$acc"_t"$divint"_dm"$divindm"_r"$regcount"
							cd ../scripts/
						fi

						./astro-accelerate.sh ../input_files/ska_ania.txt > profile_results/"$trialString"_"$paramString".dat

						if [ $trial -eq 1 ]; then
							cp ../lib/AstroAccelerate/params.h profile_results/$paramString.h
						fi
						
						echo "unrolls: $unroll	acc: $acc    divint: $divint    divindm: $divindm    reg: $regcount"
						speedupThisTrial=$(grep "Real" ./profile_results/"$trialString"_"$paramString".dat | awk -F" " '{print $4}')
						totalSpeedup=$(echo $speedupThisTrial+$totalSpeedup | bc)

				
					done
				averageSpeedup=$(echo "scale=5;$totalSpeedup/$numTrials" | bc)
				bestTrial=$(grep "Real" ./profile_results/*"$paramString"* | awk -F" " '{print $4" "$1}' | sort -n | tail -1 | awk -F" " '{print $1}')
				worstTrial=$(grep "Real" ./profile_results/*"$paramString"* | awk -F" " '{print $4" "$1}' | sort -n | head -1 | awk -F" " '{print $1}')
				echo params: $paramString avg: $averageSpeedup max: $bestTrial min: $worstTrial >> profile_optimum_averaged/stats.txt

				fi
			done
		done
	done
done

#optimum=$(grep "Real" ./profile_results/* | awk -F" " '{print $4" "$1}' | sort -n | tail -1 | awk -F" " '{print $2}' | awk -F":" '{print $1}' | awk -F"." '{print "."$2}')
#echo "optimum: "$optimum

optimumParams=$(cat profile_optimum_averaged/stats.txt | awk -F" " '{print $4" "$2}' | sort -n | tail -1 | awk -F" " '{print $2}')
optimumParamString=$(grep "$optimumParams" profile_optimum_averaged/stats.txt)
echo BEST: $optimumParamString >> profile_optimum_averaged/stats.txt

cp ./profile_results/$optimumParams.h ../lib/AstroAccelerate/params.h
cp ./profile_results/$optimumParams.h profile_optimum_averaged/params_averaged.h

for trial in $(seq "$numTrials")
do
	cp ./profile_results/trial${trial}_$optimumParams.dat profile_optimum_averaged/trial${trial}_params_averaged.dat
done


#cd ../lib
#pwd
#make clean
#regcount=$(make -j 16 2>&1 | grep -A2 shared_dedisperse_kernel | tail -1 | awk -F" " '{print $5}')
#cd ../scripts/


