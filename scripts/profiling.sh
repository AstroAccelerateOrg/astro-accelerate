#!/bin/bash

maxThreads=1024
numTrials=5

rm -rf profile_optimum
mkdir profile_optimum

maxRange=${1-7}

for range in `seq 0 $maxRange` 
do
	rm profile_optimum/range${range}Stats.txt
	rm -rf profile_results
	mkdir profile_results
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
								rm headers/params.h

								echo "#define UNROLLS $unroll" >> ./params.txt
								echo "#define SNUMREG $acc" >> ./params.txt
								echo "#define SDIVINT $divint" >> ./params.txt
								echo "#define SDIVINDM $divindm" >> ./params.txt
								echo "#define SFDIVINDM $divindm.0f" >> ./params.txt
								echo "#define FILTER_OUT_RANGES 1" >> ./params.txt
								echo "#define RANGE_TO_KEEP $range" >> ./params.txt

								cat header_template.h >> ./params.txt 
								mv params.txt headers/params.h
				
								make clean
						
								regcount=$(make -j 16 2>&1 | grep -A2 shared_dedisperse_kernel | tail -1 | awk -F" " '{print $5}')
								paramString=u"$unroll"_a"$acc"_t"$divint"_dm"$divindm"_r"$regcount"
								cd ../scripts/
							fi

							./astro-accelerate.sh ../input_files/ska_ania.txt > profile_results/"$trialString"_"$paramString".dat

							if [ $trial -eq 1 ]; then
								cp ../lib/headers/params.h profile_results/"$paramString".h
							fi

							echo "unrolls: $unroll  acc: $acc    divint: $divint    divindm: $divindm    reg: $regcount"
							speedupThisTrial=$(grep "SPEEDUP FACTOR" ./profile_results/"$trialString"_"$paramString".dat | awk -F" " '{print $6}')
							totalSpeedup=$(echo $speedupThisTrial+$totalSpeedup | bc)
						done

					averageSpeedup=$(echo "scale=5;$totalSpeedup/$numTrials" | bc)
					bestTrial=$(grep "SPEEDUP FACTOR" ./profile_results/*"$paramString"* | awk -F" " '{print $6}' | sort -n | tail -1)
					worstTrial=$(grep "SPEEDUP FACTOR" ./profile_results/*"$paramString"* | awk -F" " '{print $6}' | sort -n | head -1)
					echo params: $paramString avg: $averageSpeedup max: $bestTrial min: $worstTrial >> profile_optimum/range${range}Stats.txt
					
					fi
				done
			done
		done
	done

	optimumParams=$(cat profile_optimum/range"$range"Stats.txt | awk -F" " '{print $4" "$2}' | sort -n | tail -1 | awk -F" " '{print $2}')
	optimumParamString=$(grep "$optimumParams" profile_optimum/range"$range"Stats.txt)
	echo BEST: $optimumParamString >> profile_optimum/range${range}Stats.txt

	cp ./profile_results/$optimumParams.h ../lib/headers/params.h
	cp ./profile_results/$optimumParams.h profile_optimum/params_$range.h

	for trial in $(seq "$numTrials")
	do
		cp ./profile_results/trial${trial}_$optimumParams.dat profile_optimum/trial${trial}_params_$range.dat
	done

done

#cd ../lib
#pwd
#make clean
#regcount=$(make -j 16 2>&1 | grep -A2 shared_dedisperse_kernel | tail -1 | awk -F" " '{print $5}')
#cd ../scripts/


