################################################################################
# unit_test_and_results.sh
#
# Description: Runs astro-accelerate over all files in a directory.
# Usage:       The input is the path to search for input_files.
# Notice:      N/A.
################################################################################

#!/usr/bin/env bash 

if [ $# -eq 0 ]
	then
    		echo "No arguments supplied."
		echo "You must supply a directory location."
		exit
fi

# Directory of list of input_files to use.
INPUT_FILES_PATH=${1}

rm -rf ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/*
rm -rf ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/tests/*
rm -rf ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/input_files/tests

mkdir -p ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/input_files/tests
mkdir -p ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/tests


for i in $(ls $INPUT_FILES_PATH/*.fil)
do
	echo $i
	cat ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/input_files/header > input.txt
	echo "file $i" >> input.txt
	j=$(basename $i)
	mv input.txt ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/input_files/tests/$j.txt
done

for i in $(ls ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/input_files/tests/*)
do
	echo $i
	${ASTRO_ACCELERATE_SCRIPTS_PATH}/astro-accelerate.sh $i
done

for i in $(ls ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/*/gl*frb*)
do
	j=$(echo $i | awk -F"/" '{print $3}'| awk -F"." '{print $1}')
	echo $j
	echo "set terminal png" > $j.plt
	echo "set output \"$j.png\"" >> $j.plt
	echo "set title \"$j\"" >> $j.plt
	echo "splot \"./$i\" binary format=\"%float%float%float%float\" u 1:2:3 palette" >> $j.plt
	mv $j.plt ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/tests/$j.plt
done

for i in $(ls ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/tests/*.plt)
do
	gnuplot $i
done

convert -delay 125 -loop 0 *.png ani.gif

mv *.png ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/tests/
mv *.gif ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/tests/

google-chrome-stable ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/tests/ani.gif &
