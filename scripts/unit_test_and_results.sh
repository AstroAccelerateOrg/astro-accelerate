################################################################################
# unit_test_and_results.sh
#
# Description: Runs astro-accelerate over all files in a directory.
# Usage:       The input is the path to search for input_files.
# Notice:      This script requires ImageMagick and google-chrome-stable.
################################################################################

#!/usr/bin/env bash 

echo "=== Running unit_test_and_results.sh script ==="

if [ $# -eq 0 ]
then
    echo "No arguments supplied."
    echo "You must supply a directory location."
    exit
fi

# Directory of list of input_files to use.
FIL_FILES_PATH=${1}

# Delete previous output directories
rm -rf ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/*
rm -rf ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/tests/*
rm -rf ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/input_files/tests

# Create output directories if not already created
mkdir -p ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/input_files/tests
mkdir -p ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/tests

# Loop over all fil files in the folder and create input_files for each one
for i in $(ls $FIL_FILES_PATH/*.fil)
do
    echo $i
    cat ${ASTRO_ACCELERATE_REPOSITORY_PATH}/input_files/header > input.txt
    echo "file $i" >> input.txt
    j=$(basename $i)
    mv input.txt ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/input_files/tests/$j.txt
done

# Loop over all newly created input_files and run astro-accelerate.sh on each one
for i in $(ls ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/input_files/tests/*)
do
    echo $i
    ${ASTRO_ACCELERATE_SCRIPTS_PATH}/astro-accelerate.sh $i ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}
done

# astro-accelerate.sh has merged the output data into single data files.
# Loop over all output files for each job, 
for i in $(ls ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/*/gl*frb*)
do
    # The extraction of the base name is not safe, it relies on a pre-defined format for the path
    parentname="$(basename "$(dirname "${i}")")"
    filename=$(basename -- "$parentname")
    filename="${filename%.*}"
    echo "Filename is "${filename}
    echo "set terminal png" > ${filename}.plt
    echo "set output \"${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/tests/${filename}.png\"" >> ${filename}.plt
    echo "set title \"$filename\"" >> ${filename}.plt
    echo "splot \"$i\" binary format=\"%float%float%float%float\" u 1:2:3 palette" >> ${filename}.plt
    mv ${filename}.plt ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/tests/${filename}.plt
done

for i in $(ls ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/tests/*.plt)
do
    gnuplot ${i}
done

# ImageMagick must be available in the environment
convert -delay 125 -loop 0 ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/tests/*.png ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/tests/ani.gif

# Enable / disable depending on whether google-chrome-stable is available on the system
#google-chrome-stable ${ASTRO_ACCELERATE_SCRIPTS_UNIT_TESTS_OUTPUT_PATH}/output/tests/ani.gif &
