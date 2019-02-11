################################################################################
# astro-accelerate.sh
#
# Description: A shell script that runs astro-accelerate
# Usage:       Please run the script from the same directory as the executable.
# Usage:       The input is the name of the input_file to use.
# Usage:       The second input is the directory where the output should go.
# Notice:      N/A.
################################################################################

#!/usr/bin/env bash

echo "=== Running astro-accelerate.sh script ==="
# Parse the parameter that contains the path to the input_file
input_file=$(readlink -f $1)

# Name of the input_file, without a path and without file extension
job_name=$(basename "$1" .txt)

# The current directory is assumed to contain the executable
executable_directory=${ASTRO_ACCELERATE_EXECUTABLE_PATH}

# Job output will be stored within a subfolder of the specified directory
output_dir=${2}/output/${job_name}

# Create the job directory
mkdir -p ${output_dir}

# The executable output defaults to the current directory
# So, switch to the output directory
cd ${output_dir}

# Remove existing output files, if they exist
rm -f analysed* 
rm -f acc*
rm -f global*
rm -f fourier*
rm -f harmonic*
rm -f candidate*
rm -f peak*

# Run and time the executable
echo "The executable path is " ${executable_directory}
echo "The input_file path is " ${input_file}
echo "The job name is " ${job_name}
echo "The output directory path is " ${output_dir}

# Go to output_dir so that output is stored there
cd ${output_dir}

# Run the executable
time ${executable_directory}/astro-accelerate ${input_file}

# Test if output file exists, and combine into a single file
# If not found, do nothing, and write the return code of ls to /dev/null
if ls analysed* 1> /dev/null 2>&1; then
    cat analysed* > global_analysed_frb.dat
fi

if ls fourier-* 1> /dev/null 2>&1; then
    cat fourier-* > global_periods.dat
fi

if ls fourier_inter* 1> /dev/null 2>&1; then
    cat fourier_inter* > global_interbin.dat
fi

if ls harmo* 1> /dev/null 2>&1; then
    cat harmo* > global_harmonics.dat
fi

if ls candidate* 1> /dev/null 2>&1; then
    cat candidate* > global_candidates.dat
fi

if ls peak* 1> /dev/null 2>&1; then
    cat peak* > global_peaks.dat
fi

echo "Finished. Output is located in "${output_dir}
