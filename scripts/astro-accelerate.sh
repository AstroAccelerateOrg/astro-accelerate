#!/bin/bash

job_name=$(basename "$1" .txt)

# find out where this script is running from
scripts_dir=$(dirname "$0");
cd $scripts_dir;
scripts_dir=$PWD

# define where to place output and find out executables
project_top=$(dirname "$scripts_dir")
output_dir="$project_top/output"
job_dir=$output_dir/$job_name
build_dir=$project_top/build

input_file=$(readlink -f $1)

mkdir -p $job_dir

cd $job_dir

rm -f analysed* 
rm -f acc*
rm -f global*
rm -f fourier*
rm -f harmonic*
rm -f candidate*
rm -f peak*

time $build_dir/dedisperse-gpu $input_file

cat analysed* > global_analysed_frb.dat
cat fourier-* > global_periods.dat
cat fourier_inter* > global_interbin.dat
cat harmo* > global_harmonics.dat
cat candidate* > global_candidates.dat
cat peak* > global_peaks.dat

echo "Written output to $job_dir"
