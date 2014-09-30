#!/bin/bash


job_name=$(basename "$1" .txt)
project_top=$(readlink -f $( dirname $(dirname "$0")))
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

time $build_dir/dedisperse-gpu $input_file

cat analysed* > global_analysed_frb.dat
cat fourier-* > global_periods.dat
cat fourier_inter* > global_interbin.dat
cat harmo* > global_harmonics.dat
cat candidate* > global_candidates.dat

echo "Written output to $job_dir"
