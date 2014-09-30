#!/bin/bash

rm analysed*
rm acc*
rm global*
rm fourier*
rm harmonic*
rm candidate*

time ./dedisperse-gpu $1

cat analysed* > global_analysed_frb.dat
cat fourier-* > global_periods.dat
cat fourier_inter* > global_interbin.dat
cat harmo* > global_harmonics.dat
cat candidate* > global_candidates.dat
