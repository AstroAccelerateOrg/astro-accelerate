/*************************************************************************
    This is GPU implementation of a polyphase filter. 
    Copyright (C) 2015  Adamek Karel, Novotny Jan, Armour Wes

    This file is part of Astro-Accelerate PolyPhase Filter (AAPPF).

    AAPPF is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AAPPF is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AAPPF.  If not, see <http://www.gnu.org/licenses/>.

**************************************************************************/

#include <fftw3.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

void Normalize_FFT(float2 *output, int nChannels, int nSpectra, double factor){
	int s,c;
	for(s=0;s<nSpectra;s++){
		for(c=0;c<nSpectra;c++){
			output[s*nChannels+c].x=output[s*nChannels+c].x/factor;
			output[s*nChannels+c].y=output[s*nChannels+c].y/factor;
		}
	}
}


void Perform_FFT(float2 *spectra, int nChannels){
	int N=nChannels;
	fftwf_plan p;
	
	p = fftwf_plan_dft_1d(N, (fftwf_complex *) spectra, (fftwf_complex *) spectra, FFTW_FORWARD, FFTW_ESTIMATE);
	fftwf_execute(p);
	fftwf_destroy_plan(p);
}


void FIR_check(float2 *input_data, float2 *spectra_GPU, float *coeff, int nTaps, int nChannels, int nSpectra, double *cumulative_error, double *mean_error){
	float2 ftemp;
	double etemp=0;
	double error=0;
	for(int bl=0;bl<nSpectra;bl++){
		for(int c=0;c<nChannels;c++){
			ftemp.x=0;ftemp.y=0;
			for(int t=0;t<nTaps;t++){
				ftemp.x+=coeff[t*nChannels + c]*input_data[bl*nChannels + t*nChannels + c].x;
				ftemp.y+=coeff[t*nChannels + c]*input_data[bl*nChannels + t*nChannels + c].y;
			}//nTaps
			etemp=abs(ftemp.x-spectra_GPU[bl*nChannels + c].x);
			error+=etemp;
			etemp=abs(ftemp.y-spectra_GPU[bl*nChannels + c].y);
			error+=etemp;
		}//nChannels
	} //nBlocks
	*cumulative_error=error;
	*mean_error=error/(double) (nChannels*nSpectra);
}

void FIR_FFT_check(float2 *input_data, float2 *spectra_GPU, float *coeff, int nTaps, int nChannels, int nSpectra, double *cumulative_error, double *mean_error){
	float2 *spectra;
	spectra	= (float2 *)malloc(nChannels*sizeof(float2));
	double etemp=0;
	double error=0;
	for(int bl=0;bl<nSpectra;bl++){
		for(int c=0;c<nChannels;c++){
			spectra[c].x=0.0;spectra[c].y=0.0;
			for(int t=0;t<nTaps;t++){
				spectra[c].x+=coeff[t*nChannels + c]*input_data[bl*nChannels + t*nChannels + c].x;
				spectra[c].y+=coeff[t*nChannels + c]*input_data[bl*nChannels + t*nChannels + c].y;
			}//nTaps
		}//nChannels
		Perform_FFT(spectra, nChannels);
		
		for(int c=0;c<nChannels;c++){	
			etemp=abs(spectra[c].x-spectra_GPU[bl*nChannels + c].x);
			error+=etemp;
			etemp=abs(spectra[c].y-spectra_GPU[bl*nChannels + c].y);
			error+=etemp;
		}
		
	} //nSpectra
	*cumulative_error=error;
	*mean_error=error/(nChannels*nSpectra);
	delete[] spectra;
}
