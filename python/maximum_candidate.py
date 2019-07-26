#!/usr/bin/env python3
##
# @package use_python use_python.py
#

from __future__ import print_function

# Check Python version on this machine
import sys
import struct
import time
if (sys.version_info < (3, 0)):
    print("ERROR: Python version less than 3.0. Exiting...")
    sys.exit()
from py_astro_accelerate import *


# Open filterbank file for reading metadata and signal data
filedata = str(sys.argv[1])
print(filedata)
#sigproc_input = aa_py_sigproc_input("/home/novotny/filterbank/injectfrb_nfrb500_DM20-1992_5183sec_20190610-1710.fil")
sigproc_input = aa_py_sigproc_input(filedata)
metadata = sigproc_input.read_metadata()
if not sigproc_input.read_signal():
    print("ERROR: Invalid .fil file path. Exiting...")
    sys.exit()
input_buffer = sigproc_input.input_buffer()

# ddtr_plan settings
# settings: aa_py_dm(low, high, step, inBin, outBin)
dm1 = aa_py_dm(0, 150, 0.1, 1, 1)
dm2 = aa_py_dm(150, 300, 0.2, 1, 1)
dm3 = aa_py_dm(300, 500, 0.25, 1, 1)
dm4 = aa_py_dm(500, 900, 0.4, 2, 2)
dm5 = aa_py_dm(900, 1200, 0.6, 4, 4)
dm6 = aa_py_dm(1200, 1500, 0.8, 4, 4)
dm7 = aa_py_dm(1500, 2000, 1.0, 4, 4)
dm8 = aa_py_dm(2000, 3000, 2.0, 8, 8)
dm_list = np.array([dm1, dm2, dm3, dm4, dm5, dm6, dm7, dm8],dtype=aa_py_dm)

# Create ddtr_plan
ddtr_plan = aa_py_ddtr_plan(dm_list)
enable_msd_baseline_noise=True
ddtr_plan.set_enable_msd_baseline_noise(enable_msd_baseline_noise)
ddtr_plan.print_info()

# Set up pipeline components
pipeline_components = aa_py_pipeline_components()
pipeline_components.dedispersion = True
pipeline_components.analysis = True

# Create analysis plan
sigma_cutoff = 6
sigma_constant = 4.0
max_boxcar_width_in_sec = 0.5
candidate_algorithm = False
enable_msd_baseline_noise = True

analysis_plan = aa_py_analysis_plan(sigma_cutoff, sigma_constant, max_boxcar_width_in_sec, candidate_algorithm, enable_msd_baseline_noise)
analysis_plan.print_info()

# Set up pipeline component options
pipeline_options = aa_py_pipeline_component_options()
pipeline_options.output_dmt = False

# Select GPU card number on this machine
card_number = 0

# Create pipeline
pipeline = aa_py_pipeline(pipeline_components, pipeline_options, metadata, input_buffer, card_number)
pipeline.bind_ddtr_plan(ddtr_plan) # Bind the ddtr_plan
pipeline.bind_analysis_plan(analysis_plan) # Bind the analysis plan

# Run the pipeline with AstroAccelerate
DM = []
SNR = []
WIDTH = []
TS = []
TIME = []

while (pipeline.run()):
    print("NOTICE: Python script running over next chunk")
    # during the run get the Candidates
    if pipeline.status_code() == 1:
        start = time.time()
        (nCandidates, dm, ts, snr, width, c_range, c_tchunk, ts_inc)=pipeline.get_candidates()
        end = time.time()
        print(bcolors.WARNING + "Time to read: " + str(end - start) + bcolors.ENDC)
        print("Number of candidates", nCandidates)
        start = time.time()
        if (nCandidates > 0):
            (tmp_dm, tmp_snr, tmp_time_sample, tmp_time, tmp_width) = SPD.maximum(metadata, pipeline, ddtr_plan, ts_inc, nCandidates, dm, ts, snr, width, c_range, c_tchunk)
        DM.append(tmp_dm)
        SNR.append(tmp_snr)
        TS.append(tmp_time_sample)
        TIME.append(tmp_time)
        WIDTH.append(tmp_width)
        end = time.time()
        # Write the candidates to disk
        print(bcolors.WARNING + "Time to find maximum: " + str(end - start) + bcolors.ENDC)
    if pipeline.status_code() == -1:
        print("ERROR: Pipeline status code is {}. The pipeline encountered an error and cannot continue.".format(pipeline.status_code()))
        break

#write the Maximum candidate to the disk and print to the console
SPD.write_maximum(DM, SNR, TIME, TS, WIDTH)


# Cleaning the plans
del ddtr_plan
del analysis_plan
