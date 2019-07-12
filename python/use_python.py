#!/usr/bin/env python3
##
# @package use_python use_python.py
#

from __future__ import print_function

# Check Python version on this machine
import sys
if (sys.version_info < (3, 0)):
    print("ERROR: Python version less than 3.0. Exiting...")
    sys.exit()
from py_astro_accelerate import *

# Open filterbank file for reading metadata and signal data (please provide your filterbank file0)
sigproc_input = aa_py_sigproc_input("<input_some_filterbank_data.fil>")
metadata = sigproc_input.read_metadata()
if not sigproc_input.read_signal():
    print("ERROR: Invalid .fil file path. Exiting...")
    sys.exit()
input_buffer = sigproc_input.input_buffer()

# ddtr_plan settings
# settings: low, high, step, inBin, outBin.
dm1 = aa_py_dm(0, 370, 0.307, 1, 1)
dm2 = aa_py_dm(370, 740, 0.652, 2, 2)
dm3 = aa_py_dm(740, 1480, 1.266, 4, 4)
dm_list = np.array([dm1, dm2, dm3], dtype=aa_py_dm)
power = 2.0
enable_msd_baseline_noise=False

# Create ddtr_plan
ddtr_plan = aa_py_ddtr_plan(dm_list)
ddtr_plan.set_enable_msd_baseline_noise(enable_msd_baseline_noise)
ddtr_plan.print_info()

# Create analysis plan
sigma_cutoff = 6
sigma_constant = 4.0
max_boxcar_width_in_sec = 0.05
candidate_algorithm = True
enable_msd_baseline_noise = False

analysis_plan = aa_py_analysis_plan(sigma_cutoff, sigma_constant, max_boxcar_width_in_sec, candidate_algorithm, enable_msd_baseline_noise)
analysis_plan.print_info()

# Create periodicity plan
nHarmonics = 32
export_powers = 2
periodicity_plan = aa_py_periodicity_plan(sigma_cutoff, sigma_constant, nHarmonics, export_powers, candidate_algorithm, enable_msd_baseline_noise)

# Create fdas plan
fdas_plan = aa_py_fdas_plan(sigma_cutoff, sigma_constant, enable_msd_baseline_noise)

# Set up pipeline components
pipeline_components = aa_py_pipeline_components()
pipeline_components.dedispersion = True
pipeline_components.analysis = True
pipeline_components.periodicity = False
pipeline_components.fdas = False

# Set up pipeline component options
pipeline_options = aa_py_pipeline_component_options()
pipeline_options.zero_dm = False
pipeline_options.zero_dm_with_outliers = False
pipeline_options.old_rfi = False
pipeline_options.msd_baseline_noise = enable_msd_baseline_noise
pipeline_options.output_dmt = True
pipeline_options.output_ffdot_plan = False
pipeline_options.output_fdas_list = False
pipeline_options.candidate_algorithm = False
pipeline_options.fdas_custom_fft = True
pipeline_options.fdas_inbin = False
pipeline_options.fdas_norm = False

# Select GPU card number on this machine
card_number = 0

# Create pipeline
pipeline = aa_py_pipeline(pipeline_components, pipeline_options, metadata, input_buffer, card_number)
pipeline.bind_ddtr_plan(ddtr_plan) # Bind the ddtr_plan
pipeline.bind_analysis_plan(analysis_plan)
pipeline.bind_periodicity_plan(periodicity_plan)
pipeline.bind_fdas_plan(fdas_plan)

while (pipeline.run()):
    print("NOTICE: Python script running over next chunk")
    if pipeline.status_code() == 1:      
        #this is the example hot to get SPD candidates to the python
        (nCandidates, dm, ts, snr, width, c_range, c_tchunk, ts_inc)=pipeline.get_candidates()
    if pipeline.status_code() == -1:
        print("ERROR: Pipeline status code is {}. The pipeline encountered an error and cannot continue.".format(pipeline.status_code()))
        break
        
if pipeline.status_code() == -1:
    print("ERROR: Pipeline status code is {}. The pipeline encountered an error and cannot continue.".format(pipeline.status_code()))

del fdas_plan
del periodicity_plan
del analysis_plan
del ddtr_plan
