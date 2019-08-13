#!/usr/bin/env python3
##
# @package use_python use_python.py
#

from __future__ import print_function

# Check Python version on this machine
import sys
import struct
if (sys.version_info < (3, 0)):
    print("ERROR: Python version less than 3.0. Exiting...")
    sys.exit()
from py_astro_accelerate import *


dDMs      = [0.2, 0.3, 0.5, 1.0]
dsubDMs   = [4.8, 7.2, 12.0, 24.0]
downsamps = [1, 2, 4, 8]
subcalls  = [8, 3, 4, 4]
startDMs  = [0.0, 38.4, 60.0, 108.0]
dmspercall = 24
nsub = 32
numout = 72000
basename = "aa_dm90"
# cannot read more than one filterbank now, we can not use wildcards
# how multiple files are then interpreted?
rawfiles = basename+".fil"

# Open filterbank file for reading metadata and signal data
sigproc_input = aa_py_sigproc_input(rawfiles)
metadata = sigproc_input.read_metadata()
if not sigproc_input.read_signal():
    print("ERROR: Invalid .fil file path. Exiting...")
    sys.exit()
input_buffer = sigproc_input.input_buffer()

# ddtr_plan settings
# settings: aa_py_dm(low, high, step, inBin, outBin)
temp_list = []
for i in range(len(dDMs)):
    highDM = startDMs[i] + dmspercall*subcalls[i]*dDMs[i]
    tmp = aa_py_dm(startDMs[i],highDM,dDMs[i],downsamps[i],1)
    temp_list.append(tmp)
dm_list = np.array(temp_list,dtype=aa_py_dm)

# Create ddtr_plan
ddtr_plan = aa_py_ddtr_plan(dm_list)
enable_msd_baseline_noise=False
ddtr_plan.set_enable_msd_baseline_noise(enable_msd_baseline_noise)
ddtr_plan.print_info()

# Set up pipeline components
pipeline_components = aa_py_pipeline_components()
pipeline_components.dedispersion = True

# Set up pipeline component options
pipeline_options = aa_py_pipeline_component_options()
pipeline_options.output_dmt = True
#Need to be enabled otherwise there will be no data copy from GPU memory to host memory
pipeline_options.copy_ddtr_data_to_host = True;

# Select GPU card number on this machine
card_number = 0

# Create pipeline
pipeline = aa_py_pipeline(pipeline_components, pipeline_options, metadata, input_buffer, card_number)
pipeline.bind_ddtr_plan(ddtr_plan) # Bind the ddtr_plan

# Run the pipeline with AstroAccelerate
while (pipeline.run()):
    print("NOTICE: Python script running over next chunk")
    if pipeline.status_code() == -1:
        print("ERROR: Pipeline status code is {}. The pipeline encountered an error and cannot continue.".format(pipeline.status_code()))
        break

# Get the data of DDTR to python
(ts_inc, ddtr_output) = pipeline.get_buffer()

for pos_range in range(pipeline.ddtr_range()):
    list_ndms = pipeline.ddtr_ndms()
    for n_dms in range(list_ndms[pos_range]):
        DM = pipeline.dm_low(pos_range) + dDMs[pos_range]*n_dms
        filename = basename + "_DM" + "{:07.2f}".format(DM)
        result_file = filename + ".dat"
        print("Writing results to: " + result_file, end='\r')
        newfile = open(result_file, "wb")
        nsamp_for_range = int(ts_inc/downsamps[pos_range])
        header.information_file(filename,nsamp_for_range, DM, downsamps[pos_range], metadata)
        for samples_pos in range(nsamp_for_range):
            newfile.write(struct.pack('f',ddtr_output[pos_range][n_dms][samples_pos]*downsamps[pos_range]))

if pipeline.status_code() == -1:
    print("ERROR: Pipeline status code is {}. The pipeline encountered an error and cannot continue.".format(pipeline.status_code()))

# Cleaning the plans
pipeline.cleanUp()
del ddtr_plan

#Presto part
#for dDM, dsubDM, downsamp, subcall, startDM in \
#        zip(dDMs, dsubDMs, downsamps, subcalls, startDMs):
#    print("StartDM: ", startDM, dDM)
#    for ii in range(subcall):
#        subDM = startDM + (ii + 0.5)*dsubDM
#        loDM = startDM + ii*dsubDM
#        print(ii,loDM)
#

#print(pipeline.dm_low(1))
