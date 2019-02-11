from python_wrapper import *

with aa_py_sigproc_input("/mnt/data/AstroAccelerate/filterbank/BenMeerKAT.fil") as sigproc_input:
    data = sigproc_input.read_metadata()


tstart = 1.0
tsamp = 1.0
nbits = 8
nsamples = 1024
fch1 = 1.0
foff = 1.0
nchans = 1024
with aa_py_filterbank_metadata(tstart, tsamp, nbits, nsamples, fch1, foff, nchans) as filterbank_metadata:
    print("aa_py_filterbank_metadata tstart {}".format(filterbank_metadata.tstart()))
    filterbank_metadata.pointer()

# ddtr_plan settings
dm1 = aa_py_dm(1, 2, 1, 1, 1)
dm2 = aa_py_dm(2, 3, 1, 1, 1)
dm_list = np.array([dm1, dm2], dtype=aa_py_dm)
power = 5.0
enable_msd_baseline_noise=False

# Create ddtr_plan
ddtr_plan = aa_py_ddtr_plan(dm_list)
ddtr_plan.set_power(power)
ddtr_plan.set_enable_msd_baseline_noise(enable_msd_baseline_noise)
ddtr_plan.print_info()

sigma_cutoff = 1.0
sigma_constant = 1.0
max_boxcar_width_in_sec = 1.0
candidate_algorithm = False
enable_msd_baseline_noise = False

with aa_py_analysis_plan(sigma_cutoff, sigma_constant, max_boxcar_width_in_sec, candidate_algorithm, enable_msd_baseline_noise) as analysis_plan:
    analysis_plan.print_info()


# Create an aa_py_pipeline
with aa_py_pipeline(data, sigproc_input.read_signal(), 0) as pipeline:
    pipeline.bind_ddtr_plan(ddtr_plan) # Bind the ddtr_plan
    
    # Run the pipeline
    #pipeline.run()
    
    # Obtain the dedispersed output data buffer
    #ptr = pipeline.get_buffer()
    #print(ptr[1][1][2])


del ddtr_plan
