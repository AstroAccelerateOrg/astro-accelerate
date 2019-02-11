import ctypes
import numpy as np

lib = ctypes.CDLL('/home/carels/astro-accelerate/build/libastroaccelerate.so')

# Define ctypes for float pointers
FLOAT = ctypes.c_float
PFLOAT = ctypes.POINTER(FLOAT)
PPFLOAT = ctypes.POINTER(PFLOAT)
PPPFLOAT = ctypes.POINTER(PPFLOAT)

class filterbank_metadata_struct (ctypes.Structure):
    _fields_ = [
        ("m_tstart",ctypes.c_double),
        ("m_tsamp",ctypes.c_double),
        ("m_fch1", ctypes.c_double),
        ("m_foff",ctypes.c_double),
        ("m_nbits",ctypes.c_int),
        ("m_nsamples",ctypes.c_int),
        ("m_nchans",ctypes.c_int)
    ]

##
# \brief Python class that reads a sigproc file.
# \details Please see include/aa_sigproc_input.hpp for library implementation.
# \author Cees Carels.
# \date 05 February 2019.
#
class aa_py_sigproc_input:
    def __init__(self, path: str):
        lib.aa_py_sigproc_input.argtypes = [ctypes.c_char_p]
        lib.aa_py_sigproc_input.restype = ctypes.c_void_p
        c_string = ctypes.c_char_p(path.encode('utf-8'))
        self.m_obj = lib.aa_py_sigproc_input(c_string)
        
        print("Constructed aa_py_sigproc_input")
        # Call into library to construct object, and get the metadata and input data

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        print("Destructed aa_py_sigproc_input")

    def read_metadata(self):
        lib.aa_py_sigproc_input_read_metadata.argtypes = [ctypes.c_void_p]
        lib.aa_py_sigproc_input_read_metadata.restype = filterbank_metadata_struct
        return lib.aa_py_sigproc_input_read_metadata(self.m_obj)
        print("Read metadata")
        # Call into library to get filterbank_metadata.

    def read_signal(self):
        lib.aa_py_sigproc_input_read_signal.argtypes = [ctypes.c_void_p]
        lib.aa_py_sigproc_input_read_signal.restype = ctypes.POINTER(ctypes.c_ushort)
        return lib.aa_py_sigproc_input_read_signal(self.m_obj)
    def input_buffer():
        print("input_buffer")
        # Call into library to get input_buffer.

##
# \brief Python class for creating a filterbank_metadata object.
# \details Please see include/aa_filterbank_metadata.hpp for library implementation.
# \author Cees Carels.
# \date 05 February 2019.
#

class aa_py_filterbank_metadata:
    def __init__(self, tstart: float, tsamp: float, nbits: int, nsamples: int, fch1: float, foff: float, nchans: int):
        self.m_tstart = tstart
        self.m_tsamp = tsamp
        self.m_nbits = nbits
        self.m_nsamples = nsamples
        self.m_fch1 = fch1
        self.m_foff = foff
        self.m_nchans = nchans
        lib.aa_py_filterbank_metadata.argtypes = []
        lib.aa_py_filterbank_metadata.restype = ctypes.c_void_p
        self.m_obj = lib.aa_py_filterbank_metadata(ctypes.c_double(self.m_tstart), ctypes.c_double(self.m_tsamp), ctypes.c_int(self.m_nbits), ctypes.c_int(self.m_nsamples), ctypes.c_double(self.m_fch1), ctypes.c_double(self.m_foff), ctypes.c_int(self.m_nchans))
        
    def __exit__(self, exc_type, exc_value, traceback):
        lib.aa_py_filterbank_metadata_delete.argtypes = [ctypes.c_void_p]
        lib.aa_py_filterbank_metadata_delete(self.m_obj)
        print("Destructed aa_py_filterbank_metadata")

    def __enter__(self):
        return self

    def pointer(self):
        return self.m_obj
        
    def tstart(self):
        lib.aa_py_filterbank_metadata_tstart.argtypes = [ctypes.c_void_p]
        lib.aa_py_filterbank_metadata_tstart.restype = ctypes.c_double
        self.m_tstart = ctypes.c_double(lib.aa_py_filterbank_metadata_tstart(self.m_obj)).value
        return self.m_tstart
        
##
# \brief Python class for creating dm settings that can be added to an aa_py_ddtr_plan.
# \details Please see include/aa_ddtr_plan.hpp `struct dm` for library implementation.
# \author Cees Carels.
# \date 05 February 2019.
#
class aa_py_dm:
    def __init__(self, low: float, high: float, step: int, inBin: int, outBin: int):
        self.m_low    = low
        self.m_high   = high
        self.m_step   = step
        self.m_inBin  = inBin
        self.m_outBin = outBin

    def low(self):
        return self.m_low

    def high(self):
        return self.m_high

    def step(self):
        return self.m_step

    def inBin(self):
        return self.m_inBin

    def outBin(self):
        return self.m_outBin

    def __exit__(self, exc_type, exc_value, traceback):
        print("Destructed aa_py_dm")

##
# \brief Python class for creating a ddtr_plan.
# \details Please see include/aa_ddtr_plan for library implementation.
# \author Cees Carels.
# \date 05 February 2019.
#
class aa_py_ddtr_plan:
    def __init__(self, dm: np.array):
        lib.aa_py_ddtr_plan.argtypes = []
        lib.aa_py_ddtr_plan.restype = ctypes.c_void_p

        lib.aa_py_ddtr_plan_add_dm.argtypes = [ctypes.c_void_p, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_int, ctypes.c_int]
        lib.aa_py_ddtr_plan_add_dm.restype = ctypes.c_bool
        self.m_set_power = 0.0
        self.m_set_enable_msd_baseline_noise = False
        self.m_obj = lib.aa_py_ddtr_plan()
        
        if(dm.size):
            if(type(dm[0]) is not aa_py_dm):
                print("ERROR: Supplied dm is the wrong type, {}".format(type(dm[0]).__name__))
            else:
                self.m_dm = dm
                for dm in self.m_dm:
                    lib.aa_py_ddtr_plan_add_dm(self.m_obj, dm.low(), dm.high(), dm.step(), dm.inBin(), dm.outBin())
        else:
            print("ERROR: The array is empty.")

    def __exit__(self, exc_type, exc_value, traceback):
        lib.aa_py_ddtr_plan_delete.argtypes = [ctypes.c_void_p]
        lib.aa_py_ddtr_plan_delete(self.m_obj)
        print("Destructed aa_py_ddtr_plan")

    def __enter__(self):
        return self

    def pointer(self):
        return self.m_obj
            
    def set_power(self, power: float):
        lib.aa_py_ddtr_plan_set_power.argtypes = [ctypes.c_void_p, ctypes.c_float]
        lib.aa_py_ddtr_plan_set_power.restype = ctypes.c_bool
        self.m_set_power = power
        return lib.aa_py_ddtr_plan_set_power(self.m_obj, self.m_set_power)

    def power(self):
        lib.aa_py_ddtr_plan_power.argtypes = [ctypes.c_void_p]
        lib.aa_py_ddtr_plan_power.restype = ctypes.c_float
        self.m_power = ctypes.c_float(lib.aa_py_ddtr_plan_power(self.m_obj)).value
        return self.m_power

    def set_enable_msd_baseline_noise(self, enable_msd_baseline_noise: bool):
        lib.aa_py_ddtr_plan_set_enable_msd_baseline.argtypes = [ctypes.c_void_p, ctypes.c_bool]
        lib.aa_py_ddtr_plan_set_enable_msd_baseline.restype = ctypes.c_bool
        self.m_enable_msd_baseline_noise = enable_msd_baseline_noise
        return lib.aa_py_ddtr_plan_set_enable_msd_baseline(self.m_obj, self.m_set_enable_msd_baseline_noise)

    def enable_msd_baseline_noise(self):
        lib.aa_py_ddtr_plan_enable_msd_baseline_noise.argtypes = [ctypes.c_void_p]
        lib.aa_py_ddtr_plan_enable_msd_baseline_noise.restype = ctypes.c_bool
        self.m_enable_msd_baseline_noise = ctypes.c_float(lib.aa_py_ddtr_plan_enable_msd_baseline_noise).value
        return self.m_enable_msd_baseline_noise

    def print_info(self):
        print("AA_PY_DDTR_PLAN INFORMATION:")
        if(self.m_dm.size):
            for i in range(self.m_dm.size):
                print("     aa_py_ddtr_plan range {}: low {}, high {}, step {}, inBin {}, outBin {}".format(i, self.m_dm[i].m_low, self.m_dm[i].m_high, self.m_dm[i].m_step, self.m_dm[i].m_inBin, self.m_dm[i].m_outBin))
        else:
            print("No dm ranges have been provided.")
        print("     aa_py_ddtr_plan power {}".format(self.m_set_power))
        print("     aa_py_ddtr_plan enable_msd_baseline_noise: {}".format(self.m_enable_msd_baseline_noise))

##
# \brief Class for configuring a ddtr_strategy.
# \details Please see include/aa_ddtr_strategy.hpp for library implementation.
# \details This class is used to provide a Python object from the library when the user binds an aa_py_ddtr_plan.
# \author Cees Carels.
# \date 05 February 2019.
#
class aa_py_ddtr_strategy():
    def __init__(self, aa_py_ddtr_plan):
        print("Constructed aa_py_ddtr_strategy")
    def __exit__(self, exc_type, exc_value, traceback):
        print("Destructed aa_py_ddtr_strategy")
        
##
# \brief Class for configuring an analysis_plan.
# \details Please see include/aa_analysis_plan.hpp for library implementation.
# \author Cees Carels.
# \date 05 February 2019.
#
class aa_py_analysis_plan():
    def __init__(self, sigma_cutoff: float, sigma_constant: float, max_boxcar_width_in_sec: float, candidate_algorithm: bool, enable_msd_baseline_noise: bool):
        self.m_sigma_cutoff = sigma_cutoff
        self.m_sigma_constant = sigma_constant
        self.m_max_boxcar_width_in_sec = max_boxcar_width_in_sec
        self.m_candidate_algorithm = candidate_algorithm
        self.m_enable_msd_baseline_noise = enable_msd_baseline_noise
        # Call into library using the existing aa_ddtr_strategy.

    def __exit__(self, exc_type, exc_value, traceback):
        print("Destructed aa_py_analysis_plan")

    def __enter__(self):
        return self

    def sigma_cutoff(self):
        return self.m_sigma_cutoff

    def sigma_constant(self):
        return self.m_sigma_constant

    def max_boxcar_width_in_sec(self):
        return self.m_max_boxcar_width_in_sec

    def candidate_algorithm(self):
        return self.m_candidate_algorithm

    def enable_msd_baseline_noise(self):
        return self.m_enable_msd_baseline_noise

    def print_info(self):
        print("AA_PY_ANALYSIS_PLAN INFORMATION:")
        print("     aa_py_analysis_plan sigma_cutoff {}".format(self.m_sigma_cutoff))
        print("     aa_py_analysis_plan sigma_constant {}".format(self.m_sigma_constant))
        print("     aa_py_analysis_plan max_boxcar_width_in_sec {}".format(self.m_max_boxcar_width_in_sec))
        print("     aa_py_analysis_plan candidate_algorithm {}".format(self.m_candidate_algorithm))
        print("     aa_py_analysis_plan enable_msd_baseline_noise {}".format(self.m_enable_msd_baseline_noise))
              

##
# \brief Class for configuring an analysis_strategy
# \details Please see include/aa_analysis_strategy.hpp for library implementation.
# \author Cees Carels.
# \date 05 February 2019.
#
class aa_py_analysis_strategy():
    def __init__(self):
        print("Constructed aa_py_analysis_strategy")

    def __exit__(self, exc_type, exc_value, traceback):
        print("Destructed aa_py_analysis_strategy")

##
# \brief Class for configuring a periodicity_plan.
# \details Please see include/aa_periodicity_plan.hpp for library implementation.
# \author Cees Carels.
# \date 05 February 2019.
#
class aa_py_periodicity_plan():
    def __init__(self, sigma_cutoff: float, sigma_constant: float, nHarmonics: int, export_powers: int, candidate_algorithm: bool, enable_msd_baseline_noise: bool):
        self.m_sigma_cutoff = sigma_cutoff
        self.m_sigma_constant = sigma_constant
        self.m_nHarmonics = nHarmonics
        self.m_export_powers = export_powers
        self.m_candidate_algorithm = candidate_algorithm
        self.m_enable_msd_baseline_noise = enable_baseline_noise

    def __exit__(self, exc_type, exc_value, traceback):
        print("Destructed aa_py_periodicity_plan")

    def sigma_cutoff(self):
        return self.m_sigma_cutoff

    def sigma_constant(self):
        return self.m_sigma_constant

    def nHarmonics(self):
        return self.m_nHarmonics

    def export_powers(self):
        return self.m_export_powers

    def candidate_algorithm(self):
        return self.m_candidate_algorithm

    def enable_msd_baseline_noise(self):
        return self.m_enable_msd_baseline_noise

##
# \brief Class for configuring a periodicity_strategy
# \details Please see include/aa_periodicity_strategy.hpp for library implementation.
# \author Cees Carels.
# \date 05 February 2019.
#
class aa_py_periodicity_strategy():
    def __init__(self):
        print("Constructed aa_py_periodicity_strategy")

    def __exit__(self, exc_type, exc_value, traceback):
        print("Destructed aa_py_periodicity_strategy")

        
##
# \brief Class for configuring an fdas_plan.
# \details Please see include/aa_fdas_plan.hpp for library implementation.
# \author Cees Carels.
# \date 05 February 2019.
#
class aa_py_fdas_plan():
    def __init__(self, sigma_cutoff, float, sigma_constant: float, num_boots: int, num_trial_bins: int, navdms: int, narrow: float, wide: float, nsearch: int, aggression: float, enable_msd_baseline_noise: bool):
        self.m_sigma_cutoff = sigma_cutoff
        self.m_sigma_constant = sigma_constant
        self.m_num_boots = num_boots
        self.m_num_trial_bins = num_trial_bins
        self.m_navdms = navdms
        self.m_narrow = narrow
        self.m_wide = wide
        self.m_nsearch = nsearch
        self.m_aggression = aggression
        self.m_enable_msd_baseline_noise = enable_msd_baseline_noise

    def __exit__(self, exc_type, exc_value, traceback):
        print("Destructed aa_pu_fdas_plan")

##
# \brief Class for configuring an fdas_strategy
# \details Please see include/aa_fdas_strategy.hpp for library implementation.
# \author Cees Carels.
# \date 05 February 2019.
#
class aa_py_fdas_strategy():
    def __init__(self):
        print("Constructed aa_py_fdas_strategy")
        # Call into library to obtain aa_fdas_strategy.

    def __exit__(self, exc_type, exc_value, traceback):
        print("Destructed aa_py_fdas_strategy")

##
# \brief Class for interacting with aa_pipeline_api objects from the library.
# \details Please see include/aa_pipeline_api.hpp for library implementation.
# \author Cees Carels.
# \date 05 February 2019.
#
class aa_py_pipeline():
    def __init__(self, metadata: filterbank_metadata_struct, input_data: ctypes.POINTER(ctypes.c_ushort), card_number: int):
        lib.aa_py_pipeline_api.argtypes = [filterbank_metadata_struct, ctypes.POINTER(ctypes.c_ushort), ctypes.c_int]
        lib.aa_py_pipeline_api.restype = ctypes.c_void_p
        self.m_obj = lib.aa_py_pipeline_api(metadata, input_data, ctypes.c_int(card_number))

    def __exit__(self, exc_type, exc_value, traceback):
        print("Destructed aa_py_pipeline")

    def __enter__(self):
        return self
        
    def bind_ddtr_plan(self, plan):
        lib.aa_py_pipeline_api_bind_ddtr_plan.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
        lib.aa_py_pipeline_api_bind_ddtr_plan.restype = ctypes.c_bool
        lib.aa_py_pipeline_api_bind_ddtr_plan(self.m_obj, plan.pointer())
        # Call into library to bind plan

    def bind_analysis_plan(self, plan):
        lib.aa_py_pipeline_api_bind_analysis_plan.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
        lib.aa_py_pipeline_api_bind_analysis_plan.restype = ctypes.c_bool
        return lib.aa_py_pipeline_api_bind_analysis_plan(self.m_obj, plan.pointer())
        
    def ddtr_strategy(self):
        lib.aa_py_pipeline_api_ddtr_strategy.argtypes = [ctypes.c_void_p]
        lib.aa_py_pipeline_api_ddtr_strategy.restype = ctypes.c_void_p
        self.m_ddtr_strategy = lib.aa_py_pipeline_api_ddtr_strategy(self.m_obj)
        print("ddtr_strategy")
        # Call into library to retrieve aa_ddtr_strategy and then set to aa_py_ddtr_strategy.

    def analysis_strategy(self):
        print("analysis_strategy")
        # Call into library to retrieve aa_analysis_strategy and then set to aa_py_analysis_strategy.

    def periodicity_strategy(self):
        print("periodicity_strategy")
        # Call into library to retrieve aa_periodicity_strategy and then set to aa_py_periodicity_strategy.

    def fdas_strategy(self):
        print("fdas_strategy")
        # Call into library to retrieve aa_fdas_strategy and then set to aa_py_fdas_strategy.

    ## \brief Runs the entire pipeline end-to-end, saving to disk. #
    def run_save_to_disk(self):
        print("run_save_to_disk")
        # Call into library to run pipeline.

    ## \brief Runs the entire pipeline end-to-end. #
    def run(self):
        print("run")
        # Call into library to run pipeline.

    def analysis_output_store(self):
        print("analysis_output_store")
        # Call into library to obtain all analysis output from the output_store.

    def periodicity_output_store(self):
        print("periodicity_output_store")
        # Call into library to obtain all periodicity output from the output_store.

    def fdas_output_store(self):
        print("fdas_output_store")
        # Call into library to obtain all fdas output from the output_store.

    ## \brief Returns a pointer to the dedispersed output_buffer in the library. #
    def get_buffer(self):
        #lib.my_class_get_buffer.argtypes = [ctypes.c_void_p]
        #lib.my_class_get_buffer.restype = PPPFLOAT
        print("Getting buffer")
        print(self.obj)
        return lib.my_class_get_buffer(self.obj)

    def get_filterbank_metadata():
        print("get_filterbank_metadta")
        # Call into library to obtain the filterbank_metadata.
