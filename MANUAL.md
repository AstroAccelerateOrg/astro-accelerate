# User Manual
## Hardware and software support status

AstroAccelerate supported GPU status

| Compute Architecture | Architecture Name | Supported (Y/N) |
|----------------------|-------------------|-----------------|
| 2.0                  | Fermi             | N               |
| 2.1                  | Fermi             | N               |
| 3.0                  | Fermi             | N               |
| 3.2                  | Kepler            | N               |
| 3.5                  | Kepler            | Y               |
| 3.7                  | Kepler            | Y               |
| 5.0                  | Maxwell           | Y               |
| 5.2                  | Maxwell           | Y               |
| 5.3                  | Maxwell           | Y               |
| 6.0                  | Pascal            | Y               |
| 6.1                  | Pascal            | Y               |
| 6.2                  | Pascal            | -               |
| 7.0                  | Volta             | Y               |

AstroAccelerate supported CUDA&reg; SDK status

| SDK Version | Supported (Y/N) |
|-------------|-----------------|
| 8.0         | Y               |
| 9.0         | Y               |
| 10.0        | N               |

AstroAccelerate supported build system versions
* AstroAccelerate can be compiled via `make` or via `CMake` (> 2.6.0) which then generates a Makefile that can be used to create the shared object library and an executable which links against the library.

AstroAccelerate supported compilers
* The compiler for all device code is [nvcc](https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html).
* For C++, the minimum recommended compiler is indicated on a per version basis.
For more information, please see the [NVIDIA CUDA Installation Guide](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html).

AstroAccelerate supported operating systems
* In general, AstroAccelerate can be compiled on any system where the build system and compiler requirements are met.
Only Linux distributions are officially supported at this time.
Please see the [NVIDIA CUDA Installation Guide](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html) for information about the operating system support for CUDA.

## How to obtain, compile, and run AstroAccelerate
The code can be obtained from the [astro-accelerate GitHub repository](https://github.com/AstroAccelerateOrg/astro-accelerate).

```
git clone git@github.com:AstroAccelerateOrg/astro-accelerate.git
```

Or
```
git clone https://github.com/AstroAccelerateOrg/astro-accelerate.git
```

Then go to into the code directory
```
cd astro-accelerate/
```

Please see the [`README.md`](README.md) in the top-level directory of the repository for further instructions on how to compile and run AstroAccelerate.
## How to configure the AstroAccelerate standalone using the input_file
Please see the [`README.md`](README.md) in the top-level directory of the repository
for instructions on how to obtain, compile, and run AstroAccelerate.

## How to use the AstroAccelerate API in your own code
In addition to being able to run as a standalone application, AstroAccelerate also offers an API, enabling it to be used in your own application code by linking against the AstroAccelerate shared object library object (`.so`) that is generated at compile-time.

### Selecting components
The AstroAccelerate API offers the user the option to run a number of combinations of analysis components. The user is free to choose which components are to be run, but the API will validate that the combination is valid, and the API will determine the order in which those components are executed. The API executes the components as a `pipeline`. The user will therefore mostly interact with the `aa_pipeline_api` class, and may find it useful to refer to the Doxygen documentation for this class, in addition to the guidance provided here.

A `pipeline` is nothing more than a collection of `components`. There are also `pipeline_option`s which consist of `component_option`s which further specify how the components run. The components and component options are documented in `aa_pipeline.hpp`.

When a user first creates an `aa_pipeline_api` object, they must provide a `pipeline` which contain the components the user wishes to run. They must also provide the `pipeline_option` which contains the component options. In addition, the user must provide an `aa_filterbank_metadata` object, which contains the information about the time series data from a telescope (such as sampling rate, number of samples, and frequency binning information). Furthermore, a pointer to the raw input data must be provided (for a telescope `.fil` this can be obtained alongside the metadata by using the `aa_sigproc_input` class, but the user can also provide their own array data). Lastly, the user must configure a GPU card to run the pipeline on, using the `aa_device_info` class (which they can use to query the GPUs on the machine).

### Configuring the plan and strategy
Once the `aa_pipeline_api` object is constructed, the user must still `bind` appropriate `plan` objects for each of the components they wish to run. There is one `plan` object for each `component`. For example, in order to run the `analysis` component, the user must `bind` an `aa_analysis_plan`. The `plan` contains the user's desired settings for how to run the component (such as the binning interval, and the dedispersion measure ranges to search for). In reality, user's settings may be sub-optimal from a performance perspective, so AstroAccelerate will optimise the user's `plan` by creating what is called a `strategy`. The `strategy` is created from the user's `plan`, and is the best compromise between the desired settings and good performance on the GPU. Strategy objects are required in order to run a pipeline, and the user is advised to verify and review the `strategy` when performing data analysis.

When `bind`ing the `plan`, the API `bind` call will return a boolean to indicate whether the `bind` was successful. The `strategy` can be obtained by calling a method by the component name and `_strategy()` added to it. For example, the `periodicity_strategy()` methods returns the currently calculated `aa_periodicity_strategy` based on the currently bound `aa_periodicity_plan` that the user provided.

All `strategy` objects have a `ready` state associated to them. All `strategy` objects must have a ready state set to `true` or else the API will not run the pipeline. The ready state of any strategy can be queried using the object's `ready()` method. Furthermore, the `aa_pipeline_api`'s `ready()` method will check whether all plan objects were supplied and whether all strategy objects are ready.

The `aa_pipeline_api` offers a `run()` method, which will run the pipeline end-to-end.

There is a complete end-to-end example of using the API in `include/aa_pipeline_generic.hpp`.

There are also C-style wrapper functions (located in `include/aa_pipeline_wrapper.hpp`) for the `include/aa_pipeline_generic.hpp` code, which simplify the API interaction even further. An example is located in `examples/src/filterbank_dedispersion.cpp`.

Alternatively, a more advanced library user may wish to bypass the API, and interact directly with the pipeline classes that execute the components. This comes with more control (such as obtaining the pipeline output iteratively for each processed time chunk). However, this approach contains fewer safety checks than when using the API. The main differences are that the strategy objects are passed directly to the pipeline (bypassing the API), and the user must run the pipeline over all time chunks themselves. A good example code is located in `examples/src/periodicity.cpp`.

### Checking for errors
Boolean return values for functions return `true` for `success` or `ready`, whereas `false` indicates `failure` or `not ready`. It is advisable to check for errors at each step of configuring `plan` and `strategy` objects, and the API will help with this as well.

## Component inputs and outputs
The inputs and outputs of the each of the components are tabulated below.

| Component       | Description                                                      | Configuration objects                         | **plan** input / **strategy** input                                                                                                                                                          | Component output                                                                                                                   |
|--------------|------------------------------------------------------------------|-----------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------|
| dedispersion | Performs dedispersion of the input data.                         | aa_ddtr_plan / aa_ddtr_strategy               | **plan**: dedispersion measure low, high, step, inBin, outBin, enable_msd_baseline_noise. <br> **strategy**: aa_ddtr_plan, aa_filterbank_metadata, amount of gpu memory, flag to indicate if analysis will be used or not. | dedispersed time chunk data (`std::vector<unsigned short>`)                                                                       |
| analysis     | Performs single pulse search (SPS) analysis of dedispersed data. | aa_analysis_plan / aa_analysis_strategy       | **plan**: aa_ddtr_strategy, sigma_cutoff, sigma_constant, max_boxcar_width_in_sec, candidate_algorithm (flag), enable_msd_baseline_noise (flag). <br> **strategy**: aa_analysis_plan.           | `analysis_output` struct (containing the data, dm_low, and dm_high for the time chunk) - see [include/aa_device_analysis.hpp](include/aa_device_analysis.hpp) for further details. |
| periodicity  | Performs a periodic peak finding analysis of dedispersed data.   | aa_periodicity_plan / aa_periodicity_strategy | **plan**: sigma_cutoff, sigma_constant, nHarmonics, export_powers, candidate_algorithm (flag), enable_msd_baseline_noise (flag). <br> **strategy**: aa_periodicity_plan.                     | Output saved to disk in bespoke format.                                                                                         |
| fdas         | Fourier Domain Accelerated Search of the dedispersed data.       | aa_fdas_plan / aa_fdas_strategy               | **plan**: sigma_cutoff, sigma_constant, enable_msd_baseline_noise. <br> **strategy**: aa_fdas_plan.                                                                   | Output saved to disk in bespoke format.                                                                                         |

Where more than one component takes an input by the same name, the standalone code will set the numerical values to be equal for all modules. The library user may, if they wish, configure this flag separately for all modules.

### Obtaining the output
The default for output is to save to disk. If interacting directly with the pipeline, the user may call any one of the overloaded `run` methods on the pipeline (using flags to indicate whether the output will be `dump`ed to disk and/or to the user at runtime), providing a `std::vector` to obtain dedispersed output from the `dedispersion` component, and an `aa_analysis_output` struct for the analysis single pulse search (SPS) output from the `analysis` component.

## Using the built-in logging facility
Users and developers alike may wish to use the logging facility provided by the `aa_log` class in order to log their activities to the `console` or to save it to a `file` on disk.

AstroAccelerate uses the logging class internally to provide output to the screen by default.

The logging facility provides several levels of logs. Please use it as described below.

| Level           | Name                 | Use                                          |
|-----------------|----------------------|----------------------------------------------|
| Developer debug | log_level::dev_debug | Developer debugging information              |
| User debug      | log_level::debug     | User debugging information                   |
| Notice          | log_level::notice    | Informational notices                        |
| Warning         | log_level::warning   | Behaviour that does not terminate execution. |
| Error           | log_level::error     | Behaviour that terminates execution.         |

By default, the executable is compiled to include only `log_level::notice` messages or higher, so only `notice`, `warning`, and `error` are recorded. This behaviour can be changed at compilation by manually changing the `LOG_THRESHOLD` setting in `aa_log.hpp` to one of the named `log_level`s defined in the table above (default is `log_level::notice`). All logs of any kind can also be completely silenced at the build stage by providing CMake the `-DLOGGING_FACILITY=OFF` flag (default is ON).

To use logging in your code, include the logging class in your code

```
#include "aa_log.hpp"
```

Here are several examples of how to log a message

```
LOG(log_level::debug, "A user debug message");
LOG(log_level::notice, "A notice message");
LOG(log_level::warning, "A warning message");
LOG(log_level::error, "An error message");
```

By default, all messages are printed to the console. To save them to the disk instead, write:

```
FILE* pFile = fopen("/full/path/to/destination/logfile.log", "a");
FILElog::stream() = pFile;
```

All subsequent logs will then only be written to disk.

All logs will print the system date and time, the log level, and the log message.
