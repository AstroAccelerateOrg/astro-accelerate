Developers Guide
=====

Dependencies
===
Astro-Accelerate depends on the CUDA libraries, and in particular the following libraries
   * cuFFT
   * cuRAND

The minimum device architecture supported is `3.5` and above.

Do not introduce any new dependencies on any other product without consultation.

Repository Conventions
===

### Solving Issues
* Create a new issue on GitHub. Fill in the issue template. Assign the issue to someone.
* Checkout a new branch named `abc_##_description` where `abc` are your initials, `##` is the issue number assigned by GitHub, and `description` is a very brief description of the branch.
* Do not use special characters and in particular no slash characters when naming branches.

### Committing and Pushing
* Work and commit incrementally.
* Commit messages to include an extended description, nominally 100 characters long or longer are permissible to provide detail.
* Push the changes regularly.
   
### Pull Requests
   * When a branch is ready to commit, create a pull request and use the pull request template.
   * Fill in the details of the pull request template. In particular, tag which issues are being fixed by using `fixes ##`, where `##` is the issue number of any issues being fixed. If multiple issues are being fixed, separate the issues with a comma, adding `fixes` for each one of them.
   * Add a meaningful description of the changes. For instance, detail each file that was changed and why it was changed. Provide simple instructions to test whether the changes work.
   * Assign the pull request to a relevant person (either the person who created the issue, or someone who has permissions to merge the branch). Tag relevant people for review, and/or tag people in the pull request description for commenting as needed.
   * Use `Draft:` and/or `[Draft]` and/or mark a merge request as a draft if the branch should not be merged yet.
   * Once the branch is merged, the issue it relates to is either closed automatically if it was tagged, or it will be manually closed by the person who merged the branch.
   * The branch may only be deleted if the pull request has indicated this explicitly. If this is the case, the branch will be deleted by the person who merged the branch.
   
### Release Tagging
   * When a new tagged release is to be made, a new issue should be created to request it.
   * This issue should be tagged with appropriate labels, and assigned to someone who can create tagged releases.
   * The tagged release will be created, and the person who did this will then close the issue.
   * The reasons for requesting a new tagged release should be added to the issue description.
   * The person creating the new tag will inspect the pull requests since the last tag, and assess what the semantic version number increment should be for the new tagged release.
   * A brief description for the tagged release should be added. Suggestions can be made in the issue description for the issue that requested the tag, but the previous pull requests should be taken into consideration when naming and describing the new tagged release.
   
### Forking
   * Changes from forked repositories will not be merged into the main repository code.


### `master` Branch Guidelines

* The `master` branch is protected. No one shall commit to it directly.
* The `master` branch should always be stable and ready to be released and tagged.
* The `master` branch is tightly controlled.
* The `master` branch shall compile and pass all tests at all times.
* Pull requests to the `master` branch must only be done after code review.
* Pull requests to the `master` branch must be done from a working development branch that passes all tests.
* Developers should keep their own branch up to date with the branch from which they are working.
* Branches that introduce new features should only be merged if they are up to date with the development branch from which they are branched.
* Commits to a development branch for which the intention is to merge to the `master` branch should be code reviewed and also function properly before being merged.
* Any failure to adhere to the above should be remedied with high priority.

Repository Structure Astro-Accelerate
===
The following is an outline of the file structure of the repository.
When the repository structure is changed, this document should be updated to reflect it.
* `/`
    * Top-level files
    * Setup scripts
    * Basic documentation
    * `.github/`
        * Templates for issues, bugs, and pull requests.
    * `input_files/`
        * Sample input configuration files.
    * `lib/`
        * Source code of the project.
    * `scripts/`
        * Script to run the Astro-Accelerate routines

Coding Conventions
===
Please adhere to the following where possible. Exceptions to coding conventions should be justified as appropriate.

### Naming Conventions
* Methods and functions shall be all lower case, with word seperation with an underscore character e.g. `do_something()`.
* Variable names shall be all lower case.
* File names shall follow this naming convention: aa_<optional: device/host>_<mandatory: name>.<mandatory: cpp/hpp/cu/cuh>.
* Header guards shall follow this naming convention: ASTRO_ACCELERATE_<full_file_name_with_extension>.

### Use of Indentation
* Use 2 spaces instead of tab.
* Namespaces shall not be indented

### Use of Carriage Return, Line Feed
* Use LF in preference to CRLF

Documentation
===
* Add Doxygen style comments for each method and each parameter in every header file.
* Where relevant, comments should provide example usage.
* Document the bounds of any parameter
* Document untested use cases/parameters (i.e. where you have not provided an explicit unit test).
* Add comments as appropriate in the `.cpp` file and non-public headers.

Using CUDA
===
If CUDA is available and it has been explicitly activated in the build system, then the `ENABLE_CUDA` flag will be set. Use this flag to ensure code compiles and that code tests are run with CUDA as appropriate.

# Build System
The build system uses CMake and is configured in the top-level file `CMakeLists.txt`.
## How to add a new architecture setting
To add a new architecture setting

1. Open `CMakeLists.txt`
2. Find the line that says `if(ARCH MATCHES ALL|[Aa]ll)`, and add a new line to append to the list, in the same way as the existing ones: `list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_70,code=sm_70)` - then also add it to the comma separated list for `ASTRO_ACCELERATE_CUDA_ARCH_VERSION` and `ASTRO_ACCELERATE_CUDA_SM_VERSION`.

3. Then, add the setting to a new `elseif` statement for `elseif(ARCH MATCHES X.Y)` (where `X.Y` are integers reflecting the architecture major.minor version being added) in the same way as the existing architecture settings.

The architecture setting is used when running CMake by adding a flag `-DCUDA_ARCH="X.Y"` (where `X.Y` are substituted with integers reflecting the architecture major.minor version).
## How to add a new unit test
It is good practice that developers and users create and run unit tests as part of regular project development.

To add a unit test:

1. From the top-level directory, go into the `tests` folder.
2. Add a new source file with the `.cpp` extension and give it a sensible name (e.g. `test_a_new_feature` in this example).
3. Write the unit test in the source file.
4. Copy one of the `.cmake` files and rename it to the same name as the source file, keeping the `.cmake` extension.
5. Edit the newly created `.cmake` file and set the `TEST_NAME` to the same name as the file name but exclude the `.cmake` extension.
6. Under `target_sources`, the `PRIVATE` sources and `PUBLIC` headers can be edited as per the existing examples.
7. The `set_tests_properties` can be edited as per CMake instructions. Several of the existing examples have a `"Runs"` string as the `PASS_REGULAR_EXPRESSION`, which means the test passes if the test executable will print `Runs`, and fails otherwise. More information about [`set_tests_properties` can be found here](https://cmake.org/cmake/help/v3.0/command/set_tests_properties.html).
8. Edit the `CMakeLists.txt` file located inside the `tests` directory (`tests/CMakeLists.txt`) and include the newly created `.cmake` file `include(tests/test_a_new_feature.cmake)`

## How to add a new example
AstroAccelerate includes several example codes, which help new users to develop their applications more quickly by learning by example.

To add a new example:

1. From the top-level directory, go into the `examples/src` folder.
2. Add a new source file with the `.cpp` extension and give it a sensible name (e.g. `my_example.cpp` in this example).
3. Write the example code.
4. Edit the `CMakeLists.txt` file located inside the `examples` directory (`examples/CMakeLists.txt`) and include the newly created `.cpp` file as an executable `add_executable(my_example src/my_example.cpp)` and add AstroAccelerate as the target link library to the executable `target_link_libraries(my_example astroaccelerate)`

## How to add a new field in `aa_version.hpp`
`aa_version.hpp` is a file that is automatically generated by the `CMake` build system from a template file.

From the top-level directory, the template file is located in `cmake/version.h.in`.

The file can be edited just like any file, and the contents will be copied into the `aa_version.hpp` when running `cmake`.

Ordinary text is copied verbatim, whereas a variable (variable names are surrounded by the @ symbol, e.g. @my_variable@) configured via CMake will be substituted by CMake.

CMake knows about the `version.h.in` file via the line in `CMakeLists.txt` that says `configure_file("${PROJECT_SOURCE_DIR}/cmake/version.h.in" "${PROJECT_SOURCE_DIR}/include/aa_version.hpp")`.

# Testing
## How to run unit tests
The compilation of unit tests is enabled/disabled when running CMake by adding a flag `-DENABLE_TESTS=ON` (to enable unit tests) or `-DENABLE_TESTS=OFF` (to disable unit tests).

To run all unit tests, write `make test`, the tests will run and the console will print the test results to the screen, and to a number of log files in the build directory.

## How to run examples
When compiling and building via `CMake`, all example codes are built automatically and their executables will be placed inside a folder called `examples` which will be inside the build directory.

# Developing new components
To indicate the availability of a new component to the user, add it to `enum class components` in `include/aa_pipeline.hpp`, and likewise for `component_option`s.

The component must allow configuration via a `plan`, a `strategy`, and be implemented in a `permitted_pipeline`. Finally, the component must be provisioned via `aa_config` (to read the new settings from an `input_file`) and `aa_pipeline_api` so that the `permitted_pipeline` can be run.

## How to add a new plan
Name the class, appending `_plan` to the file name. Provision the class with a trivial constructor that sets all values to 0 (or similar). Add another constructor that takes the necessary parameters. Ideally, the class only offers getter methods, so that all configuration is done at construction.

## How to add a new strategy
Name the class, appending `_strategy` to the file name.
Provision the class with a trivial constructor that sets all values to 0 (or similar). Also implement a constructor that takes a `plan` as an input parameter, and sets its values based on the `plan` object. Remember to adhere to the `aa_strategy` interface defined in the base class. Ideally, the class only offers getter methods, so that all configuration is done at construction.

## How to add a new permitted_pipeline
1. Create a new `permitted_pipeline` from one of the existing templates. Add (templated) constructors like for the existing examples. The constructors should include only the required and `strategy` objects, and any other settings required for configuration at construction. Remember to adhere to the `aa_pipeline_runner` interface, so that the class offers the expected functionality for the user. Additional functionality can optionally be added. The `permitted_pipeline` is responsible for its own cleanup, and the user may not clean up their input data themselves, so plan resource allocation and de-allocation accordingly.

2. Once the `permitted_pipeline` exists, it must be added to the `permitted_pipelines` in `aa_permitted_pipelines.hpp`. Do this by adding a new `static const aa_pipeline::pipeline` declaration (which you implement in `aa_permitted_pipelines.cpp`), and add a statement that returns `true` for the `is_permitted` method.

3. In the `aa_permitted_pipelines.cpp` source file, implement the new pipeline by listing the components that it will run. If the component runs in combination with other components, create additional pipelines for each combination.

4. The `permitted_pipeline` must still be provisioned in `aa_pipeline_api.hpp`, and any new module components must be added to `pipeline::component` in `aa_pipeline.hpp`. This is done via a lookup in the `ready()` method. The pipeline instance must be created via a `std::unique_ptr` instance of the same type as the class name of the newly created `permitted_pipeline` class, and then it must be assigned to `m_runner` (which will run the pipeline by calling the pipeline's `run` method). The lookup itself is well-documented and contains plenty of examples showing how this is implemented for each `permitted_pipeline` with `if` and `else` statements.

5. Add `bind` calls to the `aa_pipeline_api` class so that a `plan` object can be bound for this new component. Similarly, add a method to return the `strategy` object. The existing code shows how to implement this functionality. Ensure that the ready state of the pipeline is falsified whenever a new `bind` call occurs. Equally, when a valid strategy is obtained, it must be returned.

At this stage, a library user has everything they need to use the component in their own application code. However, it has not yet been provisioned for the standalone code. To do this, the new component must be configurable via an `input_file`, and any user flags that can be used to configure it must be offered via the `aa_config` class.

1. Add any configuration flags for the component to `struct aa_config_flags`, and in the `get_user_input()` method implement that the `input_file` is checked if the user provided them in the `input_file`.
2. The component must now be added to `src/aa_main.cpp` which implements the standalone executable functionality. Create a plan object for the new component, passing the flags from the `aa_config` that read the user input. Surround the `bind` call by a simple check to see if the plan is required, for example `if(pipeline.find(aa_pipeline::component::fdas) != pipeline.end())` - this makes sure the bind call only happens if the user actually requested a particular component.
