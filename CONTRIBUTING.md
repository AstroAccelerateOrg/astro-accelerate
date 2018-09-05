Developers Guide
=====

Dependencies
===
Do not introduce any new dependencies on any other product without consultation.


Repository `master` Branch Guidelines
===
### Workflow and Branches

* The `master` branch should always be stable and ready to be released and tagged.
* The `master` branch is tightly controlled.
* The `master` branch shall compile and pass all tests at all times.
* Commits to the `master` branch must only be done after code review.
* Commits to the `master` branch must be done from a working development branch that passes all tests.
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
