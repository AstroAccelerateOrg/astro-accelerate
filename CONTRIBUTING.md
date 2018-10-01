Developers Guide
=====

Dependencies
===
Astro-Accelerate depends on the CUDA libraries, and in particular the following libraries
   * CUFFT
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
   * Prepend `[WIP]` and/or `[DO NOT MERGE]` if the branch is `Work In Progress` or if the branch should not be merged yet. Remove these flags when no longer needed.
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
