= Developers Guide

== Dependencies Do not introduce any new dependencies on any other product without consultation

= Repository == Workflow and Branches - Never work directly on the master branch - this is the "stable release" and is tightly controlled - The master shall compile and pass all unit tests at all times. - Commits to the master must only be done after code review and from a working dev branch - devs should keep their own branch in sync with the dev branch. - The dev should compile and pass all unit tests. Any build breaks should receive top priority. - Only merge from a feature branch that is already synced with the dev branch. - Commits to the dev should usually only be done after code review and from a working feature branch

== Repository Structure AstroAccelerate 
+ + astrolib Archive of the library containing cuda modules / main routines 
	+ doc Documentation, contains developers guide 
	+ input_files Input files 
	+ lib Source code of the project 
	+ scripts Script to run the AstroAccelerate routine 
	
= Coding Conventions Coding conventions are used in this project to ensure conformity in style and to help reduce any of the very many pitfalls c++ programming can introduce, Remember that these are just general guidelines. Try to stick to them as much as possible, but there will always be cases where doing it differently is justified. Document these exceptions with comments in the code.

== Naming Conventions
    method and functions shall be all lower case, with word seperation with an "_" e.g. do_something()
    variable names shall be all lower case

== Indendation & CRLF of source files

    use 2 spaces instead of tab.
    use LF in preference to CRLF
    namespaces shall not be indented

== Documentation

    Add doxygen style comments for each method and each parameter in every public header file.
    Think about adding an example of the use of each method
    Document the bounds of any parameter
    Document untested use cases/parameters (i.e. where you have not provided an explicit unit test)
    Add comments as appropriate in the cpp file and non-public headers.

== Using CUDA If CUDA is available and it has been explicitly activated in the build system, the ENABLE_CUDA flag will be set. Use this flag to ensure your code compiles and your unit tests run with or without CUDA availability.

