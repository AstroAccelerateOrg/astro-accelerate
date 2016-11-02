= Developers Guide

== Dependencies
Do not introduce any new dependencies on any other product without consultation

== Build Products
The project will build the pipeline executable, two c++ libraries, and their associated header files.
One library will support the main product and the second a library is a set of utilitites to aid in unit testing code
that uses the main library.

= Repository
== Workflow and Branches
    - Use the gitlab workflow (https://repo.oerc.ox.ac.uk/help/workflow) as much as possible.
    - Never work directly on the master branch - this is the "stable release" and is tightly controlled
        - The master shall compile and pass all unit tests at all times.
        - Commits to the master must only be done after code review and from a working dev branch
    - devs should keep their own branch in sync with the dev branch.
        - The dev shuold compile and pass all unit tests. Any build breaks should receive top priority.
    - Only merge from a feature branch that is already synced with the dev branch.
        - Commits to the dev should usually only be done after code review and from a working feature branch
        - use the gitlab merge request feature for feature->dev branch merges

== Repository Structure
cmake          : project specific files for configuring the build system
tools          : useful things to aid the developer such as a class template generator
doc            : documentation & doc templates, with the exception of the API doc of the classess. Also see module doc directoy.
cheetah        : this is the src code for the cheetah library. It is organised as a set of modules with
                 the following minimal structure (following the boost conventions for templates)

                 cheetah +
                         + <module_name> +                      In this directory only put header files that you wish exposed to module users
                                         + doc                  Add additional doxygen documentation for tutorials etc here.
                                         + examples             src code example uses
                                         + src                  All .cpp files for the library and internal headers
                                         + detail               All .cpp implementations details for template classes that cannot be precompiled
                                         + test_utils +         Include files for the supporting test library
                                                      + src     Implementation details (.cpp files) for the test library
                                                      + detail  Template implementation details for the test support library
                                         + test                 The unit test directory. Header files for unit tests live here
                                                      + src     Implementation details (.cpp files) for unit tests


cheetah/test     : the location of tests and any public headers for the cheetah integration tests
cheetah/test/src : the location of all test implementation coded

= The Build System
This project uses cmake as is build system. See the notes in the CMakeLists.txt file in the top directory

== Adding a src file/class
    - Each time you add a new class, in order for it to build you will need to add a reference to it in the appropriate CMakeLists.txt file
      (either in cheetah/module_name/CMakeLists.txt or cheetah/module_name/test_utils/CMakeLists.txt)
      - place source names in alphabetical order when adding to variables
    - You can use the class tool (in the tools directory) to autogenerate a suitable class template

== Adding a submodule
    - To add a submodule use the SUBPACKAGE macro (see doc in cmake/subpackage.cmake).
      To simplify, the purpose of this macro is to add the submodules contribution to the lib_src_${arch} variables which
      are used to define the contents of the cheetah libraries.
      (n.b ensure the variables defined in the submodule use PARENT_SCOPE qualifier)

== Including a dependency in the build system
    - unless essential (and approved to be so) make the dependency optional, and off by default. Allow the ENABLE_<feature> flag to enable.
    - add a new file cmake file in the cmake directory for that product/feature and include this from cmake/dependencies.cmake
    - Stick to a standard convention for naming the product location.
      TODO - decide what these should be
      i.e. <NAME>_ROOT, <NAME>_INCLUDE_DIRECTORIES, <NAME>_LIBRARIES, FOUND_<NAME>
    - Test these can be overriden succesfully from the command line and non-system versions can be specified easily


= Coding Conventions
Coding conventions are used in this project to ensure conformity in style and to help
reduce any of the very many pitfalls c++ programming can introduce,
Remember that these are just general guidelines. Try and stick to them as much as possible,
but there will always be cases where doing it differently is justified. Document these exceptions with comments in the code.

== Namespaces
- all classes shall be in the ska::cheetah::<module_name> or ska::cheetah::<module_name>::test namespace as appropriate.
- avoid using where possible "using namespace ..". Avoid completely in public header files.

== Naming Conventions
- class names shall start with a capital letter and be CamelCased
- method and functions shall be all lower case, with word seperation with an "_" e.g. do_something()
- variable and member names shall be all lower case
- all members shall be prefixed with an underscore "_"

== Indendation & CRLF of source files
- use 4 spaces instead of tab.
- use LF in preference to CRLF
- namespaces shall not be indented

== C++ Files
Each class shall have its own separate cpp and header file. The name of the file should be the same as the class with the .h or .cpp extension
Do not mix two classes together in the same file.
An exception to this rule can be made where there is a private Internal Implementation class unique to the
class being implemented AND it is not confusing to do so.

- All classes must be in the ska, cheetah namespace
- Header files must start with.
'''
#ifdef SKA_CHEETAH_<MODULE_NAME>_<CLASS_NAME>_H
#define SKA_CHEETAH_<MODULE_NAME>_<CLASS_NAME>_H
... body of header
#endif // SKA_CHEETAH_<MODULE_NAME>_<CLASS_NAME>_H
'''
  In the tool directory you will find a class template generator that you may find helpful.

== Documentation
- Add doxygen style comments for each method and each parameter in every public header file.
- Think about adding an example of the use of each method
- Document the bounds of any parameter
- Document untested use cases/parameters (i.e. where you have not provided an explicit unit test)
- Add comments as appropriate in the cpp file and non-public headers.

== Constructors
- Never use new inside a constructor parameter list e.g. dodgy = MyClass(new BadThingToDo());.
  This could easily lead to a memory leak if the constructor of MyClass should fail.
- Do not create long lists of parameters with the same types. Use strong typeing as much as possible to
  get your compiler to do the debugging for you.
    e.g. thing(int x_dimension, int y_dimension, int subject_index)
         could become:
         thing( Dimension x, Dimension y, SubjectIndex i );
         or even:
         thing( DimensionX x, DimensionY y, SubjectIndex i );
- align commas in the initialiser list on seperate lines
    e.g.
      MyClass(A a, B b)
         : _a(a)  // indent : by 4 spaces on its own line
         , _b(b)  // align comma under the :
      {
      }


== Class Declarations
- all members shall be prefixed with an underscore "_"
- all members shall be made either private or protected as appropriate
- indent private, protected, amd public keywords by 4 spaces
- always use the override keyword when you inherit from a method
- do not mix members, typdefs, or functions together, but keep them in their own public:, protected:, or private: blocks.
- all public methods and typedef should be documented
- header shall be ordered as much as possible as follows, each type in their own public/private/protected block:
  e.g.

  class MyClass
  {
      public:
        // typdefs

      public:
        // public methods before protected

      protected:
        // protected methods before private

      private:
        // private methods

      protected:
        // protected members

      private:
        // private members
  };

== Exceptions
- Prefer exceptions, but use return codes where it is more appropriate
- No exceptions to be propagated outside of its thread - think about passing it along with a std::exception_ptr.
- Use std::exceptions as appropriate, or panda::Error if you want to add messages with runtime information.
- If you need your own exceptions for some reason, try and avoid any memory allocation of the exception. i.e follow the std::exception model
  It should also be OK to inherit from panda::Error.

== Unit Tests
- gtest shall be used as the unit test framework
- each class should have its own unit test suite which shall be in files <CLASS_NAME>Test.cpp and <CLASS_NAME>Test.h in the unit test directory
- any common functionality should be considered for inclusion in the test utility library, and provded with appropriate headers and documentation
- all unit tests should be verified to pass before making a merge reuest
- unit test names shall follow the function naming convention of the c++ irrespective of whether it is implemented as a class or not
- unit test groups shall be CamelCased
  e.g TEST_F(MyTestGroup, test_something_really_well)
- to run all the tests
'''
make test
'''
- to run all the tests and generate a test coverage report (only in 'coverage' type builds)
'''
make coverage_report_test
'''

== Using CUDA
If CUDA is available and it has been explicitly activated in the build system, the ENABLE_CUDA flag will be set.
Use this flag to ensure your code compiles wnd your unit tests run with or without CUDA availability.

Setting specific compile flags for subpackage products shall be done through the cuda_subpackage_compile macro:
e.g.

cuda_subpackage_compile(${lib_src_cuda} OPTIONS "-arch compute_35")

== Extending std classes
While this style guide requests CamelCase typedefs, there are cases where the standard convention is to use lower_case
names for typedefs. This is true in the case of classes that extend or can be used in place of classes from the
standard library (or as in the case of thrust, boost, etc. libraries that follow conventions from the standard library) that conform to specific concepts. Places where this is done should be clearly documented.
