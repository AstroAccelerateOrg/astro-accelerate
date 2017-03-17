= Developers Guide

== Dependencies
Do not introduce any new dependencies on any other product without consultation

= Repository
== Workflow and Branches
    - Never work directly on the master branch - this is the "stable release" and is tightly controlled
        - The master shall compile and pass all unit tests at all times.
        - Commits to the master must only be done after code review and from a working dev branch
    - devs should keep their own branch in sync with the dev branch.
        - The dev should compile and pass all unit tests. Any build breaks should receive top priority.
    - Only merge from a feature branch that is already synced with the dev branch.
        - Commits to the dev should usually only be done after code review and from a working feature branch
        
== Repository Structure
      AstroAccelerate +
                      + astrolib        static library, contains cuda kernels
                      + cmake           Build system
                      + doc             Documentation, contains developers guide
                      + input_files     Input files, user can use it to  define which 					features he wants to use
                      + lib             Source code of the project.
                      + scripts         Script to run the AstroAccelerate routine
                      + thirdparty      
                        + googletest    Google unit test
                      + tools           Useful tools like a class generator

= The Build System
This project uses cmake as its build system. See the notes in the CMakeLists.txt file in the top directory

== Adding a src file/class
    - You can use the class tool (in the tools directory) to autogenerate a suitable class template
    - Each time you add a new class, in order for it to build you will need to add a reference to it in the appropriate CMakeLists.txt file (either in AstroAccelerate/module_name/CMakeLists.txt or AstroAccelerate/module_name/test_utils/CMakeLists.txt)
    - place source names in alphabetical order when adding to variables

= Coding Conventions
Coding conventions are used in this project to ensure conformity in style and to help
reduce any of the very many pitfalls c++ programming can introduce,
Remember that these are just general guidelines. Try to stick to them as much as possible,
but there will always be cases where doing it differently is justified. Document these exceptions with comments in the code.

== Namespaces
- all classes shall be in the AstroAccelerate::<module_name> or AstroAccelerate::<module_name>::test namespace as appropriate.
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

- All classes must be in the astroaccelerate namespace
- Header files must start with.
'''
#ifdef ASTROACCELERATE_<CLASS_NAME>_H
#define ASTROACCELERATE_<CLASS_NAME>_H
... body of header
#endif // ASTROACCELERATE_<CLASS_NAME>_H
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
- indent private, protected, and public keywords by 4 spaces
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
- Use std::exceptions as appropriate
- If you need your own exceptions for some reason, try and avoid any memory allocation of the exception. i.e follow the std::exception model

== Unit Tests
- gtest shall be used as the unit test framework
- each class should have its own unit test suite which shall be in files <CLASS_NAME>Test.cpp and <CLASS_NAME>Test.h in the unit test directory
- any common functionality should be considered for inclusion in the test utility library, and provided with appropriate headers and documentation
- all unit tests should be verified to pass before making a merge request
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
Use this flag to ensure your code compiles and your unit tests run with or without CUDA availability.

Setting specific compile flags for subpackage products shall be done through the cuda_subpackage_compile macro:
e.g.
cuda_subpackage_compile(${lib_src_cuda} OPTIONS "-arch compute_35")

== Extending std classes
While this style guide requests CamelCase typedefs, there are cases where the standard convention is to use lower_case
names for typedefs. This is true in the case of classes that extend or can be used in place of classes from the
standard library (or as in the case of thrust, boost, etc. libraries that follow conventions from the standard library) that conform to specific concepts. Places where this is done should be clearly documented.