#ifndef ASTRO_ACCELERATE_AA_FDAS_UTIL_HPP
#define ASTRO_ACCELERATE_AA_FDAS_UTIL_HPP

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <getopt.h>
#include <string.h>
#include <sys/stat.h>

namespace astroaccelerate {

  /**
   * \struct cmd_args
   * \brief Struct that contains command line arguments for configuring fdas and is used by AstroAccelerate internally.
   * \brief The user should not use this struct for configuring fdas. Instead the user should use aa_fdas_plan and aa_fdas_strategy.
   */
  typedef struct {
    int search;
    int wsig;
    int basic;
    int kfft;
    char * afname;
    char * nfname;
    float freq0;
    float thresh;
    int mul;
    int zval;
    int writef;
    float duty;
    int iter;
    int nsig;
    int nharms;
    int devid;
    int inbin;
    int norm;
    float sigamp;
  } cmd_args;

  //Function declarations
  extern "C" {
    void print_usage(); 

    void read_command_line(int argc, char *argv[], cmd_args *args);

    int read_input_file(char *fname, float **array);

    void read_input_file_cpx(char *fname, float2 *array);

    void write_output_file(char *fname, float **array, unsigned long len);
  }

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_FDAS_UTIL_HPP
