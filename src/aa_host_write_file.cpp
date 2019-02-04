#include <stdio.h>
#include <stdlib.h>

namespace astroaccelerate {

  /**
   * \brief This function writes the transformed space (dm,t) out to file in binary format.
   */
  void write_output(int i, int t_processed, int ndms, float *output_buffer, float *dm_low, float *dm_high) {
    FILE *fp_out;
    char filename[200];
    
    sprintf(filename, "dm_%.2f-%.2f.dat", dm_low[i], dm_high[i]);
    if (( fp_out = fopen(filename, "w") ) == NULL) {
      fprintf(stderr, "Error opening output file!\n");
      return;
    }

    for (int k = 0; k < ndms; k++) {
      for (int l = 0; l < t_processed; l++) {
	fprintf(fp_out, "%f ", output_buffer[k * t_processed + l]);
      }
      fprintf(fp_out, "\n");
    }
  }
} //namespace astroaccelerate
