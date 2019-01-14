#include <cstdio>
#include <fstream>
#include <vector>

namespace astroaccelerate {

  void Export_data_raw(float *h_data, size_t primary_dimension, size_t secondary_dimension, const char *base_filename, int sd_per_file) {
    char final_filename[200];
    size_t inner_sec_dim_shift;
    if(sd_per_file<0) sd_per_file = secondary_dimension;
    size_t export_size = primary_dimension*sd_per_file;
	
    float *h_export_data;
	
    int nRepeats = secondary_dimension/sd_per_file;
    int nRest = secondary_dimension%sd_per_file;
    std::vector<int> chunk_size;
    for(int f=0; f<nRepeats; f++) chunk_size.push_back(sd_per_file);
    if(nRest>0) chunk_size.push_back(nRest);
    printf("  Data will be exported into %d files\n", (int) chunk_size.size());
    float percent_per_chunk = 100.0/(secondary_dimension);

    h_export_data = new float[export_size];

    sprintf(final_filename,"%s_info.txt", base_filename);
    std::ofstream FILEOUT;
    FILEOUT.open (final_filename, std::ofstream::out);
    FILEOUT << primary_dimension << std::endl;
    FILEOUT << secondary_dimension << std::endl;
    FILEOUT << sd_per_file << std::endl;
    FILEOUT << (int) chunk_size.size() << std::endl;
    FILEOUT.close();	

    inner_sec_dim_shift = 0;
    for(int i=0; i<(int) chunk_size.size(); i++){
      sprintf(final_filename,"%s_%d.dat", base_filename, i);
		
      for(int d=0; d<chunk_size[i]; d++) {
	for(size_t t=0; t<primary_dimension; t++) {
	  h_export_data[d*primary_dimension + t] = h_data[(inner_sec_dim_shift + d)*primary_dimension + t];
	}
      }
		
      FILE *fp_out;
      if (( fp_out = fopen(final_filename, "wb") ) == NULL) {
	fprintf(stderr, "Error opening output file!\n");
	exit(0);
      }
      fwrite(h_export_data, primary_dimension*chunk_size[i]*sizeof(float), 1, fp_out);
      fclose(fp_out);
		
      inner_sec_dim_shift = inner_sec_dim_shift + chunk_size[i];
      printf("    Exported: %0.2f%%\n", ((float) inner_sec_dim_shift)*percent_per_chunk);
    }

    delete [] h_export_data;
  }


  void Export_data_as_list(float *h_data, size_t primary_dimension, float prim_mul, float prim_add, size_t secondary_dimension, float sec_mul, float sec_add, const char *base_filename, int sd_per_file) {
    char final_filename[200];
    size_t inner_sec_dim_shift;
    if(sd_per_file<0) sd_per_file = secondary_dimension;
    size_t export_size = primary_dimension*sd_per_file;
	
    float *h_data_to_export;
	
    int nRepeats = secondary_dimension/sd_per_file;
    int nRest = secondary_dimension%sd_per_file;
    std::vector<int> chunk_size;
    for(int f=0; f<nRepeats; f++) chunk_size.push_back(sd_per_file);
    if(nRest>0) chunk_size.push_back(nRest);
    printf("  Data will be exported into %d files\n", (int) chunk_size.size());
    float percent_per_chunk = 100.0/(secondary_dimension);
	
    h_data_to_export = new float[export_size*3];
	
    inner_sec_dim_shift = 0;
    for(int i=0; i<(int) chunk_size.size(); i++){
      sprintf(final_filename,"%s_%d.dat", base_filename, i);
		
      for(int sec=0; sec<chunk_size[i]; sec++) {
	for(size_t prim=0; prim<primary_dimension; prim++) {
	  h_data_to_export[3*(sec*primary_dimension + prim)] = ((float) sec)*sec_mul + sec_add;
	  h_data_to_export[3*(sec*primary_dimension + prim)+1] = ((float) prim)*prim_mul + prim_add;
	  h_data_to_export[3*(sec*primary_dimension + prim)+2] = h_data[sec*primary_dimension + prim];
	}
      }
		
      FILE *fp_out;
      if (( fp_out = fopen(final_filename, "wb") ) == NULL) {
	fprintf(stderr, "Error opening output file!\n");
	exit(0);
      }
      fwrite(h_data_to_export, primary_dimension*chunk_size[i]*sizeof(float), 3, fp_out);
      fclose(fp_out);
		
      inner_sec_dim_shift = inner_sec_dim_shift + chunk_size[i];
      printf("    Exported: %0.2f%%\n", (float) (inner_sec_dim_shift*percent_per_chunk));
    }
	
    delete [] h_data_to_export;
  }

} //namespace astroaccelerate
