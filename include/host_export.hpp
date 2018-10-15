#ifndef ASTRO_ACCELERATE_HOST_EXPORT_HPP
#define ASTRO_ACCELERATE_HOST_EXPORT_HPP

extern size_t Calculate_sd_per_file_from_file_size(size_t size_in_MB,
                                                   size_t primary_dimension,
                                                   size_t floats_per_elements);

extern void Export_DD_data(int         range,
                           float ***   output_buffer,
                           size_t      max_nTimesamples,
                           int *       ndms,
                           int *       inBin,
                           float *     dm_low,
                           float *     dm_high,
                           float *     dm_step,
                           const char *base_filename,
                           int *       ranges_to_export,
                           int         DMs_per_file);

extern void Export_data_raw(float *     data,
                            size_t      primary_dimension,
                            size_t      secondary_dimension,
                            const char *base_filename,
                            int         sd_per_file);

extern void Export_data_as_list(float *     h_data,
                                size_t      primary_dimension,
                                float       prim_mul,
                                float       prim_add,
                                size_t      secondary_dimension,
                                float       sec_mul,
                                float       sec_add,
                                const char *base_filename,
                                int         sd_per_file);

#endif
