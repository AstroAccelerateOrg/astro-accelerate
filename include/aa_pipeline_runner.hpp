#ifndef ASTRO_ACCELERATE_AA_PIPELINE_RUNNER_HPP
#define ASTRO_ACCELERATE_AA_PIPELINE_RUNNER_HPP

#include <iostream>
#include <vector>
#include <fstream>

#include "aa_log.hpp"
#include "aa_device_analysis.hpp"


namespace astroaccelerate {

	/**
	 * \class aa_pipeline_runner aa_pipeline_runner.hpp "include/aa_pipeline_runner.hpp"
	 * \brief Abstract base class for running a pipeline.
	 * \details In practice, this class is used for creating base class pointers for pointing to a permitted pipeline instance.
	 * \author AstroAccelerate team.
	 * \date 03 July 2019.
	 */
	class aa_pipeline_runner {
	public:

		enum class status : int {
			error = -1,
			finished = 0,
			has_more = 1,
			finished_component = 2
		};

		/** \brief Virtual destructor for aa_pipeline_runner. */
		virtual ~aa_pipeline_runner() {

		}

		/** \brief Virtual setup method to be implemented by derived classes. */
		virtual bool setup() = 0;

		/**
		 * \brief Base class virtual methods for running a pipeline.
		 * \details In case a derived class does not implement a method, this method will be called.
		 */
		virtual bool next() {
			// If a derived class does not implement this method, this method is used.
			std::cout << "ERROR:  The selected operation is not supported on this pipeline." << std::endl;
			return false;
		}

		/**
		 * \brief Base class virutal methods for running a pipeline.
		 * \details In case a derived class does not implement a method, this method will be called.
		 * \details This implementation returns a status code.
		 */
		virtual bool next(aa_pipeline_runner::status &status_code) {
			LOG(log_level::error, "The selected operation is not supported on this pipeline.");
			status_code = aa_pipeline_runner::status::finished;
			return false;
		}

                virtual bool next(std::vector<analysis_output> &, aa_pipeline_runner::status &status_code) {
                        LOG(log_level::error, "The selected operation is not supported on this pipeline.");
                        status_code = aa_pipeline_runner::status::finished;
                        return false;
                }

		/**
		 * \brief Base class virtual methods for running a pipeline.
		 * \details In case a derived class does not implement a method, this method will be called.
		 */
		virtual bool next(std::vector<std::vector<float>> &, int &, std::vector<int> &) {
			// If a derived class does not implement this method, this method is used.
			std::cout << "ERROR:  The selected operation is not supported on this pipeline." << std::endl;
			return false;
		}

		virtual float ***output_buffer(){
			LOG(log_level::error, "The selected operation is not supported on this pipeline (output_buffer).");
			return NULL;
		}

		virtual float *h_SPD_snr(){
			LOG(log_level::error, "The selected operation is not supported on this pipeline (h_SPD_snr).");
			return NULL;
		}

		virtual unsigned int* h_SPD_ts(){
			LOG(log_level::error, "The selected operation is not supported on this pipeline (h_SPD_time).");
			return NULL;
		}

		virtual unsigned int* h_SPD_dm(){
			LOG(log_level::error, "The selected operation is not supported on this pipeline (h_SPD_dm).");
			return NULL;
		}

		virtual unsigned int* h_SPD_width(){
			LOG(log_level::error, "The selected operation is not supported on this pipeline (h_SPD_width).");
			return NULL;
		}

		virtual size_t get_SPD_nCandidates(){
			LOG(log_level::error, "The selected operation is not supported on this pipeline (nCandidates).");
			return 0;
		}

		virtual int get_current_range(){
			LOG(log_level::error, "The selected operation is not supported on this pipeline (current_range).");
			return 0;
		}

                virtual int get_current_tchunk(){
                        LOG(log_level::error, "The selected operation is not supported on this pipeline (current_time_chunk).");
                        return 0;
                }

		virtual long int get_current_inc(){
			LOG(log_level::error, "The selected operation is not supported on this pipeline (current_inc).");
			return 0;
		}

		virtual int total_computed_samples(){
			LOG(log_level::error, "The selected operation is not supported on this pipeline (total_computed_samples).");
			return 0;
		}

                virtual bool cleanup(){
                        LOG(log_level::error, "The selected operation is not supported on this pipeline (cleanup).");
                        return 0;
                }

		
	protected:
		/** \brief Exports dedispersed data to disk. */
		void Export_DD_data(int range, float ***dedispersed_data, size_t max_nTimesamples, int const*const ndms, int *inBin, const std::string &base_filename, int *const ranges_to_export, int DMs_per_file) {
			char final_filename[200];

			int inner_DM_shift;
			float *h_data;


			for (int r = 0; r<range; r++) {
				if (ranges_to_export[r]) {
					size_t nTimesamples = max_nTimesamples/inBin[r];
					size_t nDMs = ndms[r];
					if (DMs_per_file<0) DMs_per_file = (int)nDMs;
					size_t export_size = nTimesamples*DMs_per_file;

					int nRepeats = nDMs/DMs_per_file;
					int nRest = nDMs%DMs_per_file;
					std::vector<int> chunk_size;
					for (int f = 0; f<nRepeats; f++) chunk_size.push_back(DMs_per_file);
					if (nRest>0) chunk_size.push_back(nRest);
					printf("  Data in range %d will be exported into %d files\n", r, (int)chunk_size.size());
					float percent_per_DMtrial = 100.0/(nDMs);

					h_data = new float[export_size];

					sprintf(final_filename, "%s_%d.txt", base_filename.c_str(), r);
					std::ofstream FILEOUT;
					FILEOUT.open(final_filename, std::ofstream::out);
					FILEOUT << nTimesamples << std::endl;
					FILEOUT << nDMs << std::endl;
					FILEOUT << DMs_per_file << std::endl;
					FILEOUT << (int)chunk_size.size() << std::endl;
					FILEOUT.close();

					inner_DM_shift = 0;
					for (int i = 0; i<(int)chunk_size.size(); i++) {
						sprintf(final_filename, "%s_%d_%d.dat", base_filename.c_str(), r, i);

						for (int d = 0; d<chunk_size[i]; d++) {
							for (size_t t = 0; t<nTimesamples; t++) {
								h_data[d*nTimesamples + t] = dedispersed_data[r][inner_DM_shift + d][t];
							}
						}

						FILE *fp_out;
						if ((fp_out = fopen(final_filename, "wb")) == NULL) {
							fprintf(stderr, "Error opening output file!\n");
							exit(0);
						}
						fwrite(h_data, nTimesamples*chunk_size[i]*sizeof(float), 1, fp_out);
						fclose(fp_out);

						inner_DM_shift = inner_DM_shift + chunk_size[i];
						printf("    Exported: %0.2f%% \n", ((float)inner_DM_shift)*percent_per_DMtrial);
					}

					delete[] h_data;
				}
			}
		}

		/** \brief Save dedispersed data to disk. (Use as a debug feature only). */
		inline size_t Calculate_sd_per_file_from_file_size(size_t size_in_MB, size_t primary_dimension, size_t floats_per_elements) {
			size_t nFloats = (size_in_MB*1024*1024)/4;
			size_t sd_per_file = (nFloats*floats_per_elements)/primary_dimension;
			return(sd_per_file);
		}
	};

} // namespace astroaccelerate
#endif // ASTRO_ACCELERATE_AA_PIPELINE_RUNNER_HPP

