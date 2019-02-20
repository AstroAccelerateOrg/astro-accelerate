#ifndef ASTRO_ACCELERATE_AA_FAKE_SIGNAL_GENERATOR_HPP
#define ASTRO_ACCELERATE_AA_FAKE_SIGNAL_GENERATOR_HPP

#include <vector>
#include "aa_ddtr_strategy.hpp"
#include "aa_filterbank_metadata.hpp"

namespace astroaccelerate {

	class aa_fake_signal_metadata {
		public:
			aa_fake_signal_metadata(const double &dm_pos, const int &signal_start, const int &signal_width) : 
				m_dm_position(dm_pos),
				m_signal_start(signal_start),
				m_signal_width(signal_width)
				{}

			double get_dm_pos() const{
				return m_dm_position;
			}

			int get_signal_start() const{
				return m_signal_start;
			}

			int get_signal_width() const{
				return m_signal_width;
			}


		private:
			double m_dm_position;
			int m_signal_start;
			int m_signal_width;
	};


	class aa_fake_signal_generator {

		public: 
	
			void produce_mask(float sigma, int width);
			void get_fake_signal(const aa_ddtr_strategy &strategy, const aa_fake_signal_metadata &m_fake, const aa_filterbank_metadata &m_signal);
			/** \returns inverse Gaussian shape for scaling the fake data */
			std::vector<float> get_mask_data() const{
				return mask_data;
			}
			int max_pos(){
				return maximum_pos;
			}
			std::vector<unsigned short> signal_data() const{
				return signal_output;
			}

		private:
			std::vector<float> mask_data;
			std::vector<unsigned short> signal_output;
			int maximum_pos = 0;
	};



} // namespace astroaccelerate
#endif // ASTRO_ACCELERATE_AA_FAKE_SIGNAL_GENERATOR_HPP

