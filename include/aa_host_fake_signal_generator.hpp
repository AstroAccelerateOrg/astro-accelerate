#ifndef ASTRO_ACCELERATE_AA_FAKE_SIGNAL_GENERATOR_HPP
#define ASTRO_ACCELERATE_AA_FAKE_SIGNAL_GENERATOR_HPP

#include <vector>
#include "aa_ddtr_strategy.hpp"
#include "aa_filterbank_metadata.hpp"

#include "aa_log.hpp"

namespace astroaccelerate {

	class aa_fake_signal_metadata {
		public:
			aa_fake_signal_metadata(const double &dm_pos, const int &signal_start, const int &signal_width, const float &sigma) : 
				m_dm_position(dm_pos),
				m_signal_start(signal_start),
				m_signal_width(signal_width),
				m_signal_sigma(sigma)
				{}

			aa_fake_signal_metadata(const double &dm_pos, const int &signal_width, const float &sigma, const int &periodicity) :
                                m_dm_position(dm_pos),
                                m_signal_width(signal_width),
                                m_signal_sigma(sigma),
				m_signal_period(periodicity)
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
		
			float get_signal_sigma() const{
				return m_signal_sigma;
			}

			int get_signal_period() const{
				return m_signal_period;
			}

		private:
			double m_dm_position;
			int m_signal_start = 0;
			int m_signal_width;
			float m_signal_sigma;
			int m_signal_period = 0;
	};

	class aa_fake_signal_generator {

		public: 
			bool produce_mask(const float sigma, const int width);

			bool create_signal(const aa_ddtr_strategy &strategy, const aa_fake_signal_metadata &m_fake, const aa_filterbank_metadata &m_signal, bool &dump_to_disk){
				produce_mask(m_fake.get_signal_sigma(), m_fake.get_signal_width());
				get_fake_signal(strategy, m_fake, m_signal, dump_to_disk);
				return 0;
			}

			bool get_fake_signal(const aa_ddtr_strategy &strategy, const aa_fake_signal_metadata &m_fake, const aa_filterbank_metadata &m_signal, bool &dump_to_disk);

			/** \returns inverse Gaussian shape for scaling the fake data */
			std::vector<float> get_mask_data() const{
				return mask_data;
			}
		
			int max_pos() const{
				return maximum_pos;
			}

			std::vector<unsigned short> signal_data() const{
				return signal_output;
			}

			bool ready() const {
				return m_ready;
			}

			/** \brief Static member function that prints member data for an aa_host_fake_signal object. */
			static bool print_info(const aa_fake_signal_metadata &m_fake) {
				LOG(log_level::dev_debug, "------------------------------------------------------------->");
				LOG(log_level::dev_debug, "FAKE SIGNAL INFORMATION:");
				LOG(log_level::dev_debug, "DM position of signal:\t\t" + std::to_string(m_fake.get_dm_pos()));
				LOG(log_level::dev_debug, "Signal start at [samples]:\t" + std::to_string(m_fake.get_signal_start()));
				LOG(log_level::dev_debug, "Signal width [samples]:\t\t" + std::to_string(m_fake.get_signal_width()));
				LOG(log_level::dev_debug, "Signal sigma selected:\t\t" + std::to_string(m_fake.get_signal_sigma()));
				LOG(log_level::dev_debug, "Periodicity selected [ms]:\t" + std::to_string(m_fake.get_signal_period()))
				LOG(log_level::dev_debug, "-------------------------------------------------------------<");
			return true;
			}

		private:
			std::vector<float> mask_data;
			std::vector<unsigned short> signal_output;
			int maximum_pos = 0;
			bool m_ready = false;
			bool m_ready_mask = false;
	};

} // namespace astroaccelerate
#endif // ASTRO_ACCELERATE_AA_FAKE_SIGNAL_GENERATOR_HPP

