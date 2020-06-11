#ifndef ASTRO_ACCELERATE_AA_TIMELOG_HPP
#define ASTRO_ACCELERATE_AA_TIMELOG_HPP

#include <map>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "aa_log.hpp"

namespace astroaccelerate {


 /**
   * \class aa_TimeLog aa_timelog.hpp "include/aa_timelog.hpp"
   * \brief Class for logging timing information of concrete components and their subcomponents. Log a message to a file on disk.
   * \author AstroAccelerate team.
   * \date 2 July 2019.
   */
class TimeLog{

        public:
		TimeLog(){};

		typedef std::map<std::pair<std::string, std::string>, double> maptype;
                typedef std::pair<std::map<std::pair<std::string, std::string>, double>::iterator,bool> ret_maptype;

		// adding new component with subcomponent. If it exists then just increase the time.
                void adding(std::string component, std::string sub_component, double time){
                        ret = pattern.insert(make_pair(make_pair(component,sub_component), time));
                        if (ret.second == false){
                                ret.first->second += time;
                        }
                };

		maptype pat(){
			return pattern;
		}

                void print_to_file(){
			time_file.open("time.log", std::ofstream::out);
                        for (const auto &pair : pattern){
				time_file << pair.first.first << " " << pair.first.second << " " << pair.second << '\n';
                                //note pair.first in the row column pair, pair.first.first is the row, pair.first.second is the column, pair.second is the string pattern
                        }
			time_file.close();		
                }

		void print(){
			LOG(log_level::notice,"--------------------------------------------------------------");
			LOG(log_level::notice,"Summary of the times for each component and subcomponent [ms]:")
			for (const auto &pair : pattern){
				LOG(log_level::notice, "\t" + pair.first.first + " " + pair.first.second + " " + "\t" +std::to_string(pair.second));
			}
			LOG(log_level::notice,"--------------------------------------------------------------");
		}
		
		void print(float processed_telescope_time){
			LOG(log_level::notice,"--------------------------------------------------------------");
			LOG(log_level::notice,"Summary of the times for each component and subcomponent [ms]:")
			double max_time = 0;
			for (const auto &pair : pattern){
				LOG(log_level::notice, "\t" + pair.first.first + " " + pair.first.second + " " + "\t" +std::to_string(pair.second));
				if(pair.second > max_time) max_time = pair.second;
			}
			char str[200];
			sprintf(str,"%0.3f", processed_telescope_time/max_time);
			LOG(log_level::notice, std::string("\tProcessed telescope time:\t") + std::to_string(processed_telescope_time*1000.0));
			LOG(log_level::notice, " ");
			LOG(log_level::notice, std::string("\tReal-time speed-up factor:\t") + std::to_string(processed_telescope_time/(max_time/1000.0)));
			LOG(log_level::notice,"--------------------------------------------------------------");
		}

        private:
		std::ofstream time_file;
		static maptype pattern;
                ret_maptype ret;
};

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_TIMELOG_HPP

