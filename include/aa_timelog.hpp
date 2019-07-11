#ifndef ASTRO_ACCELERATE_AA_TIMELOG_HPP
#define ASTRO_ACCELERATE_AA_TIMELOG_HPP

#include <map>
#include <algorithm>
#include <iostream>
#include <fstream>

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
//		TimeLog(){
//			time_file.open("time.log", std::ofstream::trunc);
//			time_file << "#Component Sub-Component Time \n";
//			time_file.close();
//		}

//		~TimeLog(){
//			print();
//		}

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

                void print(){
			time_file.open("time.log", std::ofstream::out);
                        for (const auto &pair : pattern) {
				time_file << pair.first.first << " " << pair.first.second << " " << pair.second << '\n';
                                //note pair.first in the row column pair, pair.first.first is the row, pair.first.second is the column, pair.second is the string pattern
                        }
			time_file.close();		
                }

        private:
		std::ofstream time_file;
		static maptype pattern;
                ret_maptype ret;
};

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_TIMELOG_HPP
