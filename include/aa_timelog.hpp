#ifndef ASTRO_ACCELERATE_AA_TIMELOG_HPP
#define ASTRO_ACCELERATE_AA_TIMELOG_HPP

#include <map>
#include <algorithm>
#include <iostream>

namespace astroaccelerate {


 /**
   * \class aa_TimeLog aa_timelog.hpp "include/aa_timelog.hpp"
   * \brief Class for logging timing information of concrete components and their subcomponents. Log a message to a file on disk.
   * \author AstroAccelerate team.
   * \date 2 July 2019.
   */
class TimeLog{

        public:

		// adding new component with subcomponent. If it exists then just increase the time.
                void adding(std::string component, std::string sub_component, double time){
                        ret = pattern.insert(make_pair(make_pair(component,sub_component), time));
                        if (ret.second == false){
                                ret.first->second += time;
                        }
                };

                void print(){
                        for (const auto &pair : pattern) {
                                std::cout << pair.first.first << " " << pair.first.second << " " << pair.second << '\n';
                                //note pair.first in the row column pair, pair.first.first is the row, pair.first.second is the column, pair.second is the string pattern
                        }
                }

                typedef std::map<std::pair<std::string, std::string>, double> maptype;
                typedef std::pair<std::map<std::pair<std::string, std::string>, double>::iterator,bool> ret_maptype;

        private:
                maptype pattern;
                ret_maptype ret;
};

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_TIMELOG_HPP

