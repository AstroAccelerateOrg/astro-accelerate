#include <iostream>
#include <sstream>

#include "aa_device_info.hpp"

using namespace astroaccelerate;

int main(int argc, char* argv[]) {
  //aa_device_info& device_info = aa_device_info::instance();

  if(!argc) {
    std::cout << "Fail." << std::endl;
    return 0;
  }
  
  
  std::stringstream s1;
  s1 << argv[1];
  /*
  
  size_t tmp = 0;
  s1 >> tmp;
  const size_t test_request = tmp;
  device_info.request(test_request);
  std::cout << device_info.requested() << std::endl;
  if(device_info.requested() != test_request) {
    std::cout << "Fail." << std::endl;
    return 0;
  }
  */

  std::cout << "Runs." << std::endl;
  
  return 0;
}
