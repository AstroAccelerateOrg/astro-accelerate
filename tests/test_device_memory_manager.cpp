#include <iostream>
#include <sstream>

#include "aa_device_memory_manager.hpp"

using namespace astroaccelerate;

int main(int argc, char* argv[]) {
  aa_device_memory_manager mem;

  if(!argc) {
    std::cout << "Fail." << std::endl;
    return 0;
  }
  
  std::stringstream s1;
  s1 << argv[1];
  
  size_t tmp = 0;
  s1 >> tmp;
  const size_t test_request = tmp;
  mem.request(test_request);
  std::cout << mem.requested() << std::endl;

  if(mem.requested() != test_request) {
    std::cout << "Fail." << std::endl;
    return 0;
  }


  std::cout << "Runs." << std::endl;
  
  return 0;
}
