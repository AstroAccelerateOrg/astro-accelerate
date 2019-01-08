#include <iostream>
#include "aa_device_info.hpp"
#include "aa_ddtr_strategy.hpp"

using namespace astroaccelerate;

int main() {
  std::cout << "Running test_ddtr_strategy_1.cpp" << std::endl;
  aa_device_info* device_info = aa_device_info::instance();
  aa_device_info::CARD_ID selected_card = 0;
  aa_device_info::aa_card_info selected_card_info;
  if(device_info->init_card(selected_card, selected_card_info)) {
    std::cout << "NOTICE: init_card complete." << std::endl;
  }
  else {
    std::cout << "ERROR: init_card incomplete." << std::endl;
    return 0;
  }
  
  aa_ddtr_plan ddtr_plan;
  ddtr_plan.add_dm(0, 370, 0.307, 1, 1);

  double az_start = 0;
  double za_start = 0;
  double src_raj = 0;
  double src_dej = 0;
  double tstart = 0;
  double tsamp = 0;
  double refdm = 0;
  double period = 0;
  double fch1 = 0;
  double foff = 0;
  double fchannel = 0;

  int telescope_id = 0;
  int machine_id = 0;
  int data_type = 0;
  int barycentric = 0;
  int pulsarcentric = 0;
  int nbits = 8;
  int nsamples = 1048576;
  int nchans = 4096;
  int nifs = 0;

  char FREQUENCY_START = '0';
  char FREQUENCY_END = '0';

  std::string rawdatafile = "unknown";
  std::string source_name = "unknown";

  
  const aa_filterbank_metadata meta(telescope_id,
				    machine_id,
				    data_type,
				    "unknown",
				    "unknown",
				    barycentric,
				    pulsarcentric,
				    az_start,
				    za_start,
				    src_raj,
				    src_dej,
				    tstart,
				    tsamp,
				    nbits,
				    nsamples,
				    fch1,
				    foff,
				    FREQUENCY_START,
				    fchannel,
				    FREQUENCY_END,
				    nchans,
				    nifs,
				    refdm,
				    period);

  const bool enable_analysis = false;
  aa_ddtr_strategy ddtr_test_1(ddtr_plan, meta, selected_card_info.free_memory, enable_analysis);
  
  std::cout << "Runs" << std::endl;
  return 0;
}
