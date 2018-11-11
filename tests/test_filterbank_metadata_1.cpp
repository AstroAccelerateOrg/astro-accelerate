#include <iostream>
#include "aa_filterbank_metadata.hpp"

using namespace astroaccelerate;

int main() {

  const double az_start = 0;
  const double za_start = 0;
  const double src_raj = 0;
  const double src_dej = 0;
  const double tstart = 0;
  const double tsamp = 0;
  const double refdm = 0;
  const double period = 0;
  const double fch1 = 0;
  const double foff = 0;
  const double fchannel = 0;

  const int telescope_id = 0;
  const int machine_id = 0;
  const int data_type = 0;
  const int barycentric = 0;
  const int pulsarcentric = 0;
  const int nbits = 8;
  const int nsamples = 1048576;
  const int nchans = 4096;
  const int nifs = 0;

  const char FREQUENCY_START = '0';
  const char FREQUENCY_END = '0';

  const std::string rawdatafile = "unknown";
  const std::string source_name = "unknown";

  const aa_filterbank_metadata meta(telescope_id,
                                    machine_id,
                                    data_type,
                                    rawdatafile,
                                    source_name,
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

  bool pass = true;
  if(meta.telescope_id() != telescope_id) {
    std::cout << "Fail test: telescope_id" << std::endl;
    pass = false;
  }

  if(meta.machine_id() != machine_id) {
    std::cout << "Fail test: machine_id" << std::endl;
    pass = false;
  }

  if(meta.data_type() != data_type) {
    std::cout << "Fail test: data_type" << std::endl;
    pass = false;
  }

  if(meta.rawdatafile() != rawdatafile) {
    std::cout << "Fail test: rawdatafile" << std::endl;
    pass = false;
  }

  if(meta.source_name() != source_name) {
    std::cout << "Fail test: source_name" << std::endl;
    pass = false;
  }

  if(meta.barycentric() != barycentric) {
    std::cout << "Fail test: barycentric" << std::endl;
    pass = false;
  }

  if(meta.pulsarcentric() != pulsarcentric) {
    std::cout << "Fail test: pulsarcentric" << std::endl;
    pass = false;
  }

  if(meta.az_start() != az_start) {
    std::cout << "Fail test: az_start" << std::endl;
    pass = false;
  }

  if(meta.za_start() != za_start) {
    std::cout << "Fail test: za_start" << std::endl;
    pass = false;
  }

  if(meta.src_raj() != src_raj) {
    std::cout << "Fail test: src_raj" << std::endl;
    pass = false;
  }

  if(meta.src_dej() != src_dej) {
    std::cout << "Fail test: src_dej" << std::endl;
    pass = false;
  }

  if(meta.tstart() != tstart) {
    std::cout << "Fail test: tstart" << std::endl;
    pass = false;
  }

  if(meta.tsamp() != tsamp) {
    std::cout << "Fail test: tsamp" << std::endl;
    pass = false;
  }

  if(meta.nbits() != nbits) {
    std::cout << "Fail test: nbits" << std::endl;
    pass = false;
  }

  if(meta.nsamples() != nsamples) {
    std::cout << "Fail test: nsamples" << std::endl;
    pass = false;
  }

  if(meta.fch1() != fch1) {
    std::cout << "Fail test: fch1" << std::endl;
    pass = false;
  }

  if(meta.foff() != foff) {
    std::cout << "Fail test: foff" << std::endl;
    pass = false;
  }

  if(meta.FREQUENCY_START() != FREQUENCY_START) {
    std::cout << "Fail test: FREQUENCY_START" << std::endl;
    pass = false;
  }

  if(meta.fchannel() != fchannel) {
    std::cout << "Fail test: fchannel" << std::endl;
    pass = false;
  }

  if(meta.FREQUENCY_END() != FREQUENCY_END) {
    std::cout << "Fail test: FREQUENCY_END" << std::endl;
    pass = false;
  }

  if(meta.nchans() != nchans) {
    std::cout << "Fail test: nchans" << std::endl;
    pass = false;
  }
  
  if(meta.nifs() != nifs) {
    std::cout << "Fail test: nifs" << std::endl;
    pass = false;
  }

  if(meta.refdm() != refdm) {
    std::cout << "Fail test: refdm" << std::endl;
    pass = false;
  }

  if(meta.period() != period) {
    std::cout << "Fail test: period" << std::endl;
    pass = false;
  }

  //Test result
  if(pass) {
    std::cout << "test_filterbank_metadata_1 pass" << std::endl;
  }
  else {
    std::cout << "test_filterbank_metadata_1 fail" << std::endl;
  }

  
  return 0;
}
