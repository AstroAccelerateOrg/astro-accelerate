#ifndef ASTRO_ACCELERATE_AA_FILTERBANK_METADATA_HPP
#define ASTRO_ACCELERATE_AA_FILTERBANK_METADATA_HPP

#include <string>

namespace astroaccelerate {

/**
 * Data description from "SIGPROC-v3.7 (Pulsar) Signal Processing Programs"
 * Source: http://sigproc.sourceforge.net/sigproc.pdf
 */

class aa_filterbank_metadata {
public:
  aa_filterbank_metadata() {
        
  }
    
  aa_filterbank_metadata(const aa_filterbank_metadata &meta) : m_az_start(meta.m_az_start),
							       m_za_start(meta.m_za_start),
							       m_src_raj(meta.m_src_raj),
							       m_src_dej(meta.m_src_dej),
							       m_tstart(meta.m_tstart),
							       m_tsamp(meta.m_tsamp),
							       m_refdm(meta.m_refdm),
							       m_period(meta.m_period),
							       m_fch1(meta.m_fch1),
							       m_foff(meta.m_foff),
							       m_fchannel(meta.m_fchannel),
							       m_telescope_id(meta.m_telescope_id),
							       m_machine_id(meta.m_machine_id),
							       m_data_type(meta.m_data_type),
							       m_barycentric(meta.m_barycentric),
							       m_pulsarcentric(meta.m_pulsarcentric),
							       m_nbits(meta.m_nbits),
							       m_nsamples(meta.m_nsamples),
							       m_nchans(meta.m_nchans),
							       m_nifs(meta.m_nifs),
							       m_FREQUENCY_START(meta.m_FREQUENCY_START),
							       m_FREQUENCY_END(meta.m_FREQUENCY_END),
							       m_rawdatafile(meta.m_rawdatafile),
							       m_source_name(meta.m_source_name) {
    
  }
  
    aa_filterbank_metadata(const int &telescope_id,
                           const int &machine_id,
                           const int &data_type,
                           const std::string &rawdatafile,
                           const std::string &source_name,
                           const int &barycentric,
                           const int &pulsarcentric,
                           const double &az_start,
                           const double &za_start,
                           const double &src_raj,
                           const double &src_dej,
                           const double &tstart,
                           const double &tsamp,
                           const int &nbits,
                           const int &nsamples,
                           const double &fch1,
                           const double &foff,
                           const char &FREQUENCY_START,
                           const double &fchannel,
                           const char &FREQUENCY_END,
                           const int &nchans,
                           const int &nifs,
                           const double &refdm,
                           const double &period
                           )
    : m_az_start(az_start),
      m_za_start(za_start),
      m_src_raj(src_raj),
      m_src_dej(src_dej),
      m_tstart(tstart),
      m_tsamp(tsamp),
      m_refdm(refdm),
      m_period(period),
      m_fch1(fch1),
      m_foff(foff),
      m_fchannel(fchannel),
      m_telescope_id(telescope_id),
      m_machine_id(machine_id),
      m_data_type(data_type),
      m_barycentric(barycentric),
      m_pulsarcentric(pulsarcentric),
      m_nbits(nbits),
      m_nsamples(nsamples),
      m_nchans(nchans),
      m_nifs(nifs),
      m_FREQUENCY_START(FREQUENCY_START),
      m_FREQUENCY_END(FREQUENCY_END),
      m_rawdatafile(rawdatafile),
      m_source_name(source_name) {
        
    }

  aa_filterbank_metadata(const double &tstart,
			 const double &tsamp,
			 const int &nbits,
			 const int &nsamples,
			 const double &fch1,
			 const double &foff,
			 const int &nchans,
			 const int &nifs) : m_az_start(0),
			     m_za_start(0),
			     m_src_raj(0),
			     m_src_dej(0),
			     m_tstart(tstart),
			     m_tsamp(tsamp),
			     m_refdm(0),
			     m_period(0),
			     m_fch1(fch1),
			     m_foff(foff),
			     m_fchannel(0),
			     m_telescope_id(0),
			     m_machine_id(0),
			     m_data_type(0),
			     m_barycentric(0),
			     m_pulsarcentric(0),
			     m_nbits(nbits),
			     m_nsamples(nsamples),
			     m_nchans(nchans),
			     m_nifs(nifs),
			     m_FREQUENCY_START(0),
			     m_FREQUENCY_END(0),
			     m_rawdatafile(""),
			     m_source_name("") {

  }
    
    ~aa_filterbank_metadata() {
        
    }
    
    double az_start() const {
        return m_az_start;
    }
    
    double za_start() const {
        return m_za_start;
    }
    
    double src_raj() const {
        return m_src_raj;
    }
    
    double src_dej() const {
        return m_src_dej;
    }
    
    double tstart() const {
        return m_tstart;
    }
    
    double tsamp() const {
        return m_tsamp;
    }
    
    double refdm() const {
        return m_refdm;
    }
    
    double period() const {
        return m_period;
    }
    
    double fch1() const {
        return m_fch1;
    }
    
    double foff() const {
        return m_foff;
    }
    
    double fchannel() const {
        return m_fchannel;
    }
    
    int telescope_id() const {
        return m_telescope_id;
    }
    
    int machine_id() const {
        return m_machine_id;
    }
    
    int data_type() const {
        return m_data_type;
    }
    
    int barycentric() const {
        return m_barycentric;
    }
    
    int pulsarcentric() const {
        return m_pulsarcentric;
    }
    
    int nbits() const {
        return m_nbits;
    }
    
    int nsamples() const {
        return m_nsamples;
    }
    
    int nchans() const {
        return m_nchans;
    }
    
    int nifs() const {
        return m_nifs;
    }
    
    char FREQUENCY_START() const {
        return m_FREQUENCY_START;
    }
    
    char FREQUENCY_END() const {
        return m_FREQUENCY_END;
    }
    
    const std::string rawdatafile() const {
        return m_rawdatafile;
    }
    
    const std::string source_name() const {
        return m_source_name;
    }
    
    int N() const {
        return (m_nifs * m_nchans * m_nsamples);
    }
    
    int array_index(const int &sample_idx, const int &channel_idx, const int &frequency_channel) const {
        return ((sample_idx * m_nifs * m_nchans)
                + (channel_idx * m_nchans + frequency_channel));
    }
    
    double sky_frequency(const int &channel_idx) const {
        return (m_fch1 + channel_idx * m_foff);
    }
    
private:
    double m_az_start;
    double m_za_start;
    double m_src_raj;
    double m_src_dej;
    double m_tstart;
    double m_tsamp;
    double m_refdm;
    double m_period;
    double m_fch1;
    double m_foff;
    double m_fchannel;
    
    int m_telescope_id;
    int m_machine_id;
    int m_data_type;
    int m_barycentric;
    int m_pulsarcentric;
    int m_nbits;
    int m_nsamples;
    int m_nchans;
    int m_nifs;
    
    char m_FREQUENCY_START;
    char m_FREQUENCY_END;
    
    std::string m_rawdatafile;
    std::string m_source_name;
};

} //namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_FILTERBANK_METADATA_HPP
