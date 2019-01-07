#ifndef ASTRO_ACCELERATE_AA_LOG_HPP
#define ASTRO_ACCELERATE_AA_LOG_HPP

#include <iomanip>
#include <sstream>
#include <chrono>

namespace astroaccelerate {
  
#ifdef ASTRO_ACCELERATE_LOGGING_FACILITY_ENABLE
#if ASTRO_ACCELERATE_LOGGING_FACILITY_ENABLE == 0
  #define LOG_ENABLE 0
#else
  #define LOG_ENABLE 1
#endif
#else
#define LOG_ENABLE 1
#endif
  
  enum class log_level : int {
			      debug = 0,
			      notice,
			      warning,
			      error
  };

#ifndef LOG_THRESHOLD
#define LOG_THRESHOLD log_level::notice
#endif

#define DEBUG_TEXT   "DEBUG   "
#define NOTICE_TEXT  "NOTICE  "
#define WARNING_TEXT "WARNING "
#define ERROR_TEXT   "ERROR   "
  
#define LOG(level, msg) if(!LOG_ENABLE) ;		   \
  else if (level < LOG_THRESHOLD || !FILElog::stream()) ;  \
  else if(level == log_level::debug) FILElog().write(DEBUG_TEXT, msg);	\
  else if(level == log_level::notice) FILElog().write(NOTICE_TEXT, msg); \
  else if(level == log_level::warning) FILElog().write(WARNING_TEXT, msg); \
  else if(level == log_level::error) FILElog().write(ERROR_TEXT, msg);
  
  template <class stream_target> class aa_log {
  public:
    static void write(const char* text, const std::string msg) {
      auto time_point = std::chrono::system_clock::now();
      std::time_t now_c = std::chrono::system_clock::to_time_t(time_point);
      std::stringstream ss;
      ss << std::put_time(std::localtime(&now_c), "%F %H:%M:%S");
      std::string s = ss.str();
      FILE* pStream = stream();
      if(!pStream) {
	printf("Returning\n");
	return;
      }
      fprintf(pStream, "%s - %s %s\n", s.c_str(), text, msg.c_str());
      fflush(pStream);
    }
    
    static FILE*& stream();
  private:
  };

  typedef aa_log<FILE> FILElog;

template <> inline FILE*& aa_log<FILE>::stream() {
    static FILE* file = stderr;
    return file;
}
  

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_LOG_HPP
