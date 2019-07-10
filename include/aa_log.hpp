#ifndef ASTRO_ACCELERATE_AA_LOG_HPP
#define ASTRO_ACCELERATE_AA_LOG_HPP

#include <iomanip>
#include <sstream>
#include <chrono>
#include <time.h>

namespace astroaccelerate {

  /**
   * \brief Preprocessor definition defined by CMake to check whether any CMake overrides should apply.
   * \details If 0 then LOG_ENABLE is set to 0 and all logs are disabled.
   * \details If 1 then LOG_ENABLE is set to 1 and logging is enabled for all logs above the LOG_THRESHOLD (defined further down).
   */  
#ifdef ASTRO_ACCELERATE_LOGGING_FACILITY_ENABLE
#if ASTRO_ACCELERATE_LOGGING_FACILITY_ENABLE
#define LOG_ENABLE 1 /** \brief Enable logging. */
#else
#define LOG_ENABLE 0 /** \brief Disable logging. */
#endif
#else 
#define LOG_ENABLE 1 /** \brief If ASTRO_ACCELERATE_LOGGING_FACILITY_ENABLE was not defined (e.g. because the build system uses make and not cmake), then enable logging by default. */
#endif
  
  /**
   * \enum log_level
   * \brief Contains the log levels.
   **/  
  enum class log_level : int {
    dev_debug = 0, /** \brief dev_debug Developer debugging information. */
      debug, /** \brief User debugging information. */
      notice, /** \brief Informational notices. */
      warning, /** \brief Behaviour that does not terminate execution. */
      error /** \brief Behaviour that terminates execution. */
      };

  /** \brief Defines the LOG_THRESHOLD minimum level for logs to be logged at all. */
#ifndef LOG_THRESHOLD
#define LOG_THRESHOLD log_level::dev_debug
#endif

#define DEV_DEBUG_TEXT "DEBUG   "
#define DEBUG_TEXT     "DEBUG   "
#define NOTICE_TEXT    "NOTICE  "
#define WARNING_TEXT   "WARNING "
#define ERROR_TEXT     "ERROR   "
  
  /** \brief Preprocessor macro that delegates and configures at compile time how parameters of how the log message is logged. */
#define LOG(level, msg) if(LOG_ENABLE == 0) ;				\
  else if (level < LOG_THRESHOLD || !FILElog::stream()) ;		\
  else if(level == log_level::dev_debug) FILElog().write(DEV_DEBUG_TEXT, msg); \
  else if(level == log_level::debug) FILElog().write(DEBUG_TEXT, msg);	\
  else if(level == log_level::notice) FILElog().write(NOTICE_TEXT, msg); \
  else if(level == log_level::warning) FILElog().write(WARNING_TEXT, msg); \
  else if(level == log_level::error) FILElog().write(ERROR_TEXT, msg);
  
  /**
   * \class aa_log aa_log.hpp "include/aa_log.hpp"
   * \brief Class for logging library information. Log a message to the console or to a file on disk.
   * \details Detailed instructions are included in MANUAL.md in the repository.
   * \author AstroAccelerate team.
   * \date 3 July 2019.
   */
  template <class stream_target> class aa_log {
  public:
    /**
     * \brief Method to write a message to the output stream for the log_level supplied.
     * \details The log level is denoted by text, the message is denoted by msg.
     */
	static void write(const char* text, const std::string msg) {
		std::stringstream ss;
		time_t rawtime;
		struct tm * timeinfo;
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		ss << asctime (timeinfo);
		std::string timedate = ss.str();
		
		std::string s(timedate, 4, timedate.size()-5);
		FILE* pStream = stream();
		if(!pStream) {
			printf("Returning\n");
			return;
		}
		fprintf(pStream, "%s - %s %s\n", s.c_str(), text, msg.c_str());
		fflush(pStream);
    }
    
    /** \brief Set the stream to a file pointer (default standard output to console). */
    static FILE*& stream();
  private:
  };

  typedef aa_log<FILE> FILElog;

  /** \brief Templated specialisation for setting the stream to a FILE* stream. */
  template <> inline FILE*& aa_log<FILE>::stream() {
    static FILE* file = stderr;
    return file;
  }
} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_LOG_HPP
