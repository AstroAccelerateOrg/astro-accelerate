#include "aa_fdas_util.hpp"

namespace astroaccelerate {

  /** \brief Print usage instructions for fdas. */
  void print_usage()
  {

    printf("\n        ----------------------------------        \n\n");
    printf("The program reads an accelerated sinusoid or optionally a pulse with a number of harmonics,\nand performs a (single harmonic) fourier domain acceleration search on this signal on a GPU\n");
    printf("\nUsage:  fdas [OPTION] \n\n");
    printf("Note:These command line options do not include all useful parameters.\nTo change correlation and  hardware related parameters you can edit the file params.h\n");
    printf("\n OPTIONS:\n\n");
    printf(" %-13s %-2s %-.55s\n","-help", ": ","Print this message\n");
    printf(" %-13s %-2s %-.55s\n","-fname [name]", ": ","Signal file name\n");
    printf(" %-13s %-2s %-.55s\n","-nname [name]", ": ","Noise data file name (with -wsig)\n");
    printf(" %-13s %-2s %-.55s\n %-16s %-.55s\n %-16s %-.55s\n","-nharm [n]",":", "Number of harmonics to add in signal (default=1).","","Only used with -wsig, it inserts the specified number", "", "of harmonics to the signal created\n");
    printf(" %-13s %-2s %-.55s\n %-16s %-.55s\n%-16s %-.55s\n","-nsig [n]",":", "Noise multiplier (integer, default is 0, i.e. no noise)."," ", "Currently only use when a noise input file is available","", " under the 'data' directory.\n");
    printf(" %-13s %-2s %-.55s\n", "-duty [f]", ":","Duty cucle of pulse (%) (default is  10).\n"); 
    printf(" %-13s %-2s %-.60s\n %-16s %-.60s\n", "-iter [n]",":", "Number of iterations to run the acceleration search."," ", "(default is 1).\n");
    printf(" %-13s %-2s %-.60s\n", "-writef [0/1]", ":", "Write results to file for plots (default is no (0)).\n");
    printf(" %-13s %-2s %-.60s\n %-16s %-.55s\n", "-zval [n]", ":", "Value of z in bins drifted when creating a signal file", "", "(default is 4).\n");

    printf(" %-13s %-2s %-.60s\n %-16s %-.55s\n", "-mul [n]", ":", "Multiplier for non-power of two signal lengths,", "", "multiplies by 8192 to get the signal length (default is 1024).\n");
    printf(" %-13s %-2s %-.60s\n %-16s %-.55s\n", "-wsig", ":", "Create and write a time series to binary file", "", "- use to create signal files (default is no).\n");
    printf(" %-13s %-2s %-.60s\n", "-search", ":", "Perform the full search loop (default is no).\n");
    printf(" %-13s %-2s %-.60s\n", "-basic", ":", "Use with -search to search using the basic method.\n");

    printf(" %-13s %-2s %-.60s\n", "-kfft", ":", "Use with -search to search using custom FFT.\n");
    printf(" %-13s %-2s %-.60s\n", "-thresh [p]", ":", "Normalized power threshold, real number (default is 10.0).\n");
    printf(" %-13s %-2s %-.60s\n %-16s %-.55s\n", "-freq0 [f]", ":", "Central frequency of signal, real number.", "", "Used with -wsig. (default is 100.5).\n");
    printf(" %-13s %-2s %-.60s\n %-16s %-.55s\n", "-sigamp [f]", ":", "Amplitude of signal to create (only with -wsig)","", "(default is 0.1).\n");
    printf(" %-13s %-2s %-.60s\n", "-devid [i]", ":", "CUDA device id to use (default is 0).\n");
    printf(" %-13s %-2s %-.60s\n", "-inbin", ":", "Perform interbinning on the complex output.\n"); 
    printf(" %-13s %-2s %-.60s\n", "-norm", ":", "PRESTO block median normalization.\n"); 

    printf("\n        ----------------------------------         \n\n");
  }

  /** \brief Read command line arguments that can be used to configure fdas. */
  void read_command_line(int argc, char *argv[], cmd_args *args)
  {
    // read from command line
    static struct option long_options[] = {
      /* name         has_arg      flag   val */
      /* ------------------------------------------ */
      {"help", no_argument,   NULL,   0 },
      {"nharm", required_argument,   NULL,   1 },
      {"nsig", required_argument,   NULL,   2 },
      {"duty", required_argument,   NULL,   3 },
      {"iter", required_argument,   NULL,   4 },
      {"writef", required_argument, NULL,   5 },
      {"zval", required_argument, NULL,   6 },
      {"wsig", no_argument, NULL,   7 },
      {"search", no_argument, NULL,   8 },
      {"thresh", required_argument, NULL,   9 },
      {"freq0", required_argument, NULL,   10 },
      {"sigamp", required_argument, NULL,   11 },
      {"fname", required_argument, NULL,   12 },
      {"nname", required_argument, NULL,   13 },
      {"devid", required_argument, NULL,   14 },
      {"basic", no_argument, NULL,   15 },
      {"kfft", no_argument, NULL,   16 },
      {"mul", required_argument, NULL,   17 },
      {"inbin", no_argument, NULL,   18 },
      {"norm", no_argument, NULL,   19 },
      {   0,            0,           0,    0  }
    };

    int long_index =0;
    int opt;

    while ((opt = getopt_long_only(argc, argv,"", long_options, &long_index )) != -1) {
      switch (opt) {

      case 0:
	print_usage();
	exit(EXIT_FAILURE);
	break;

      case 1 : 
	args->nharms = atoi(optarg); 
	printf("\n nharm = %d\n", args->nharms);
	break;

      case 2 :
	args->nsig = atoi(optarg); 
	printf("\n nsig = %d\n", args->nsig);
	break;

      case 3 : 
	args->duty = 0.01*atof(optarg); 
	printf("\n duty cycle requested: %f\n", args->duty); 
	break;

      case 4 : 
	args->iter = atoi(optarg); 
	printf("\n Number of iterations to run: %d\n", args->iter); 
	break;

      case 5 : 
	args->writef = atoi(optarg); 
	if (args->writef) 
	  printf("\n Will write result data output to file for plots \n"); 
	break;

      case 6 : 
	args->zval = atoi(optarg); 
	printf("\n Acceleration value in bins (z) is set to %d\n", args->zval); 
	break;

      case 7 : 
	args->wsig = 1; 
	printf("\n Will write time series to binary file\n"); 
	break;

      case 8 : 
	args->search = 1; 
	printf("\n Will perform full search\n"); 
	break;

      case 9 : 
	args->thresh = atof(optarg); 
	printf("\n Power threshold =  %f\n", args->thresh); 
	break;

      case 10 : 
	args->freq0 = atof(optarg); 
	printf("\n Signal central frequency =  %f\n", args->freq0); 
	break;

      case 11 : 
	args->sigamp = atof(optarg); 
	printf("\n Signal amplitude: %.2f \n", args->sigamp); 
	break;


      case 12 : 
	args->afname=optarg; 
	sprintf(args->afname,"%s",optarg);
	printf("\n accelerated signal filename : %s \n",args->afname); 
	break;

      case 13 : 
	args->nfname=optarg; 
	sprintf(args->nfname,"%s",optarg);
	printf("\n noise filename : %s \n",args->nfname); 
	break;

      case 14 : 
	args->devid = atoi(optarg);
	printf("\nDevice requested: %d \n", args->devid); 
	break;

      case 15 : 
	args->basic = 1;
	printf("\n Will run basic search \n"); 
	break;

      case 16 : 
	args->kfft = 1;
	printf("\n Will run search using custom FFT \n"); 
	break;

      case 17 : 
	args->mul = atoi(optarg); 
	printf("\n Number of points in time series is %d =  8192*%d\n", args->mul*8192, args->mul); 
	break;

      case 18 : 
	args->inbin=1; 
	printf("\n Will perform Interbinning on the complex output \n"); 
	break;

      case 19 : 
	args->norm=1; 
	printf("\n Will perform normalization \n"); 
	break;

      default:
	print_usage(); 
	exit(EXIT_FAILURE);
      }
    }
    //  printf("\noptind = %d, argc = %d\n",optind,argc);
    if (optind < argc || optind == 1) {
      print_usage(); 
      exit(EXIT_FAILURE);
    }

  }

  /** \brief Read an input file for fdas processing. */
  int read_input_file(char *fname, float **array)
  {
    FILE *data;
    struct stat fs;
    int fd, ret;
    unsigned long len;

    printf("\nReading time series from file %s\n",fname);
    //  data = fopen(fname,"rb"); //first try open for reading
    if ((data = fopen(fname,"rb")) == NULL){
      printf("\nError opening file %s for reading: %s\nRun with -wsig to create file and try again\n", fname, strerror(errno));
      exit(1);
    }

    fd = fileno(data);
    ret = fstat(fd,&fs);
    if(ret == -1){
      perror("\nfstat() error");
      exit(-1);
    }
    len = fs.st_size/sizeof(float);
    printf("\nfilelen in floats: %zu\n",len);
    *array = (float*)malloc(len*sizeof(float));
    if(*array==NULL){
      printf("\narray allocation failed! %s\n",strerror(errno));
      exit(1);
    }

    if (data!=NULL){ // file exists
      if (fread (*array, sizeof(float), len, data)!= len){
	printf("\nError reading file %s : %s\n", fname, strerror(errno));
	exit(1);
      }
    }
    printf("\nFile read successful\n");

    fclose(data);
    return len;
  }

  /** \brief Read fourier series data from disk. */
  void read_input_file_cpx(char *fname, float2 *array)
  {
    FILE *data;
    struct stat fs;
    int fd, ret;
    unsigned long len;

    printf("\nReading complex fourier series from file %s\n",fname);
    data = fopen(fname,"rb"); //first try open for reading
    if (data!=NULL){ // file exists
      fd = fileno(data);
      ret = fstat(fd,&fs);
      if(ret == -1){
	perror("\nfstat() error");
	exit(-1);
      }
      len = fs.st_size/sizeof(float2);
      if (fread (array, sizeof(float2), len, data)!= len){
	printf("\nError reading file %s : %s\n", fname, strerror(errno));
	exit(1);
      }
    }
    else{
      printf("\nError opening file %s for reading: %s\nRun with -wsig to create file and try again\n", fname, strerror(errno));
      fclose(data);
      exit(1);
    }
    printf("\nFile read successful\n");
    fclose(data);
  }

  /** \brief Write fdas output to disk. */
  void write_output_file(char *fname, float **array,  unsigned long len)
  {
    FILE *data;

    printf("\nWriting data to file %s\n",fname) ;

    data = fopen(fname, "wb");
    if (data == NULL){
      printf("\nError opening file %s for writing : %s\n", fname, strerror(errno));
      exit(1);
    }
    if (fwrite(*array, sizeof(float), len, data) != len){
      printf("\n Error writing to file %s : %s\n",fname, strerror(errno));
      exit(1);
    }
    printf("\nFile write successful\n");

    fclose(data);

  }
} //namespace astroaccelerate
