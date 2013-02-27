#include <stdio.h>
#include <stdlib.h>

void init(int argc, char *argv[], FILE **fp, int *tdms, float *dm_step, float *dm_low) {

	// Local variables
	
	int i;

	*tdms = 0;
	*dm_step = *dm_low = 0.0f;	

	//{{{ Read in the command line parameters and open the input file

	if (argc < 2){
		fprintf(stderr, "Need data file.\n");
		exit(0);
	} else if (argc == 2) {
		if ((*fp=fopen(argv[1], "rb")) == NULL) {
			fprintf(stderr, "Invalid data file!\n");
			exit(0);
		}
	} else {
		for(i = 0; i < argc; i++) {
			if (strcmp(argv[i],"-tdms") == 0) *tdms = atoi(argv[i+1]);
			if (strcmp(argv[i],"-dm_step") == 0) *dm_step = atof(argv[i+1]);
			if (strcmp(argv[i],"-dm_low") == 0) *dm_low = atof(argv[i+1]);
		//	if (strcmp(argv[i],"-L1") == 0) kernel_type = 1;
		}
		if ((*fp=fopen(argv[i-1], "rb")) == NULL) {
			fprintf(stderr, "Invalid data file!\n");
			exit(0);
		}
	}

	if(*tdms == 0) *tdms = 2000;
	if(*dm_step == 0.0f) *dm_step = 0.01;

	

	printf("\ntdms:\t\t%d", *tdms);
	printf("\ndm_step:\t%f", *dm_step);
	printf("\ndm_low:\t\t%f", *dm_low);

	//}}}


}
