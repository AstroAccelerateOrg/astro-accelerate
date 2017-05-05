
static const int UNROLLS_ARRAY[] =  {UNROLLS_0,UNROLLS_1,UNROLLS_2,UNROLLS_3};
static const int SNUMREG_ARRAY[] =  {SNUMREG_0,SNUMREG_1,SNUMREG_2,SNUMREG_3};
static const int SDIVINT_ARRAY[] =  {SDIVINT_0,SDIVINT_1,SDIVINT_2,SDIVINT_3};
static const int SDIVINDM_ARRAY[] = {SDIVINDM_0,SDIVINDM_1,SDIVINDM_2,SDIVINDM_3};

typedef void (*shared_dedisperse_kernel_FTYPE)(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep);

shared_dedisperse_kernel_FTYPE shared_dedisperse_kernel_FARRAY[] = {
	shared_dedisperse_kernel_range_0,shared_dedisperse_kernel_range_1,shared_dedisperse_kernel_range_2,shared_dedisperse_kernel_range_3
};

