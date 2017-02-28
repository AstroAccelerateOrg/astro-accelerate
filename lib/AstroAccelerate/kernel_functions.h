
static const int UNROLLS_ARRAY[] =  {UNROLLS_0,UNROLLS_1,UNROLLS_2,UNROLLS_3,UNROLLS_4,UNROLLS_5,UNROLLS_6,UNROLLS_7};
static const int SNUMREG_ARRAY[] =  {SNUMREG_0,SNUMREG_1,SNUMREG_2,SNUMREG_3,SNUMREG_4,SNUMREG_5,SNUMREG_6,SNUMREG_7};
static const int SDIVINT_ARRAY[] =  {SDIVINT_0,SDIVINT_1,SDIVINT_2,SDIVINT_3,SDIVINT_4,SDIVINT_5,SDIVINT_6,SDIVINT_7};
static const int SDIVINDM_ARRAY[] = {SDIVINDM_0,SDIVINDM_1,SDIVINDM_2,SDIVINDM_3,SDIVINDM_4,SDIVINDM_5,SDIVINDM_6,SDIVINDM_7};

typedef void (*shared_dedisperse_kernel_FTYPE)(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep);

shared_dedisperse_kernel_FTYPE shared_dedisperse_kernel_FARRAY[] = {
	shared_dedisperse_kernel_range_0,shared_dedisperse_kernel_range_1,shared_dedisperse_kernel_range_2,shared_dedisperse_kernel_range_3,shared_dedisperse_kernel_range_4,shared_dedisperse_kernel_range_5,shared_dedisperse_kernel_range_6,shared_dedisperse_kernel_range_7
};

