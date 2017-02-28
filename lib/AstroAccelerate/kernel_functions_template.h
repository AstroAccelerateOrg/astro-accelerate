
static const int UNROLLS_ARRAY[] =  {{UNROLLS_RANGES}};
static const int SNUMREG_ARRAY[] =  {{SNUMREG_RANGES}};
static const int SDIVINT_ARRAY[] =  {{SDIVINT_RANGES}};
static const int SDIVINDM_ARRAY[] = {{SDIVINDM_RANGES}};

typedef void (*shared_dedisperse_kernel_FTYPE)(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep);

shared_dedisperse_kernel_FTYPE shared_dedisperse_kernel_FARRAY[] = {
	{KERNEL_FUNCTION_RANGES}
};

