#define ACCMAX 350
#define ACCSTEP 11
#define UNROLLS 16
#define SNUMREG 10
#define SDIVINT 8
#define SDIVINDM 40
#define SFDIVINDM 40.0f
#define CARD 0
#define NOPSSHIFT 5
#define NOPSLOOP 3
#define NDATAPERLOOP 1
#define BINDIVINT 6
#define BINDIVINF 32
#define CT 256
#define CF 2
#define NOPS 4.0
#define STATST 128
#define STATSLOOP 8

//Added by Karel Adamek
#define WARP 32 // doens change
#define HALF_WARP 16 // doesn't change
#define MSD_ELEM_PER_THREAD 8 // 1-32
#define MSD_WARPS_PER_BLOCK 16 //1-32 paddle 2
#define THR_ELEM_PER_THREAD 4 // 2-32 step 1
#define THR_WARPS_PER_BLOCK 4 // 1-32
#define PD_NTHREADS 512 // 64-1024 step 32
#define PD_NWINDOWS 2 // 1-8 step 1
#define PD_MAXTAPS 16 // doesn't change
#define PD_SMEM_SIZE 1280 // delete it
#define PD_FIR_ACTIVE_WARPS 2 // 2-32
#define PD_FIR_NWINDOWS 2 // 1-8
