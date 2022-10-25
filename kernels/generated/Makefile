include ../Make.inc 
ibit=64
IFLAGS = -DINDEXTYPE=int$(ibit)_t
OMPFLAGS = -fopenmp
PTFLAGS = $(OMPFLAGS) -DPTTIME -DLDB -DNTHREADS=24
SFLAGS = 
INC=$(INCSdir)/kernels.h 
#generated headers 
gcnGINC=$(GENINCdir)/sgkernels_gcn.h  
spmmGINC=$(GENINCdir)/sgkernels_spmm.h  
sigmoidGINC=$(GENINCdir)/sgkernels_sigmoid.h  
tdistGINC=$(GENINCdir)/sgkernels_tdist.h  
ARCHIVER = ar 
ARFLAGS = r
RANLIB = echo 
SDMMlib = $(LIBdir)/slibgfusedMM_sequential.a 
SDMMlib_pt = $(LIBdir)/slibgfusedMM_pt.a 

all: slib slib_pt 

sobj = $(BINdir)/skernels.o $(GENSRCdir)/sgfusedMM_K8_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K8_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K16_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K16_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K24_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K24_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K32_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K32_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K40_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K40_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K48_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K48_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K56_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K56_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K64_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K64_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K72_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K72_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K80_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K80_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K88_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K88_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K96_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K96_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K104_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K104_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K112_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K112_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K120_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K120_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K128_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K128_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K136_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K136_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K144_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K144_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K152_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K152_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K160_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K160_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K168_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K168_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K176_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K176_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K184_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K184_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K192_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K192_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K200_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K200_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K208_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K208_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K216_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K216_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K224_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K224_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K232_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K232_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K240_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K240_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K248_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K248_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K256_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K256_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K264_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K264_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K272_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K272_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K280_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K280_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K288_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K288_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K296_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K296_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K304_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K304_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K312_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K312_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K320_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K320_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K328_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K328_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K336_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K336_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K344_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K344_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K352_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K352_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K360_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K360_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K368_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K368_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K376_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K376_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K384_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K384_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K392_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K392_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K400_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K400_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K408_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K408_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K416_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K416_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K424_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K424_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K432_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K432_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K440_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K440_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K448_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K448_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K456_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K456_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K464_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K464_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K472_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K472_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K480_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K480_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K488_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K488_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K496_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K496_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K504_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K504_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K512_gcn_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K512_gcn_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K8_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K8_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K16_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K16_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K24_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K24_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K32_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K32_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K40_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K40_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K48_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K48_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K56_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K56_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K64_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K64_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K72_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K72_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K80_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K80_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K88_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K88_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K96_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K96_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K104_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K104_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K112_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K112_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K120_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K120_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K128_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K128_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K136_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K136_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K144_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K144_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K152_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K152_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K160_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K160_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K168_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K168_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K176_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K176_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K184_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K184_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K192_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K192_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K200_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K200_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K208_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K208_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K216_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K216_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K224_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K224_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K232_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K232_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K240_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K240_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K248_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K248_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K256_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K256_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K264_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K264_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K272_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K272_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K280_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K280_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K288_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K288_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K296_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K296_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K304_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K304_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K312_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K312_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K320_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K320_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K328_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K328_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K336_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K336_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K344_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K344_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K352_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K352_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K360_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K360_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K368_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K368_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K376_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K376_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K384_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K384_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K392_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K392_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K400_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K400_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K408_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K408_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K416_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K416_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K424_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K424_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K432_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K432_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K440_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K440_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K448_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K448_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K456_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K456_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K464_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K464_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K472_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K472_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K480_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K480_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K488_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K488_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K496_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K496_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K504_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K504_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K512_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K512_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K8_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K8_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K16_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K16_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K24_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K24_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K32_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K32_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K40_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K40_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K48_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K48_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K56_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K56_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K64_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K64_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K72_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K72_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K80_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K80_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K88_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K88_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K96_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K96_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K104_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K104_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K112_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K112_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K120_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K120_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K128_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K128_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K136_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K136_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K144_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K144_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K152_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K152_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K160_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K160_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K168_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K168_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K176_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K176_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K184_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K184_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K192_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K192_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K200_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K200_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K208_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K208_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K216_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K216_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K224_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K224_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K232_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K232_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K240_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K240_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K248_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K248_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K256_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K256_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K264_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K264_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K272_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K272_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K280_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K280_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K288_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K288_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K296_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K296_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K304_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K304_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K312_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K312_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K320_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K320_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K328_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K328_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K336_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K336_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K344_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K344_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K352_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K352_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K360_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K360_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K368_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K368_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K376_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K376_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K384_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K384_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K392_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K392_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K400_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K400_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K408_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K408_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K416_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K416_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K424_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K424_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K432_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K432_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K440_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K440_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K448_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K448_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K456_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K456_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K464_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K464_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K472_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K472_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K480_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K480_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K488_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K488_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K496_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K496_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K504_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K504_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K512_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K512_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K8_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K8_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K16_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K16_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K24_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K24_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K32_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K32_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K40_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K40_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K48_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K48_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K56_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K56_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K64_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K64_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K72_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K72_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K80_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K80_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K88_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K88_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K96_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K96_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K104_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K104_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K112_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K112_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K120_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K120_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K128_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K128_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K136_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K136_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K144_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K144_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K152_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K152_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K160_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K160_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K168_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K168_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K176_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K176_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K184_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K184_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K192_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K192_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K200_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K200_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K208_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K208_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K216_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K216_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K224_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K224_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K232_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K232_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K240_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K240_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K248_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K248_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K256_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K256_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K264_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K264_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K272_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K272_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K280_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K280_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K288_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K288_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K296_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K296_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K304_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K304_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K312_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K312_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K320_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K320_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K328_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K328_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K336_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K336_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K344_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K344_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K352_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K352_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K360_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K360_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K368_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K368_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K376_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K376_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K384_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K384_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K392_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K392_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K400_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K400_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K408_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K408_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K416_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K416_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K424_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K424_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K432_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K432_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K440_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K440_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K448_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K448_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K456_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K456_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K464_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K464_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K472_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K472_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K480_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K480_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K488_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K488_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K496_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K496_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K504_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K504_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K512_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K512_tdist_b1_csr.o

slib: $(LIBdir)/slib.grd 

$(LIBdir)/slib.grd: $(sobj)
	$(ARCHIVER) $(ARFLAGS) $(SDMMlib) $(sobj) 
	$(RANLIB) $(SDMMlib)  
	touch $(LIBdir)/slib.grd 

$(BINdir)/skernels.o : $(SRCdir)/kernels.c 
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -I$(GENINCdir) -o $@ -c $(SRCdir)/kernels.c 

$(GENSRCdir)/sgfusedMM_K8_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K8_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K8_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K8_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K8_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K8_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K16_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K16_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K16_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K16_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K16_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K16_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K24_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K24_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K24_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K24_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K24_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K24_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K32_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K32_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K32_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K32_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K32_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K32_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K40_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K40_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K40_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K40_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K40_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K40_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K48_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K48_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K48_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K48_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K48_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K48_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K56_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K56_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K56_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K56_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K56_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K56_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K64_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K64_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K64_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K64_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K64_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K64_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K72_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K72_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K72_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K72_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K72_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K72_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K80_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K80_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K80_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K80_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K80_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K80_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K88_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K88_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K88_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K88_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K88_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K88_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K96_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K96_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K96_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K96_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K96_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K96_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K104_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K104_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K104_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K104_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K104_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K104_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K112_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K112_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K112_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K112_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K112_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K112_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K120_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K120_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K120_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K120_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K120_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K120_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K128_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K128_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K128_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K128_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K128_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K128_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K136_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K136_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K136_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K136_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K136_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K136_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K144_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K144_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K144_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K144_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K144_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K144_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K152_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K152_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K152_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K152_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K152_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K152_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K160_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K160_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K160_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K160_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K160_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K160_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K168_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K168_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K168_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K168_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K168_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K168_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K176_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K176_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K176_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K176_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K176_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K176_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K184_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K184_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K184_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K184_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K184_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K184_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K192_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K192_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K192_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K192_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K192_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K192_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K200_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K200_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K200_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K200_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K200_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K200_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K208_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K208_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K208_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K208_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K208_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K208_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K216_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K216_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K216_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K216_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K216_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K216_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K224_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K224_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K224_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K224_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K224_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K224_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K232_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K232_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K232_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K232_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K232_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K232_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K240_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K240_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K240_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K240_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K240_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K240_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K248_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K248_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K248_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K248_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K248_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K248_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K256_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K256_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K256_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K256_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K256_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K256_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K264_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K264_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K264_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K264_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K264_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K264_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K272_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K272_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K272_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K272_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K272_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K272_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K280_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K280_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K280_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K280_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K280_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K280_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K288_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K288_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K288_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K288_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K288_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K288_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K296_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K296_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K296_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K296_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K296_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K296_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K304_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K304_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K304_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K304_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K304_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K304_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K312_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K312_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K312_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K312_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K312_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K312_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K320_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K320_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K320_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K320_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K320_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K320_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K328_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K328_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K328_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K328_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K328_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K328_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K336_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K336_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K336_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K336_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K336_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K336_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K344_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K344_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K344_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K344_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K344_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K344_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K352_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K352_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K352_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K352_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K352_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K352_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K360_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K360_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K360_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K360_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K360_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K360_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K368_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K368_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K368_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K368_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K368_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K368_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K376_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K376_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K376_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K376_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K376_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K376_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K384_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K384_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K384_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K384_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K384_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K384_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K392_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K392_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K392_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K392_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K392_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K392_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K400_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K400_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K400_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K400_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K400_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K400_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K408_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K408_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K408_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K408_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K408_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K408_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K416_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K416_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K416_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K416_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K416_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K416_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K424_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K424_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K424_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K424_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K424_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K424_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K432_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K432_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K432_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K432_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K432_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K432_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K440_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K440_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K440_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K440_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K440_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K440_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K448_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K448_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K448_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K448_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K448_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K448_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K456_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K456_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K456_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K456_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K456_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K456_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K464_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K464_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K464_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K464_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K464_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K464_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K472_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K472_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K472_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K472_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K472_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K472_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K480_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K480_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K480_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K480_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K480_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K480_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K488_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K488_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K488_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K488_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K488_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K488_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K496_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K496_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K496_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K496_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K496_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K496_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K504_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K504_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K504_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K504_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K504_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K504_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K512_gcn_b0_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K512_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K512_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K512_gcn_b1_csr.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K512_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K512_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K8_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K8_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K8_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K8_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K8_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K8_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K16_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K16_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K16_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K16_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K16_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K16_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K24_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K24_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K24_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K24_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K24_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K24_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K32_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K32_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K32_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K32_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K32_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K32_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K40_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K40_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K40_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K40_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K40_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K40_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K48_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K48_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K48_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K48_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K48_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K48_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K56_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K56_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K56_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K56_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K56_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K56_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K64_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K64_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K64_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K64_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K64_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K64_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K72_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K72_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K72_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K72_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K72_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K72_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K80_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K80_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K80_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K80_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K80_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K80_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K88_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K88_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K88_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K88_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K88_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K88_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K96_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K96_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K96_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K96_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K96_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K96_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K104_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K104_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K104_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K104_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K104_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K104_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K112_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K112_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K112_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K112_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K112_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K112_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K120_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K120_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K120_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K120_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K120_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K120_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K128_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K128_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K128_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K128_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K128_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K128_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K136_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K136_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K136_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K136_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K136_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K136_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K144_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K144_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K144_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K144_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K144_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K144_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K152_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K152_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K152_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K152_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K152_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K152_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K160_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K160_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K160_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K160_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K160_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K160_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K168_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K168_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K168_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K168_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K168_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K168_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K176_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K176_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K176_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K176_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K176_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K176_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K184_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K184_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K184_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K184_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K184_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K184_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K192_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K192_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K192_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K192_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K192_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K192_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K200_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K200_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K200_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K200_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K200_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K200_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K208_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K208_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K208_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K208_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K208_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K208_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K216_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K216_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K216_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K216_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K216_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K216_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K224_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K224_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K224_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K224_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K224_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K224_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K232_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K232_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K232_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K232_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K232_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K232_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K240_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K240_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K240_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K240_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K240_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K240_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K248_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K248_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K248_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K248_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K248_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K248_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K256_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K256_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K256_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K256_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K256_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K256_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K264_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K264_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K264_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K264_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K264_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K264_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K272_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K272_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K272_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K272_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K272_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K272_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K280_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K280_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K280_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K280_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K280_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K280_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K288_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K288_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K288_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K288_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K288_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K288_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K296_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K296_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K296_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K296_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K296_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K296_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K304_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K304_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K304_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K304_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K304_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K304_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K312_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K312_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K312_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K312_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K312_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K312_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K320_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K320_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K320_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K320_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K320_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K320_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K328_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K328_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K328_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K328_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K328_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K328_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K336_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K336_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K336_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K336_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K336_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K336_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K344_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K344_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K344_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K344_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K344_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K344_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K352_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K352_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K352_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K352_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K352_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K352_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K360_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K360_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K360_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K360_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K360_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K360_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K368_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K368_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K368_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K368_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K368_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K368_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K376_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K376_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K376_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K376_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K376_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K376_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K384_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K384_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K384_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K384_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K384_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K384_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K392_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K392_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K392_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K392_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K392_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K392_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K400_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K400_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K400_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K400_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K400_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K400_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K408_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K408_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K408_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K408_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K408_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K408_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K416_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K416_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K416_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K416_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K416_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K416_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K424_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K424_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K424_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K424_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K424_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K424_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K432_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K432_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K432_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K432_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K432_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K432_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K440_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K440_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K440_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K440_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K440_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K440_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K448_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K448_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K448_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K448_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K448_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K448_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K456_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K456_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K456_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K456_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K456_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K456_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K464_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K464_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K464_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K464_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K464_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K464_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K472_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K472_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K472_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K472_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K472_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K472_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K480_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K480_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K480_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K480_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K480_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K480_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K488_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K488_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K488_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K488_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K488_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K488_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K496_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K496_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K496_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K496_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K496_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K496_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K504_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K504_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K504_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K504_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K504_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K504_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K512_spmm_b0_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K512_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K512_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K512_spmm_b1_csr.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K512_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K512_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K8_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K8_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K8_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K8_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K8_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K8_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K16_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K16_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K16_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K16_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K16_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K16_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K24_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K24_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K24_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K24_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K24_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K24_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K32_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K32_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K32_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K32_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K32_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K32_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K40_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K40_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K40_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K40_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K40_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K40_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K48_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K48_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K48_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K48_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K48_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K48_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K56_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K56_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K56_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K56_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K56_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K56_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K64_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K64_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K64_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K64_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K64_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K64_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K72_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K72_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K72_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K72_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K72_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K72_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K80_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K80_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K80_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K80_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K80_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K80_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K88_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K88_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K88_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K88_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K88_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K88_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K96_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K96_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K96_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K96_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K96_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K96_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K104_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K104_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K104_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K104_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K104_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K104_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K112_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K112_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K112_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K112_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K112_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K112_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K120_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K120_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K120_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K120_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K120_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K120_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K128_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K128_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K128_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K128_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K128_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K128_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K136_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K136_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K136_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K136_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K136_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K136_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K144_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K144_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K144_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K144_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K144_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K144_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K152_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K152_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K152_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K152_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K152_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K152_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K160_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K160_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K160_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K160_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K160_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K160_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K168_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K168_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K168_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K168_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K168_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K168_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K176_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K176_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K176_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K176_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K176_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K176_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K184_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K184_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K184_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K184_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K184_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K184_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K192_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K192_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K192_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K192_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K192_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K192_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K200_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K200_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K200_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K200_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K200_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K200_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K208_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K208_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K208_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K208_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K208_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K208_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K216_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K216_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K216_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K216_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K216_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K216_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K224_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K224_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K224_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K224_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K224_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K224_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K232_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K232_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K232_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K232_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K232_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K232_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K240_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K240_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K240_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K240_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K240_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K240_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K248_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K248_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K248_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K248_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K248_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K248_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K256_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K256_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K256_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K256_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K256_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K256_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K264_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K264_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K264_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K264_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K264_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K264_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K272_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K272_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K272_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K272_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K272_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K272_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K280_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K280_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K280_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K280_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K280_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K280_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K288_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K288_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K288_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K288_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K288_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K288_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K296_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K296_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K296_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K296_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K296_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K296_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K304_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K304_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K304_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K304_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K304_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K304_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K312_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K312_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K312_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K312_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K312_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K312_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K320_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K320_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K320_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K320_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K320_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K320_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K328_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K328_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K328_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K328_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K328_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K328_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K336_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K336_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K336_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K336_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K336_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K336_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K344_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K344_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K344_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K344_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K344_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K344_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K352_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K352_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K352_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K352_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K352_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K352_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K360_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K360_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K360_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K360_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K360_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K360_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K368_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K368_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K368_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K368_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K368_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K368_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K376_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K376_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K376_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K376_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K376_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K376_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K384_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K384_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K384_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K384_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K384_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K384_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K392_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K392_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K392_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K392_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K392_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K392_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K400_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K400_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K400_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K400_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K400_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K400_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K408_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K408_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K408_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K408_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K408_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K408_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K416_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K416_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K416_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K416_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K416_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K416_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K424_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K424_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K424_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K424_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K424_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K424_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K432_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K432_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K432_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K432_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K432_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K432_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K440_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K440_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K440_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K440_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K440_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K440_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K448_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K448_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K448_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K448_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K448_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K448_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K456_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K456_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K456_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K456_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K456_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K456_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K464_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K464_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K464_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K464_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K464_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K464_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K472_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K472_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K472_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K472_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K472_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K472_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K480_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K480_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K480_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K480_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K480_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K480_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K488_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K488_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K488_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K488_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K488_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K488_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K496_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K496_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K496_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K496_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K496_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K496_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K504_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K504_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K504_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K504_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K504_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K504_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K512_sigmoid_b0_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K512_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K512_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K512_sigmoid_b1_csr.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K512_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K512_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K8_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K8_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K8_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K8_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K8_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K8_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K16_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K16_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K16_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K16_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K16_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K16_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K24_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K24_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K24_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K24_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K24_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K24_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K32_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K32_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K32_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K32_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K32_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K32_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K40_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K40_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K40_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K40_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K40_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K40_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K48_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K48_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K48_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K48_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K48_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K48_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K56_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K56_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K56_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K56_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K56_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K56_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K64_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K64_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K64_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K64_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K64_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K64_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K72_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K72_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K72_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K72_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K72_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K72_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K80_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K80_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K80_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K80_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K80_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K80_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K88_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K88_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K88_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K88_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K88_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K88_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K96_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K96_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K96_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K96_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K96_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K96_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K104_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K104_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K104_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K104_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K104_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K104_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K112_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K112_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K112_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K112_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K112_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K112_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K120_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K120_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K120_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K120_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K120_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K120_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K128_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K128_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K128_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K128_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K128_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K128_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K136_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K136_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K136_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K136_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K136_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K136_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K144_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K144_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K144_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K144_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K144_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K144_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K152_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K152_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K152_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K152_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K152_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K152_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K160_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K160_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K160_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K160_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K160_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K160_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K168_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K168_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K168_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K168_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K168_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K168_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K176_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K176_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K176_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K176_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K176_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K176_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K184_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K184_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K184_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K184_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K184_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K184_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K192_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K192_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K192_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K192_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K192_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K192_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K200_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K200_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K200_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K200_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K200_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K200_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K208_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K208_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K208_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K208_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K208_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K208_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K216_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K216_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K216_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K216_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K216_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K216_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K224_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K224_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K224_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K224_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K224_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K224_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K232_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K232_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K232_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K232_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K232_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K232_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K240_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K240_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K240_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K240_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K240_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K240_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K248_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K248_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K248_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K248_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K248_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K248_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K256_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K256_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K256_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K256_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K256_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K256_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K264_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K264_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K264_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K264_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K264_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K264_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K272_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K272_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K272_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K272_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K272_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K272_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K280_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K280_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K280_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K280_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K280_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K280_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K288_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K288_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K288_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K288_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K288_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K288_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K296_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K296_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K296_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K296_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K296_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K296_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K304_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K304_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K304_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K304_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K304_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K304_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K312_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K312_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K312_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K312_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K312_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K312_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K320_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K320_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K320_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K320_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K320_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K320_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K328_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K328_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K328_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K328_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K328_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K328_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K336_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K336_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K336_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K336_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K336_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K336_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K344_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K344_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K344_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K344_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K344_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K344_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K352_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K352_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K352_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K352_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K352_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K352_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K360_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K360_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K360_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K360_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K360_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K360_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K368_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K368_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K368_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K368_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K368_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K368_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K376_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K376_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K376_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K376_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K376_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K376_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K384_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K384_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K384_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K384_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K384_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K384_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K392_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K392_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K392_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K392_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K392_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K392_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K400_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K400_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K400_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K400_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K400_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K400_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K408_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K408_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K408_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K408_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K408_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K408_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K416_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K416_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K416_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K416_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K416_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K416_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K424_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K424_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K424_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K424_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K424_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K424_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K432_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K432_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K432_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K432_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K432_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K432_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K440_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K440_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K440_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K440_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K440_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K440_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K448_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K448_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K448_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K448_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K448_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K448_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K456_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K456_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K456_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K456_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K456_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K456_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K464_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K464_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K464_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K464_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K464_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K464_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K472_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K472_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K472_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K472_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K472_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K472_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K480_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K480_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K480_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K480_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K480_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K480_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K488_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K488_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K488_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K488_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K488_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K488_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K496_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K496_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K496_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K496_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K496_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K496_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K504_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K504_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K504_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K504_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K504_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K504_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K512_tdist_b0_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K512_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K512_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K512_tdist_b1_csr.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K512_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(SFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K512_tdist_csr.c

sobj_pt = $(BINdir)/skernels_pt.o $(GENSRCdir)/sgfusedMM_K8_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K8_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K16_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K16_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K24_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K24_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K32_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K32_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K40_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K40_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K48_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K48_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K56_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K56_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K64_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K64_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K72_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K72_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K80_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K80_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K88_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K88_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K96_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K96_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K104_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K104_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K112_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K112_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K120_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K120_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K128_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K128_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K136_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K136_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K144_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K144_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K152_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K152_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K160_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K160_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K168_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K168_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K176_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K176_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K184_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K184_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K192_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K192_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K200_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K200_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K208_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K208_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K216_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K216_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K224_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K224_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K232_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K232_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K240_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K240_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K248_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K248_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K256_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K256_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K264_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K264_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K272_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K272_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K280_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K280_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K288_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K288_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K296_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K296_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K304_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K304_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K312_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K312_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K320_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K320_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K328_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K328_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K336_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K336_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K344_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K344_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K352_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K352_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K360_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K360_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K368_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K368_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K376_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K376_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K384_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K384_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K392_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K392_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K400_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K400_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K408_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K408_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K416_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K416_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K424_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K424_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K432_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K432_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K440_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K440_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K448_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K448_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K456_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K456_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K464_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K464_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K472_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K472_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K480_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K480_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K488_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K488_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K496_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K496_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K504_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K504_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K512_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K512_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K8_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K8_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K16_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K16_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K24_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K24_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K32_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K32_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K40_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K40_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K48_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K48_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K56_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K56_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K64_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K64_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K72_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K72_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K80_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K80_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K88_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K88_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K96_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K96_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K104_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K104_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K112_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K112_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K120_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K120_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K128_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K128_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K136_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K136_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K144_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K144_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K152_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K152_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K160_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K160_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K168_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K168_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K176_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K176_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K184_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K184_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K192_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K192_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K200_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K200_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K208_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K208_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K216_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K216_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K224_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K224_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K232_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K232_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K240_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K240_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K248_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K248_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K256_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K256_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K264_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K264_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K272_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K272_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K280_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K280_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K288_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K288_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K296_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K296_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K304_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K304_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K312_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K312_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K320_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K320_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K328_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K328_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K336_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K336_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K344_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K344_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K352_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K352_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K360_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K360_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K368_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K368_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K376_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K376_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K384_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K384_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K392_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K392_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K400_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K400_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K408_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K408_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K416_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K416_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K424_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K424_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K432_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K432_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K440_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K440_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K448_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K448_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K456_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K456_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K464_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K464_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K472_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K472_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K480_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K480_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K488_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K488_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K496_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K496_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K504_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K504_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K512_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K512_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K8_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K8_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K16_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K16_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K24_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K24_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K32_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K32_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K40_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K40_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K48_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K48_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K56_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K56_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K64_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K64_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K72_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K72_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K80_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K80_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K88_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K88_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K96_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K96_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K104_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K104_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K112_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K112_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K120_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K120_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K128_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K128_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K136_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K136_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K144_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K144_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K152_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K152_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K160_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K160_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K168_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K168_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K176_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K176_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K184_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K184_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K192_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K192_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K200_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K200_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K208_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K208_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K216_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K216_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K224_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K224_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K232_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K232_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K240_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K240_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K248_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K248_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K256_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K256_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K264_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K264_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K272_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K272_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K280_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K280_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K288_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K288_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K296_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K296_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K304_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K304_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K312_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K312_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K320_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K320_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K328_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K328_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K336_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K336_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K344_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K344_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K352_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K352_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K360_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K360_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K368_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K368_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K376_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K376_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K384_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K384_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K392_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K392_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K400_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K400_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K408_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K408_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K416_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K416_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K424_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K424_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K432_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K432_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K440_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K440_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K448_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K448_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K456_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K456_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K464_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K464_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K472_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K472_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K480_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K480_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K488_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K488_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K496_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K496_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K504_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K504_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K512_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K512_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K8_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K8_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K16_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K16_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K24_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K24_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K32_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K32_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K40_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K40_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K48_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K48_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K56_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K56_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K64_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K64_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K72_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K72_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K80_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K80_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K88_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K88_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K96_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K96_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K104_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K104_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K112_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K112_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K120_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K120_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K128_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K128_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K136_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K136_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K144_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K144_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K152_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K152_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K160_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K160_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K168_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K168_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K176_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K176_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K184_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K184_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K192_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K192_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K200_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K200_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K208_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K208_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K216_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K216_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K224_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K224_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K232_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K232_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K240_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K240_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K248_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K248_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K256_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K256_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K264_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K264_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K272_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K272_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K280_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K280_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K288_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K288_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K296_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K296_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K304_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K304_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K312_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K312_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K320_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K320_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K328_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K328_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K336_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K336_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K344_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K344_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K352_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K352_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K360_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K360_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K368_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K368_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K376_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K376_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K384_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K384_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K392_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K392_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K400_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K400_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K408_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K408_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K416_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K416_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K424_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K424_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K432_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K432_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K440_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K440_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K448_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K448_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K456_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K456_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K464_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K464_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K472_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K472_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K480_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K480_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K488_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K488_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K496_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K496_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K504_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K504_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K512_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K512_tdist_b1_csr_pt.o

slib_pt: $(LIBdir)/slib_pt.grd 

$(LIBdir)/slib_pt.grd: $(sobj_pt)
	$(ARCHIVER) $(ARFLAGS) $(SDMMlib_pt) $(sobj_pt) 
	$(RANLIB) $(SDMMlib_pt)  
	touch $(LIBdir)/slib_pt.grd 

$(BINdir)/skernels_pt.o : $(SRCdir)/kernels.c 
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -I$(GENINCdir) -o $@ -c $(SRCdir)/kernels.c 

$(GENSRCdir)/sgfusedMM_K8_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K8_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K8_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K8_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K8_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K8_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K16_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K16_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K16_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K16_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K16_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K16_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K24_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K24_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K24_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K24_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K24_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K24_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K32_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K32_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K32_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K32_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K32_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K32_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K40_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K40_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K40_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K40_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K40_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K40_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K48_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K48_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K48_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K48_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K48_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K48_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K56_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K56_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K56_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K56_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K56_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K56_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K64_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K64_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K64_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K64_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K64_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K64_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K72_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K72_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K72_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K72_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K72_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K72_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K80_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K80_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K80_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K80_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K80_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K80_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K88_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K88_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K88_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K88_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K88_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K88_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K96_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K96_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K96_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K96_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K96_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K96_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K104_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K104_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K104_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K104_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K104_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K104_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K112_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K112_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K112_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K112_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K112_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K112_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K120_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K120_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K120_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K120_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K120_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K120_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K128_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K128_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K128_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K128_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K128_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K128_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K136_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K136_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K136_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K136_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K136_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K136_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K144_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K144_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K144_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K144_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K144_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K144_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K152_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K152_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K152_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K152_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K152_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K152_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K160_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K160_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K160_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K160_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K160_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K160_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K168_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K168_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K168_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K168_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K168_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K168_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K176_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K176_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K176_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K176_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K176_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K176_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K184_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K184_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K184_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K184_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K184_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K184_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K192_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K192_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K192_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K192_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K192_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K192_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K200_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K200_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K200_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K200_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K200_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K200_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K208_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K208_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K208_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K208_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K208_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K208_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K216_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K216_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K216_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K216_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K216_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K216_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K224_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K224_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K224_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K224_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K224_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K224_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K232_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K232_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K232_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K232_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K232_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K232_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K240_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K240_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K240_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K240_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K240_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K240_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K248_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K248_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K248_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K248_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K248_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K248_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K256_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K256_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K256_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K256_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K256_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K256_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K264_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K264_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K264_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K264_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K264_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K264_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K272_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K272_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K272_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K272_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K272_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K272_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K280_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K280_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K280_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K280_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K280_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K280_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K288_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K288_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K288_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K288_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K288_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K288_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K296_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K296_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K296_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K296_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K296_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K296_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K304_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K304_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K304_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K304_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K304_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K304_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K312_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K312_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K312_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K312_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K312_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K312_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K320_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K320_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K320_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K320_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K320_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K320_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K328_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K328_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K328_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K328_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K328_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K328_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K336_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K336_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K336_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K336_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K336_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K336_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K344_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K344_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K344_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K344_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K344_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K344_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K352_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K352_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K352_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K352_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K352_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K352_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K360_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K360_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K360_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K360_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K360_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K360_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K368_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K368_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K368_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K368_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K368_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K368_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K376_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K376_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K376_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K376_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K376_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K376_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K384_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K384_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K384_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K384_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K384_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K384_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K392_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K392_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K392_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K392_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K392_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K392_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K400_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K400_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K400_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K400_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K400_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K400_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K408_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K408_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K408_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K408_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K408_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K408_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K416_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K416_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K416_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K416_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K416_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K416_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K424_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K424_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K424_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K424_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K424_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K424_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K432_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K432_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K432_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K432_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K432_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K432_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K440_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K440_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K440_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K440_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K440_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K440_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K448_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K448_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K448_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K448_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K448_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K448_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K456_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K456_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K456_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K456_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K456_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K456_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K464_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K464_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K464_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K464_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K464_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K464_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K472_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K472_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K472_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K472_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K472_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K472_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K480_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K480_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K480_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K480_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K480_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K480_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K488_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K488_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K488_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K488_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K488_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K488_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K496_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K496_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K496_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K496_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K496_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K496_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K504_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K504_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K504_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K504_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K504_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K504_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K512_gcn_b0_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K512_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K512_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K512_gcn_b1_csr_pt.o : $(gcnGINC) \
   $(GENSRCdir)/sgfusedMM_K512_gcn_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K512_gcn_csr.c
$(GENSRCdir)/sgfusedMM_K8_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K8_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K8_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K8_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K8_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K8_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K16_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K16_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K16_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K16_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K16_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K16_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K24_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K24_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K24_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K24_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K24_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K24_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K32_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K32_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K32_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K32_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K32_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K32_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K40_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K40_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K40_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K40_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K40_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K40_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K48_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K48_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K48_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K48_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K48_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K48_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K56_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K56_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K56_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K56_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K56_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K56_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K64_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K64_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K64_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K64_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K64_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K64_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K72_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K72_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K72_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K72_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K72_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K72_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K80_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K80_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K80_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K80_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K80_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K80_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K88_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K88_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K88_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K88_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K88_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K88_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K96_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K96_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K96_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K96_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K96_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K96_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K104_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K104_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K104_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K104_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K104_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K104_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K112_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K112_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K112_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K112_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K112_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K112_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K120_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K120_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K120_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K120_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K120_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K120_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K128_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K128_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K128_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K128_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K128_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K128_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K136_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K136_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K136_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K136_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K136_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K136_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K144_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K144_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K144_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K144_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K144_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K144_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K152_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K152_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K152_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K152_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K152_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K152_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K160_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K160_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K160_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K160_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K160_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K160_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K168_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K168_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K168_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K168_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K168_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K168_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K176_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K176_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K176_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K176_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K176_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K176_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K184_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K184_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K184_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K184_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K184_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K184_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K192_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K192_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K192_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K192_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K192_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K192_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K200_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K200_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K200_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K200_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K200_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K200_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K208_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K208_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K208_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K208_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K208_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K208_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K216_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K216_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K216_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K216_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K216_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K216_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K224_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K224_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K224_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K224_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K224_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K224_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K232_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K232_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K232_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K232_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K232_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K232_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K240_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K240_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K240_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K240_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K240_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K240_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K248_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K248_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K248_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K248_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K248_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K248_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K256_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K256_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K256_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K256_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K256_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K256_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K264_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K264_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K264_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K264_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K264_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K264_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K272_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K272_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K272_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K272_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K272_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K272_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K280_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K280_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K280_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K280_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K280_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K280_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K288_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K288_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K288_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K288_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K288_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K288_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K296_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K296_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K296_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K296_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K296_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K296_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K304_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K304_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K304_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K304_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K304_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K304_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K312_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K312_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K312_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K312_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K312_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K312_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K320_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K320_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K320_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K320_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K320_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K320_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K328_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K328_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K328_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K328_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K328_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K328_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K336_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K336_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K336_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K336_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K336_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K336_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K344_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K344_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K344_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K344_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K344_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K344_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K352_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K352_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K352_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K352_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K352_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K352_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K360_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K360_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K360_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K360_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K360_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K360_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K368_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K368_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K368_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K368_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K368_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K368_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K376_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K376_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K376_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K376_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K376_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K376_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K384_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K384_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K384_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K384_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K384_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K384_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K392_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K392_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K392_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K392_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K392_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K392_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K400_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K400_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K400_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K400_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K400_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K400_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K408_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K408_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K408_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K408_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K408_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K408_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K416_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K416_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K416_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K416_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K416_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K416_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K424_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K424_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K424_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K424_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K424_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K424_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K432_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K432_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K432_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K432_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K432_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K432_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K440_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K440_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K440_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K440_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K440_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K440_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K448_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K448_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K448_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K448_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K448_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K448_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K456_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K456_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K456_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K456_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K456_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K456_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K464_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K464_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K464_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K464_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K464_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K464_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K472_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K472_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K472_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K472_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K472_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K472_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K480_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K480_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K480_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K480_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K480_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K480_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K488_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K488_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K488_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K488_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K488_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K488_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K496_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K496_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K496_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K496_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K496_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K496_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K504_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K504_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K504_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K504_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K504_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K504_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K512_spmm_b0_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K512_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K512_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K512_spmm_b1_csr_pt.o : $(spmmGINC) \
   $(GENSRCdir)/sgfusedMM_K512_spmm_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K512_spmm_csr.c
$(GENSRCdir)/sgfusedMM_K8_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K8_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K8_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K8_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K8_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K8_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K16_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K16_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K16_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K16_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K16_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K16_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K24_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K24_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K24_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K24_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K24_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K24_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K32_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K32_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K32_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K32_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K32_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K32_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K40_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K40_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K40_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K40_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K40_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K40_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K48_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K48_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K48_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K48_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K48_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K48_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K56_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K56_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K56_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K56_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K56_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K56_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K64_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K64_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K64_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K64_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K64_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K64_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K72_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K72_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K72_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K72_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K72_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K72_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K80_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K80_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K80_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K80_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K80_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K80_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K88_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K88_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K88_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K88_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K88_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K88_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K96_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K96_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K96_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K96_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K96_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K96_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K104_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K104_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K104_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K104_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K104_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K104_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K112_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K112_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K112_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K112_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K112_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K112_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K120_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K120_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K120_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K120_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K120_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K120_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K128_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K128_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K128_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K128_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K128_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K128_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K136_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K136_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K136_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K136_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K136_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K136_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K144_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K144_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K144_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K144_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K144_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K144_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K152_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K152_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K152_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K152_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K152_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K152_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K160_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K160_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K160_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K160_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K160_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K160_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K168_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K168_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K168_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K168_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K168_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K168_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K176_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K176_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K176_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K176_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K176_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K176_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K184_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K184_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K184_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K184_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K184_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K184_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K192_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K192_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K192_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K192_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K192_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K192_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K200_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K200_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K200_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K200_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K200_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K200_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K208_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K208_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K208_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K208_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K208_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K208_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K216_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K216_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K216_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K216_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K216_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K216_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K224_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K224_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K224_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K224_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K224_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K224_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K232_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K232_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K232_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K232_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K232_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K232_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K240_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K240_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K240_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K240_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K240_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K240_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K248_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K248_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K248_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K248_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K248_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K248_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K256_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K256_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K256_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K256_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K256_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K256_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K264_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K264_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K264_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K264_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K264_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K264_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K272_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K272_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K272_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K272_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K272_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K272_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K280_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K280_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K280_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K280_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K280_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K280_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K288_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K288_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K288_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K288_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K288_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K288_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K296_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K296_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K296_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K296_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K296_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K296_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K304_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K304_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K304_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K304_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K304_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K304_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K312_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K312_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K312_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K312_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K312_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K312_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K320_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K320_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K320_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K320_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K320_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K320_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K328_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K328_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K328_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K328_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K328_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K328_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K336_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K336_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K336_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K336_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K336_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K336_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K344_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K344_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K344_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K344_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K344_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K344_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K352_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K352_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K352_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K352_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K352_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K352_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K360_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K360_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K360_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K360_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K360_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K360_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K368_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K368_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K368_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K368_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K368_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K368_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K376_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K376_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K376_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K376_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K376_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K376_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K384_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K384_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K384_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K384_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K384_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K384_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K392_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K392_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K392_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K392_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K392_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K392_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K400_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K400_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K400_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K400_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K400_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K400_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K408_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K408_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K408_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K408_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K408_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K408_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K416_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K416_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K416_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K416_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K416_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K416_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K424_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K424_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K424_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K424_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K424_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K424_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K432_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K432_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K432_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K432_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K432_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K432_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K440_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K440_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K440_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K440_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K440_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K440_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K448_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K448_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K448_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K448_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K448_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K448_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K456_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K456_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K456_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K456_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K456_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K456_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K464_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K464_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K464_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K464_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K464_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K464_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K472_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K472_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K472_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K472_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K472_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K472_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K480_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K480_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K480_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K480_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K480_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K480_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K488_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K488_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K488_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K488_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K488_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K488_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K496_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K496_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K496_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K496_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K496_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K496_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K504_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K504_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K504_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K504_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K504_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K504_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K512_sigmoid_b0_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K512_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K512_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K512_sigmoid_b1_csr_pt.o : $(sigmoidGINC) \
   $(GENSRCdir)/sgfusedMM_K512_sigmoid_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K512_sigmoid_csr.c
$(GENSRCdir)/sgfusedMM_K8_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K8_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K8_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K8_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K8_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K8_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K16_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K16_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K16_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K16_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K16_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K16_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K24_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K24_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K24_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K24_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K24_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K24_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K32_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K32_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K32_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K32_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K32_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K32_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K40_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K40_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K40_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K40_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K40_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K40_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K48_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K48_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K48_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K48_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K48_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K48_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K56_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K56_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K56_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K56_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K56_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K56_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K64_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K64_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K64_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K64_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K64_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K64_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K72_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K72_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K72_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K72_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K72_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K72_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K80_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K80_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K80_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K80_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K80_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K80_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K88_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K88_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K88_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K88_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K88_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K88_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K96_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K96_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K96_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K96_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K96_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K96_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K104_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K104_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K104_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K104_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K104_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K104_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K112_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K112_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K112_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K112_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K112_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K112_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K120_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K120_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K120_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K120_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K120_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K120_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K128_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K128_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K128_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K128_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K128_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K128_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K136_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K136_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K136_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K136_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K136_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K136_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K144_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K144_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K144_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K144_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K144_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K144_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K152_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K152_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K152_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K152_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K152_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K152_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K160_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K160_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K160_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K160_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K160_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K160_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K168_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K168_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K168_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K168_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K168_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K168_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K176_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K176_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K176_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K176_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K176_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K176_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K184_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K184_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K184_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K184_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K184_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K184_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K192_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K192_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K192_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K192_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K192_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K192_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K200_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K200_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K200_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K200_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K200_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K200_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K208_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K208_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K208_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K208_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K208_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K208_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K216_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K216_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K216_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K216_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K216_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K216_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K224_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K224_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K224_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K224_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K224_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K224_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K232_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K232_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K232_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K232_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K232_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K232_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K240_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K240_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K240_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K240_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K240_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K240_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K248_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K248_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K248_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K248_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K248_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K248_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K256_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K256_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K256_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K256_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K256_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K256_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K264_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K264_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K264_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K264_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K264_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K264_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K272_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K272_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K272_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K272_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K272_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K272_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K280_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K280_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K280_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K280_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K280_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K280_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K288_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K288_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K288_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K288_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K288_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K288_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K296_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K296_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K296_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K296_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K296_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K296_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K304_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K304_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K304_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K304_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K304_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K304_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K312_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K312_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K312_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K312_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K312_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K312_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K320_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K320_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K320_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K320_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K320_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K320_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K328_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K328_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K328_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K328_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K328_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K328_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K336_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K336_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K336_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K336_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K336_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K336_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K344_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K344_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K344_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K344_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K344_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K344_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K352_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K352_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K352_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K352_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K352_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K352_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K360_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K360_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K360_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K360_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K360_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K360_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K368_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K368_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K368_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K368_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K368_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K368_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K376_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K376_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K376_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K376_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K376_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K376_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K384_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K384_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K384_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K384_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K384_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K384_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K392_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K392_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K392_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K392_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K392_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K392_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K400_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K400_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K400_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K400_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K400_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K400_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K408_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K408_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K408_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K408_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K408_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K408_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K416_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K416_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K416_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K416_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K416_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K416_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K424_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K424_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K424_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K424_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K424_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K424_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K432_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K432_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K432_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K432_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K432_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K432_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K440_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K440_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K440_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K440_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K440_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K440_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K448_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K448_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K448_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K448_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K448_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K448_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K456_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K456_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K456_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K456_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K456_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K456_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K464_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K464_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K464_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K464_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K464_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K464_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K472_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K472_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K472_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K472_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K472_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K472_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K480_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K480_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K480_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K480_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K480_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K480_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K488_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K488_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K488_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K488_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K488_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K488_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K496_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K496_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K496_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K496_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K496_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K496_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K504_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K504_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K504_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K504_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K504_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K504_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K512_tdist_b0_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K512_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA0 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K512_tdist_csr.c
$(GENSRCdir)/sgfusedMM_K512_tdist_b1_csr_pt.o : $(tdistGINC) \
   $(GENSRCdir)/sgfusedMM_K512_tdist_csr.c
	$(KCC) $(KCCFLAGS) $(ARCHFLAGS) $(IFLAGS) $(PTFLAGS) -DSREAL \
        -DBETA1 -I$(SIMDdir) -o $@ -c $(GENSRCdir)/sgfusedMM_K512_tdist_csr.c



clean : 
	rm -f $(GENSRCdir)/s*.o $(LIBdir)/s*.a  $(LIBdir)/s*.grd 


