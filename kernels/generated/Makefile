include ../Make.inc 
ibit=64
IFLAGS = -DINDEXTYPE=int$(ibit)_t
OMPFLAGS = -fopenmp
PTFLAGS = $(OMPFLAGS) -DPTTIME -DLDB -DNTHREADS=4
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
       $(GENSRCdir)/sgfusedMM_K8_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K8_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K16_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K16_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K24_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K24_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K32_spmm_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K32_spmm_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K8_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K8_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K16_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K16_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K24_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K24_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K32_sigmoid_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K32_sigmoid_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K8_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K8_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K16_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K16_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K24_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K24_tdist_b1_csr.o \
       $(GENSRCdir)/sgfusedMM_K32_tdist_b0_csr.o \
       $(GENSRCdir)/sgfusedMM_K32_tdist_b1_csr.o

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

sobj_pt = $(BINdir)/skernels_pt.o $(GENSRCdir)/sgfusedMM_K8_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K8_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K16_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K16_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K24_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K24_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K32_gcn_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K32_gcn_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K8_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K8_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K16_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K16_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K24_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K24_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K32_spmm_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K32_spmm_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K8_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K8_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K16_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K16_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K24_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K24_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K32_sigmoid_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K32_sigmoid_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K8_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K8_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K16_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K16_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K24_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K24_tdist_b1_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K32_tdist_b0_csr_pt.o \
          $(GENSRCdir)/sgfusedMM_K32_tdist_b1_csr_pt.o

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



clean : 
	rm -f $(GENSRCdir)/s*.o $(LIBdir)/s*.a  $(LIBdir)/s*.grd 


