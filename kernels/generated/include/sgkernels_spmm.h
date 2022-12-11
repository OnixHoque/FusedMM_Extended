#ifndef DG_SPMM_KERNEL_H
#define DG_SPMM_KERNEL_H
#define GVLEN 8 /* generated kernels VLEN, check with simd.h */
#define MAXDIM_SPMM 32 /* max value of dimension (k) */
#define KRUNTIME_SPMM 0
/*
 * NOTE: put the best K value after tuning here, needed when kruntime = 1 
 */
#define BESTK_SPMM 32
/*
 * function pointer type for generated kernels 
 */
/*
 * Kernels for beta, b0
 */

typedef void (*kern_sgfusedMM_spmm_b0_t) ( const char transa, const INDEXTYPE m, 
      const INDEXTYPE n, const INDEXTYPE k,const float alpha, const INDEXTYPE nnz, 
      const INDEXTYPE rows, const INDEXTYPE cols, const float *val, const INDEXTYPE *indx, 
      const INDEXTYPE *pntrb, const INDEXTYPE *pntre, const float *A, 
      const INDEXTYPE lda, const float *B, const INDEXTYPE ldb, 
      const float beta, float *C, const INDEXTYPE ldc);
/*
 * FIXME: add a check for VLEN in generator with simd.h 
 * VLEN = 8, MAX DIM(K) = 32
 */
void sgfusedMM_K8_spmm_b0_csr (const char transa, const INDEXTYPE m, const INDEXTYPE n, 
      const INDEXTYPE k,const float alpha, const INDEXTYPE nnz, 
      const INDEXTYPE rows, const INDEXTYPE cols, const float *val, const INDEXTYPE *indx, 
      const INDEXTYPE *pntrb, const INDEXTYPE *pntre, const float *A, 
      const INDEXTYPE lda, const float *B, const INDEXTYPE ldb, 
      const float beta, float *C, const INDEXTYPE ldc);
void sgfusedMM_K16_spmm_b0_csr (const char transa, const INDEXTYPE m, const INDEXTYPE n, 
      const INDEXTYPE k,const float alpha, const INDEXTYPE nnz, 
      const INDEXTYPE rows, const INDEXTYPE cols, const float *val, const INDEXTYPE *indx, 
      const INDEXTYPE *pntrb, const INDEXTYPE *pntre, const float *A, 
      const INDEXTYPE lda, const float *B, const INDEXTYPE ldb, 
      const float beta, float *C, const INDEXTYPE ldc);
void sgfusedMM_K24_spmm_b0_csr (const char transa, const INDEXTYPE m, const INDEXTYPE n, 
      const INDEXTYPE k,const float alpha, const INDEXTYPE nnz, 
      const INDEXTYPE rows, const INDEXTYPE cols, const float *val, const INDEXTYPE *indx, 
      const INDEXTYPE *pntrb, const INDEXTYPE *pntre, const float *A, 
      const INDEXTYPE lda, const float *B, const INDEXTYPE ldb, 
      const float beta, float *C, const INDEXTYPE ldc);
void sgfusedMM_K32_spmm_b0_csr (const char transa, const INDEXTYPE m, const INDEXTYPE n, 
      const INDEXTYPE k,const float alpha, const INDEXTYPE nnz, 
      const INDEXTYPE rows, const INDEXTYPE cols, const float *val, const INDEXTYPE *indx, 
      const INDEXTYPE *pntrb, const INDEXTYPE *pntre, const float *A, 
      const INDEXTYPE lda, const float *B, const INDEXTYPE ldb, 
      const float beta, float *C, const INDEXTYPE ldc);
/*
 * keep a global array of function pointer to select the correct one 
 */
   kern_sgfusedMM_spmm_b0_t sgenkernels_spmm_b0[4] = 
   {
      sgfusedMM_K8_spmm_b0_csr,      /*  0 */
      sgfusedMM_K16_spmm_b0_csr,      /*  1 */
      sgfusedMM_K24_spmm_b0_csr,      /*  2 */
      sgfusedMM_K32_spmm_b0_csr      /*  3 */
   };
/*
 * Kernels for beta, b1
 */

typedef void (*kern_sgfusedMM_spmm_b1_t) ( const char transa, const INDEXTYPE m, 
      const INDEXTYPE n, const INDEXTYPE k,const float alpha, const INDEXTYPE nnz, 
      const INDEXTYPE rows, const INDEXTYPE cols, const float *val, const INDEXTYPE *indx, 
      const INDEXTYPE *pntrb, const INDEXTYPE *pntre, const float *A, 
      const INDEXTYPE lda, const float *B, const INDEXTYPE ldb, 
      const float beta, float *C, const INDEXTYPE ldc);
/*
 * FIXME: add a check for VLEN in generator with simd.h 
 * VLEN = 8, MAX DIM(K) = 32
 */
void sgfusedMM_K8_spmm_b1_csr (const char transa, const INDEXTYPE m, const INDEXTYPE n, 
      const INDEXTYPE k,const float alpha, const INDEXTYPE nnz, 
      const INDEXTYPE rows, const INDEXTYPE cols, const float *val, const INDEXTYPE *indx, 
      const INDEXTYPE *pntrb, const INDEXTYPE *pntre, const float *A, 
      const INDEXTYPE lda, const float *B, const INDEXTYPE ldb, 
      const float beta, float *C, const INDEXTYPE ldc);
void sgfusedMM_K16_spmm_b1_csr (const char transa, const INDEXTYPE m, const INDEXTYPE n, 
      const INDEXTYPE k,const float alpha, const INDEXTYPE nnz, 
      const INDEXTYPE rows, const INDEXTYPE cols, const float *val, const INDEXTYPE *indx, 
      const INDEXTYPE *pntrb, const INDEXTYPE *pntre, const float *A, 
      const INDEXTYPE lda, const float *B, const INDEXTYPE ldb, 
      const float beta, float *C, const INDEXTYPE ldc);
void sgfusedMM_K24_spmm_b1_csr (const char transa, const INDEXTYPE m, const INDEXTYPE n, 
      const INDEXTYPE k,const float alpha, const INDEXTYPE nnz, 
      const INDEXTYPE rows, const INDEXTYPE cols, const float *val, const INDEXTYPE *indx, 
      const INDEXTYPE *pntrb, const INDEXTYPE *pntre, const float *A, 
      const INDEXTYPE lda, const float *B, const INDEXTYPE ldb, 
      const float beta, float *C, const INDEXTYPE ldc);
void sgfusedMM_K32_spmm_b1_csr (const char transa, const INDEXTYPE m, const INDEXTYPE n, 
      const INDEXTYPE k,const float alpha, const INDEXTYPE nnz, 
      const INDEXTYPE rows, const INDEXTYPE cols, const float *val, const INDEXTYPE *indx, 
      const INDEXTYPE *pntrb, const INDEXTYPE *pntre, const float *A, 
      const INDEXTYPE lda, const float *B, const INDEXTYPE ldb, 
      const float beta, float *C, const INDEXTYPE ldc);
/*
 * keep a global array of function pointer to select the correct one 
 */
   kern_sgfusedMM_spmm_b1_t sgenkernels_spmm_b1[4] = 
   {
      sgfusedMM_K8_spmm_b1_csr,      /*  0 */
      sgfusedMM_K16_spmm_b1_csr,      /*  1 */
      sgfusedMM_K24_spmm_b1_csr,      /*  2 */
      sgfusedMM_K32_spmm_b1_csr      /*  3 */
   };
#endif
