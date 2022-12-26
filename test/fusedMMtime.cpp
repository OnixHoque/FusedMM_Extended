/*
 * This file provides examples to demonstrate how to use/call our FusedMM kernel 
 * library. It also consists of tester and timer for different applications 
 * which use both FusedMM and trusted implementations (some directly copied from 
 * the application repository) and compare the results. 
 */
extern "C" {

#include <cstdio>
#include <cstdint>
#include <cstdlib>
// #include <random>
#include <cmath>

#include <cassert>
#include <omp.h>
/*
 * Header files for I/O utilities to convert the mtx dataset into CSR or CSC 
 * format. This implementation is taken from previous HipGraph projects managed
 * by Dr. Azad.
 */
// #include "include/CSC.h"
// #include "include/CSR.h"
// #include "include/commonutility.h"
// #include "include/utility.h"

/*
 * NOTE: please select data types for both VALUETYPPE and INDEXTYPE in Makefile
 * using following make var 
 *    pre=[s,d]
 *    ityp=[int64_t,int32_t]
 */

#ifdef DREAL
   #define VALUETYPE double
#else
   #define VALUETYPE float
#endif

/*
 * Added header file for general fusedMM 
 */
#include "../fusedMM.h"

/*
 * Check whether the system supports the desire int data type  
 */   
#ifdef INT64
   #ifndef INT64_MAX 
      #error "64bit integer not supported in this architecture!!!"
   #endif
#endif
#ifdef INT32
   #ifndef INT32_MAX 
      #error "32bit integer not supported in this architecture!!!"
   #endif
#endif
/*
 * some misc definition for timer, this definitions are taken from ATLAS 
 * (math-atlas) project 
 */
#define ATL_Cachelen 64
   #define ATL_MulByCachelen(N_) ( (N_) << 6 )
   #define ATL_DivByCachelen(N_) ( (N_) >> 6 )

#define ATL_AlignPtr(vp) (void*) \
        ATL_MulByCachelen(ATL_DivByCachelen((((size_t)(vp))+ATL_Cachelen-1)))


#define SM_TABLE_SIZE 2048
#define SM_BOUND 5.0
#define SM_RESOLUTION SM_TABLE_SIZE/(2.0 * SM_BOUND)
/*=============================================================================
 * Test kernels: 
 *    In our new design, We will always call fusedMM, analyzing patterns it may 
 *    call optimized kernel from there
 *
 *============================================================================*/
/* 
 * Accessory funcitons to compute sigmoid 
 */
VALUETYPE *SM_TABLE;
inline VALUETYPE uscale_SM(VALUETYPE val)
{
   VALUETYPE sval;
   sval = (val > SM_BOUND) ? SM_BOUND : val;
   sval = (val < -SM_BOUND) ? -SM_BOUND : val;
   return(sval); 
}
void uinit_SM_TABLE()
{
   VALUETYPE x;
   SM_TABLE = (VALUETYPE*)malloc(SM_TABLE_SIZE*sizeof(VALUETYPE));
   assert(SM_TABLE);
   for(INDEXTYPE i = 0; i < SM_TABLE_SIZE; i++)
   {
      x = 2.0 * SM_BOUND * i / SM_TABLE_SIZE - SM_BOUND;
      SM_TABLE[i] = 1.0 / (1 + exp(-x));
   }
}

VALUETYPE ufast_SM(VALUETYPE v)
{
   if (v > SM_BOUND) return 1.0;
   else if (v < -SM_BOUND) return 0.0;
   return SM_TABLE[(INDEXTYPE)((v + SM_BOUND) * SM_RESOLUTION)];
}

VALUETYPE tscale(VALUETYPE v)
{
   if(v > SM_BOUND) return SM_BOUND;
   else if(v < -SM_BOUND) return -SM_BOUND;
   return v;
}

/*
 * NOTE: The implementation of User defined functions differ from different 
 * models. We need to enable/disable it compile time!!!!
 */


extern "C" int SOP_UDEF_FUNC(VALUETYPE val, VALUETYPE *out);
#ifdef SIGMOID_UDEF 
// USER DEFINED FUNCTION for SOP with Sigmoid calc 
int  SOP_UDEF_FUNC(VALUETYPE val, VALUETYPE *out)
{
   *out = 1.0 - ufast_SM(val);
   return FUSEDMM_SUCCESS_RETURN;
}
#elif defined(FR_UDEF)
// SOP_UDEF for FR model
int SOP_UDEF_FUNC(VALUETYPE val, VALUETYPE *out)
{
   *out = 1.0 + 1.0 / val;
   return FUSEDMM_SUCCESS_RETURN;
}
#elif defined(TDIST_UDEF)
// SOP_UDEF for t-distribution  model
int SOP_UDEF_FUNC(VALUETYPE val, VALUETYPE *out)
{  
   *out = tscale(-2.0 / (1.0 + val));
   return FUSEDMM_SUCCESS_RETURN;
}
#elif defined(LL_UDEF)
// SOP_UDEF for LL model
int SOP_UDEF_FUNC(VALUETYPE val, VALUETYPE *out)
{
   *out = log2(1 + sqrt(val));;
   return FUSEDMM_SUCCESS_RETURN;
}
#elif defined(FA_UDEF)
// SOP_UDEF for FA model
int SOP_UDEF_FUNC(VALUETYPE val, VALUETYPE *out)
{
   *out = sqrt(val) + 1.0 / val;;
   return FUSEDMM_SUCCESS_RETURN;
} 
#else 
/*
 * NOTE: other kernels don't use SOP funciton (NOOP or COPY)
 * However, since we enable SOP_UDEF_IMPL in fusedMM.h, we need a dummy func.
 * Normally, users should disable the macro if they don't want to provide any 
 * implementation. We are using this dummy since we use same source for all 
 * the different executables.
 */
int SOP_UDEF_FUNC(VALUETYPE val, VALUETYPE *out)
{
   *out = val;;
   return FUSEDMM_SUCCESS_RETURN;
} 
#endif
#if 0
/*
 * NOTE:  Example of a user-defined ROP and VSC funcitons for tdistribution 
 * application. We use imsg to set the computation    
 */
// tdist 
int ROP_UDEF_FUNC(INDEXTYPE lhs_dim, const VALUETYPE *lhs, INDEXTYPE rhs_dim,
      const VALUETYPE *rhs, VALUETYPE &out)
{
   out = 0.0;
   for (INDEXTYPE i = 0; i < rhs_dim; i += 1)
   {
      out += rhs[i] * rhs[i];
   }  
   return FUSEDMM_SUCCESS_RETURN;
}
// tdist 
int VSC_UDEF_FUNC(INDEXTYPE rhs_dim, const VALUETYPE *rhs, VALUETYPE scal,
      INDEXTYPE out_dim, VALUETYPE *out)
{
   for (INDEXTYPE i = 0; i < rhs_dim; i += 1)
   {
      out[i] = scale(scal * rhs[i]);
   }
   return FUSEDMM_SUCCESS_RETURN;
}
#endif

   
void mytest_csr
(
   const char tkern,       // kernel variations
   const INDEXTYPE m,      // rows of A 
   const INDEXTYPE n,      // rows of B
   const INDEXTYPE k,      // dimension: col of A and B
   const VALUETYPE alpha,  // not used yet  
   const INDEXTYPE nnz,    // nonzeros  
   const INDEXTYPE rows,   // number of rows for sparse matrix 
   const INDEXTYPE cols,   // number of columns for sparse matrix 
   const VALUETYPE *val,   // NNZ value  
   const INDEXTYPE *indx,  // colids -> column indices 
   const INDEXTYPE *pntrb, // starting index for rowptr
   const INDEXTYPE *pntre, // ending index for rowptr
   const VALUETYPE *a,     // Dense A (X) matrix
   const INDEXTYPE lda,    // leading dimension of A (col size since A row-major)  
   const VALUETYPE *b,     // Dense B matrix
   const INDEXTYPE ldb,    // leading dimension of B (col size since B row-major)  
   const VALUETYPE beta,   // beta value 
   VALUETYPE *c,           // Dense matrix c
   const INDEXTYPE ldc     // leading dimension of c (col size since C row-major) 
)
{
   int32_t imsg; 
   switch(tkern)
   {
      case 't' : // t-dist 
	 imsg = VOP_SUBR | ROP_NORMR | SOP_UDEF | VSC_MUL | AOP_ADD;
	 fusedMM_csr(imsg, m, n, k, alpha, nnz, rows, cols, val, indx, pntrb,
               pntre, a, lda, b, ldb, beta, c, ldc);
         break;
      case 'f': // fr model 
	 imsg = VOP_SUBR | ROP_NORMR | SOP_UDEF | VSC_MUL | AOP_ADD;
	 fusedMM_csr(imsg, m, n, k, alpha, nnz, rows, cols, val, indx, pntrb,
               pntre, a, lda, b, ldb, beta, c, ldc);
	 break;
      case 's' : // sigmoid
         uinit_SM_TABLE();    // create sigmoid table to use it from SOP_UDEF
         imsg = VOP_COPY_RHS | ROP_DOT | SOP_UDEF | VSC_MUL | AOP_ADD;
         fusedMM_csr(imsg, m, n, k, alpha, nnz, rows, cols, val, indx, pntrb, 
               pntre, a, lda, b, ldb, beta, c, ldc);
         break;
      case 'm' : // spmm
         imsg = VOP_COPY_RHS | ROP_NOOP | SOP_COPY | VSC_MUL | AOP_ADD;
         
         fusedMM_csr(imsg, m, n, k, alpha, nnz, rows, cols, val, indx, pntrb, 
               pntre, a, lda, b, ldb, beta, c, ldc);
         break;
      case 'g' : // gcn 
         imsg = VOP_COPY_RHS | ROP_NOOP | SOP_NOOP | VSC_NOOP | AOP_ADD;
         
         fusedMM_csr(imsg, m, n, k, alpha, nnz, rows, cols, val, indx, pntrb, 
               pntre, a, lda, b, ldb, beta, c, ldc);
         break;
      default:
         printf("unknown trusted kernel\n");
         break;
   }
}



void SpMM
(
   const INDEXTYPE m,      // rows of A 
   const INDEXTYPE n,      // rows of B
   const INDEXTYPE k,      // dimension: col of A and B
   const VALUETYPE alpha,  // not used yet  
   const INDEXTYPE nnz,    // nonzeros  
   const INDEXTYPE rows,   // number of rows for sparse matrix 
   const INDEXTYPE cols,   // number of columns for sparse matrix 
   const VALUETYPE *val,   // NNZ value  
   const INDEXTYPE *indx,  // colids -> column indices 
   const INDEXTYPE *pntrb, // starting index for rowptr
   const INDEXTYPE *pntre, // ending index for rowptr
   const VALUETYPE *a,     // Dense A (X) matrix
   const INDEXTYPE lda,    // leading dimension of A (col size since A row-major)  
   const VALUETYPE *b,     // Dense B matrix
   const INDEXTYPE ldb,    // leading dimension of B (col size since B row-major)  
   const VALUETYPE beta,   // beta value 
   VALUETYPE *c,           // Dense matrix c
   const INDEXTYPE ldc     // leading dimension of c (col size since C row-major) 
)
{
   // printf("=============Before calling:\n");
   // printf("S:\n");
   // for (int i = 0; i<6; i++)
   //    printf("%f\t", val[i]);

   // printf("\n\nB:\n");
   // for (int i = 0; i<6; i++)
   //    printf("%f\t", b[i]);
   
   // printf("\n\nC:\n");
   // for (int i = 0; i<4; i++)
   //    printf("%f\t", c[i]);
   
   // printf("===============================\n\n");

   mytest_csr('m', m, n, k, alpha, nnz, rows, cols, val, indx, pntrb, pntre, a, lda, b, ldb, beta, c, ldc);

   // printf("=============After calling:\n");
   // printf("S:\n");
   // for (int i = 0; i<6; i++)
   //    printf("%f\t", val[i]);

   // printf("\n\nB:\n");
   // for (int i = 0; i<6; i++)
   //    printf("%f\t", b[i]);
   
   // printf("\n\nC:\n");
   // for (int i = 0; i<4; i++)
   //    printf("%f\t", c[i]);
   
   // printf("===============================\n\n");
}


void performDummySpMM()
{
   printf("Performing Dummy SpMM...\n");
   INDEXTYPE M = 2, N = 3, K = 2;
   VALUETYPE alpha = 1, beta = 0;
   int option = 1, csKB = 25344, nrep = 1, isTest = 0, skHd = 0;
   char tkern = 'm';
   int nerr, szAligned; 
   size_t i, j, szA, szB, szC, lda, ldc, ldb; 
   VALUETYPE *pb, *b, *pc0, *c0, *pc, *c, *pa, *a, *values;
   lda = ldb = ldc = K; // both row major, K multiple of VLEN 


/*
 * NOTE: not sure about system's VLEN from this user code. So, make it cacheline
 * size aligned ....
 */
   szAligned = ATL_Cachelen / sizeof(VALUETYPE);
   szA = ((M*ldb+szAligned-1)/szAligned)*szAligned;  // szB in element
   szB = ((N*ldb+szAligned-1)/szAligned)*szAligned;  // szB in element
   szC = ((M*ldc+szAligned-1)/szAligned)*szAligned;  // szC in element 
   
   pa = (VALUETYPE*)calloc(szA,sizeof(VALUETYPE)+2*ATL_Cachelen);
   assert(pa);
   a = (VALUETYPE*) ATL_AlignPtr(pa);
   
   pb = (VALUETYPE*)calloc(szB,sizeof(VALUETYPE)+2*ATL_Cachelen);
   assert(pb);
   b = (VALUETYPE*) ATL_AlignPtr(pb);

   pc0 = (VALUETYPE*)calloc(szC,sizeof(VALUETYPE)+2*ATL_Cachelen);
   assert(pc0);
   c0 = (VALUETYPE*) ATL_AlignPtr(pc0); 
      
   pc = (VALUETYPE*)calloc(szC,sizeof(VALUETYPE)+2*ATL_Cachelen);
   assert(pc);
   c = (VALUETYPE*) ATL_AlignPtr(pc); 
   

   b[0] = 10; b[1] = 11; b[2] = 20;
   b[3] = 21; b[4] = 30; b[5] = 31;

   
   // if (M > S.rows) M = S.rows; // M can't be greater than A.rows  
/*
 *    csr may consists all 1 as values... init with random values
 */

   //====Sparse Matrix Declaration====
   INDEXTYPE S_rows = 2;	
	INDEXTYPE S_cols = 3;
	INDEXTYPE S_nnz = 6; // number of nonzeros
    
   INDEXTYPE S_rowptr[] = {0, 3, 6};
   INDEXTYPE S_colids[] = {0, 1, 2, 0, 1, 2};
   VALUETYPE S_values[] = {1, 2, 3, 4, 5, 6};

   // values = (VALUETYPE*)calloc(S.nnz*sizeof(VALUETYPE));

   // assert(values);
   
   
   // printf("\nNNZ in S: $d %d %d\n", szA, szB, S.nnz);

   // fprintf(stdout, "Applying test kernel\n");
   SpMM(M, N, K, alpha, S_nnz, S_rows, S_cols, S_values, S_colids, S_rowptr, S_rowptr+1, a, lda, b, ldb, beta, c, ldc);   
   printf("RESULT Matrix C: (Dimsension %d)\n", szC);
   for (i=0; i < szC; i++)
   {
      printf("%f\t", c[i]);
   }


}


// int main(int narg, char **argv)
// {
//    // performSpMM();
//    // printf("Status: %d\n", ENABLE_OPT_FUSEDMM);
//    performDummySpMM();
//    // GetSpeedup(inputfile, option, M, K, csKB, nrep, isTest, skHd, alpha, beta, tkern);
//    return 0;
// }

}