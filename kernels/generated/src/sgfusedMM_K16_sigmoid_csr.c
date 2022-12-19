#include<stdint.h>
#ifdef PTTIME
   #include<omp.h>
#endif
#include<stdio.h>
#include<math.h>
#define SREAL 1
#include"../../simd/simd.h"
#if VLEN != 8
   #error "ARCH VLEN doesn't match with generator's VLEN, see simd.h " 
#endif

/*
 * Register block  A,C and B(innermost loop) will require most registers, works 
 * better on small value of k
 */
#ifndef SOP_INHOUSE 
extern int SOP_UDEF_FUNC(float val, float *out);  
#endif
/* external declaration of misc functions  */
#ifdef BETA0 
void sgfusedMM_K16_sigmoid_b0_csr
#else /* BETA1 version */
void sgfusedMM_K16_sigmoid_b1_csr
#endif
(
   const char tkern,  	   // 's' 't' 'm'
   const INDEXTYPE m,      // rows of dense A matrix 
   const INDEXTYPE n,      // rows of dense B matrix
   const INDEXTYPE k,      // cols of A or dimension. not used since K compile time   
   const float alpha,     // const to scale, not use yet  
   const INDEXTYPE nnz,    // nonzeros of the sparse matrix 
   const INDEXTYPE rows,   // number of rows of the sparse matrix  
   const INDEXTYPE cols,   // number of columns of the sparse matrix 
   const float *val,       // value of  the sparse matrix 
   const INDEXTYPE *indx,  // colids -> column indices of sparse matrix 
   const INDEXTYPE *pntrb, // starting index for rowptr of csr of sparse matrix
   const INDEXTYPE *pntre, // ending index for rowptr of csr of sparse matrix 
   const float *a,        // Dense A matrix
   const INDEXTYPE lda,    // leading dimension of a (col size since row-major)  
   const float *b,        // Dense B matrix
   const INDEXTYPE ldb,    // leading dimension of b (col size since row-major)  
   const float beta,      // beta value, compile time not used  
   float *c,              // Dense matrix c
   const INDEXTYPE ldc     // leading dimension size of c (col size since roa-major) 
)
{
#ifdef SOP_INHOUSE
   const float sm_bound = 5.0;
   const int sm_table_size = 2048;
   const float sm_resolution = sm_table_size/(2.0 * sm_bound);

   float *sm_table = (float*)malloc(sizeof(float)*sm_table_size);
   if (!sm_table)
   {
      fprintf(stderr, 
      "Not enough memory to allocate SM TABLE in kernel, SM_TABLE_SIZE = %d!!!\n", 
              sm_table_size);
      exit(0);
   }
   { // init_sm_table 
      for(INDEXTYPE i = 0; i < sm_table_size; i++)
      {
         float x;
         x = 2.0 * sm_bound * i / sm_table_size - sm_bound;
         sm_table[i] = 1.0 / (1 + exp(-x));
      }
   }
#endif
#if defined(PTTIME) && defined(LDB)
   omp_set_num_threads(NTHREADS);
   #pragma omp parallel
   {
      INDEXTYPE RowPerThd, tt;
      INDEXTYPE i, rowb, rowe;
      INDEXTYPE Mnnz = 0; /* non-zero count in M rows  */
      INDEXTYPE deg, cumRow, curRow;
      INDEXTYPE id = omp_get_thread_num();
      INDEXTYPE nthreads = omp_get_num_threads(); 
      
      for (i=0; i < m; i++)
         Mnnz += (pntre[i] - pntrb[i]); 
      RowPerThd = Mnnz / nthreads; 
      
      curRow = cumRow = 0; 
      tt = 1; 
      rowe = -1;  /* init */
      /* set rowstart for 1st thread */ 
      if (id == 0) 
         rowb = 0;
      for (i=0; i < m; i++)
      {
         deg = pntre[i] - pntrb[i]; 
         cumRow += deg;
         curRow += deg;
         if (curRow > RowPerThd)
         {
            if (tt == id)
               rowb = i; 
            else if (tt == id+1)
               rowe = i; 
            curRow = 0;
            RowPerThd = (Mnnz - cumRow) / (nthreads - tt);
            tt += 1; 
         }
      }
      if (tt == id+1)
         rowe = m; 

      for (i=rowb; i < rowe; i++)
#else /* not LBD or not PTTIME */
   #ifdef PTTIME
      #ifdef NTHREADS
      omp_set_num_threads(NTHREADS);
      #endif
      #ifdef DYNAMIC 
         #pragma omp parallel for schedule(dynamic)
      #else
         #pragma omp parallel for schedule(static)
      #endif
   #endif
   for (INDEXTYPE i = 0; i < m; i++)
#endif
   {
      register VTYPE Va0, Vc0, Va1, Vc1;
      INDEXTYPE iindex = i * 16; 
      const float *Ai = a + iindex; 
      float *Ci = c + iindex; 
      VTYPE VMAXBOUND, VMINBOUND; 
#ifdef SOP_INHOUSE
      BCL_vset1(VMAXBOUND, sm_bound); 
      BCL_vset1(VMINBOUND, -sm_bound); 
#endif
#ifdef BETA0
/*
 * NO need to load C, just zerod Vector register    
 */
      BCL_vzero(Vc0); 
      BCL_vzero(Vc1); 

#else /* beta1 */
      // load Vc 
      BCL_vldu(Vc0, Ci+VLEN*0); 
      BCL_vldu(Vc1, Ci+VLEN*1); 
#endif
      // load Va 
      BCL_vldu(Va0, Ai+VLEN*0); 
      BCL_vldu(Va1, Ai+VLEN*1); 

      for (INDEXTYPE j = pntrb[i]; j < pntre[i]; j++)
      {
         VTYPE Vb0, Vb1;
         VTYPE Vd0, Vd1; 
         float d1;
         VTYPE Vatt0, Vatt1;
         float attrc = 0;
         INDEXTYPE colidj = indx[j];
         INDEXTYPE jindex = colidj*16;
         const float *Bj = b + jindex; 
         // load Vxj 
         BCL_vldu(Vb0, Bj+VLEN*0); 
         BCL_vldu(Vb1, Bj+VLEN*1); 
      // init Vatt  
         BCL_vzero(Vatt0);
         BCL_vzero(Vatt1);

      // vmac 
         BCL_vmac(Vatt0, Va0, Vb0);
         BCL_vmac(Vatt1, Va1, Vb1);
         // binary tree reduction 
         BCL_vadd(Vatt0, Vatt0, Vatt1);

         BCL_vrsum1(attrc, Vatt0);
#ifdef SOP_INHOUSE
         /* Calculating Sigmoid value */
         { // fast_SM 
            //d1 = fast_SM(attrc, sm_table);
            if (attrc > sm_bound) d1 = 1.0;
            else if (attrc < -sm_bound) d1 = 0.0;
            else d1 = sm_table[(INDEXTYPE) ((attrc+sm_bound)*sm_resolution)];
         }
         //d1 = STEP * degi * (1.0 - d1);
         d1 = (1.0 - d1);
#else
         SOP_UDEF_FUNC(attrc, &d1);
#endif
         BCL_vset1(Vd1, d1);
         // vmac 
         BCL_vmac(Vc0, Vd1, Vb0);
         BCL_vmac(Vc1, Vd1, Vb1);
      }
      BCL_vstu(Ci + VLEN*0, Vc0); 
      BCL_vstu(Ci + VLEN*1, Vc1); 
   }
#if defined(PTTIME) && defined(LDB)
   }
#endif
#ifdef SOP_INHOUSE
   free(sm_table);
#endif
}
