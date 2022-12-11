#include<stdint.h>
#ifdef PTTIME
   #include<omp.h>
#endif
#define SREAL 1
#include"../../simd/simd.h"
#if VLEN != 8
   #error "ARCH VLEN doesn't match with generator's VLEN, see simd.h " 
#endif

/*
 * Register block  A,C and B(innermost loop) will require most registers, works 
 * better on small value of k
 */
extern int SOP_UDEF_FUNC(float val, float *out);  
/*extern INDEXTYPE MAXBOUND ;*/
#ifdef BETA0
void sgfusedMM_K32_tdist_b0_csr
#else /* BETA1 version */
void sgfusedMM_K32_tdist_b1_csr
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
#if 0
   const float maxbound = 5.0;
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
      register VTYPE Va0, Vc0, Va1, Vc1, Va2, Vc2, Va3, Vc3;
      INDEXTYPE iindex = i * 32; 
      const float *Ai = a + iindex; 
      float *Ci = c + iindex; 
      VTYPE VMAXBOUND, VMINBOUND; 
#if 0
      BCL_vset1(VMAXBOUND, maxbound); 
      BCL_vset1(VMINBOUND, -maxbound); 
#endif
#ifdef BETA0
/*
 * NO need to load C, just zerod Vector register    
 */
      BCL_vzero(Vc0); 
      BCL_vzero(Vc1); 
      BCL_vzero(Vc2); 
      BCL_vzero(Vc3); 

#else /* beta1 */
      // load Vc 
      BCL_vldu(Vc0, Ci+VLEN*0); 
      BCL_vldu(Vc1, Ci+VLEN*1); 
      BCL_vldu(Vc2, Ci+VLEN*2); 
      BCL_vldu(Vc3, Ci+VLEN*3); 
#endif
      // load Va 
      BCL_vldu(Va0, Ai+VLEN*0); 
      BCL_vldu(Va1, Ai+VLEN*1); 
      BCL_vldu(Va2, Ai+VLEN*2); 
      BCL_vldu(Va3, Ai+VLEN*3); 

      for (INDEXTYPE j = pntrb[i]; j < pntre[i]; j++)
      {
         VTYPE Vb0, Vb1, Vb2, Vb3;
         VTYPE Vd0, Vd1, Vd2, Vd3, Vt;
         VTYPE Vatt0, Vatt1, Vatt2, Vatt3;
         float attrc = 0;
         INDEXTYPE colidj = indx[j];
         INDEXTYPE jindex = colidj*32;
         const float *Bj = b + jindex; 
         // load Vxj 
         BCL_vldu(Vb0, Bj+VLEN*0); 
         BCL_vldu(Vb1, Bj+VLEN*1); 
         BCL_vldu(Vb2, Bj+VLEN*2); 
         BCL_vldu(Vb3, Bj+VLEN*3); 
      // init Vatt  
         BCL_vzero(Vatt0);
         BCL_vzero(Vatt1);
         BCL_vzero(Vatt2);
         BCL_vzero(Vatt3);

      // vsub
         BCL_vsub(Vd0, Va0, Vb0);
         BCL_vsub(Vd1, Va1, Vb1);
         BCL_vsub(Vd2, Va2, Vb2);
         BCL_vsub(Vd3, Va3, Vb3);
      // vmac 
         BCL_vmac(Vatt0, Vd0, Vd0);
         BCL_vmac(Vatt1, Vd1, Vd1);
         BCL_vmac(Vatt2, Vd2, Vd2);
         BCL_vmac(Vatt3, Vd3, Vd3);
         // binary tree reduction 
         BCL_vadd(Vatt0, Vatt0, Vatt1);
         BCL_vadd(Vatt2, Vatt2, Vatt3);
         BCL_vadd(Vatt0, Vatt0, Vatt2);

         BCL_vrsum1(attrc, Vatt0);
#if 0
         BCL_vset1(Vatt0, attrc); // a = a
         BCL_vset1(Vt, 1.0f); // t = 1.0
         BCL_vadd(Vatt0, Vatt0, Vt); // a = 1.0 + a
         BCL_vrcp(Vatt0, Vatt0); // a = 1/a
         BCL_vset1(Vt, -2.0f); // t = -2
         BCL_vmul(Vatt0, Vatt0, Vt); // a = -2 * a
         BCL_vmax(Vatt0, Vatt0, VMINBOUND);
         BCL_vmin(Vatt0, Vatt0, VMAXBOUND); 
#else
         SOP_UDEF_FUNC(attrc, &attrc);
         BCL_vset1(Vatt0, attrc); 
#endif
      // vmul 
         BCL_vmul(Vd0, Vatt0, Vd0);
         BCL_vmul(Vd1, Vatt0, Vd1);
         BCL_vmul(Vd2, Vatt0, Vd2);
         BCL_vmul(Vd3, Vatt0, Vd3);
      // vadd 
         BCL_vadd(Vc0, Vc0, Vd0);
         BCL_vadd(Vc1, Vc1, Vd1);
         BCL_vadd(Vc2, Vc2, Vd2);
         BCL_vadd(Vc3, Vc3, Vd3);
      }
      BCL_vstu(Ci + VLEN*0, Vc0); 
      BCL_vstu(Ci + VLEN*1, Vc1); 
      BCL_vstu(Ci + VLEN*2, Vc2); 
      BCL_vstu(Ci + VLEN*3, Vc3); 
   }
#if defined(PTTIME) && defined(LDB)
   }
#endif
}
