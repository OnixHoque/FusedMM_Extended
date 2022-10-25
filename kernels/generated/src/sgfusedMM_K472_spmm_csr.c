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
#ifdef BETA0 
void sgfusedMM_K472_spmm_b0_csr
#else /* BETA1 version */
void sgfusedMM_K472_spmm_b1_csr
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
      register VTYPE Va0, Vc0, Va1, Vc1, Va2, Vc2, Va3, Vc3, Va4, Vc4, Va5, Vc5,
                     Va6, Vc6, Va7, Vc7, Va8, Vc8, Va9, Vc9, Va10, Vc10, Va11,
                     Vc11, Va12, Vc12, Va13, Vc13, Va14, Vc14, Va15, Vc15, Va16,
                     Vc16, Va17, Vc17, Va18, Vc18, Va19, Vc19, Va20, Vc20, Va21,
                     Vc21, Va22, Vc22, Va23, Vc23, Va24, Vc24, Va25, Vc25, Va26,
                     Vc26, Va27, Vc27, Va28, Vc28, Va29, Vc29, Va30, Vc30, Va31,
                     Vc31, Va32, Vc32, Va33, Vc33, Va34, Vc34, Va35, Vc35, Va36,
                     Vc36, Va37, Vc37, Va38, Vc38, Va39, Vc39, Va40, Vc40, Va41,
                     Vc41, Va42, Vc42, Va43, Vc43, Va44, Vc44, Va45, Vc45, Va46,
                     Vc46, Va47, Vc47, Va48, Vc48, Va49, Vc49, Va50, Vc50, Va51,
                     Vc51, Va52, Vc52, Va53, Vc53, Va54, Vc54, Va55, Vc55, Va56,
                     Vc56, Va57, Vc57, Va58, Vc58;
      INDEXTYPE iindex = i * 472; 
      const float *Ai = a + iindex; 
      float *Ci = c + iindex; 
      VTYPE VMAXBOUND, VMINBOUND; 
#ifdef BETA0
/*
 * NO need to load C, just zerod Vector register    
 */
      BCL_vzero(Vc0); 
      BCL_vzero(Vc1); 
      BCL_vzero(Vc2); 
      BCL_vzero(Vc3); 
      BCL_vzero(Vc4); 
      BCL_vzero(Vc5); 
      BCL_vzero(Vc6); 
      BCL_vzero(Vc7); 
      BCL_vzero(Vc8); 
      BCL_vzero(Vc9); 
      BCL_vzero(Vc10); 
      BCL_vzero(Vc11); 
      BCL_vzero(Vc12); 
      BCL_vzero(Vc13); 
      BCL_vzero(Vc14); 
      BCL_vzero(Vc15); 
      BCL_vzero(Vc16); 
      BCL_vzero(Vc17); 
      BCL_vzero(Vc18); 
      BCL_vzero(Vc19); 
      BCL_vzero(Vc20); 
      BCL_vzero(Vc21); 
      BCL_vzero(Vc22); 
      BCL_vzero(Vc23); 
      BCL_vzero(Vc24); 
      BCL_vzero(Vc25); 
      BCL_vzero(Vc26); 
      BCL_vzero(Vc27); 
      BCL_vzero(Vc28); 
      BCL_vzero(Vc29); 
      BCL_vzero(Vc30); 
      BCL_vzero(Vc31); 
      BCL_vzero(Vc32); 
      BCL_vzero(Vc33); 
      BCL_vzero(Vc34); 
      BCL_vzero(Vc35); 
      BCL_vzero(Vc36); 
      BCL_vzero(Vc37); 
      BCL_vzero(Vc38); 
      BCL_vzero(Vc39); 
      BCL_vzero(Vc40); 
      BCL_vzero(Vc41); 
      BCL_vzero(Vc42); 
      BCL_vzero(Vc43); 
      BCL_vzero(Vc44); 
      BCL_vzero(Vc45); 
      BCL_vzero(Vc46); 
      BCL_vzero(Vc47); 
      BCL_vzero(Vc48); 
      BCL_vzero(Vc49); 
      BCL_vzero(Vc50); 
      BCL_vzero(Vc51); 
      BCL_vzero(Vc52); 
      BCL_vzero(Vc53); 
      BCL_vzero(Vc54); 
      BCL_vzero(Vc55); 
      BCL_vzero(Vc56); 
      BCL_vzero(Vc57); 
      BCL_vzero(Vc58); 

#else /* beta1 */
      // load Vc 
      BCL_vldu(Vc0, Ci+VLEN*0); 
      BCL_vldu(Vc1, Ci+VLEN*1); 
      BCL_vldu(Vc2, Ci+VLEN*2); 
      BCL_vldu(Vc3, Ci+VLEN*3); 
      BCL_vldu(Vc4, Ci+VLEN*4); 
      BCL_vldu(Vc5, Ci+VLEN*5); 
      BCL_vldu(Vc6, Ci+VLEN*6); 
      BCL_vldu(Vc7, Ci+VLEN*7); 
      BCL_vldu(Vc8, Ci+VLEN*8); 
      BCL_vldu(Vc9, Ci+VLEN*9); 
      BCL_vldu(Vc10, Ci+VLEN*10); 
      BCL_vldu(Vc11, Ci+VLEN*11); 
      BCL_vldu(Vc12, Ci+VLEN*12); 
      BCL_vldu(Vc13, Ci+VLEN*13); 
      BCL_vldu(Vc14, Ci+VLEN*14); 
      BCL_vldu(Vc15, Ci+VLEN*15); 
      BCL_vldu(Vc16, Ci+VLEN*16); 
      BCL_vldu(Vc17, Ci+VLEN*17); 
      BCL_vldu(Vc18, Ci+VLEN*18); 
      BCL_vldu(Vc19, Ci+VLEN*19); 
      BCL_vldu(Vc20, Ci+VLEN*20); 
      BCL_vldu(Vc21, Ci+VLEN*21); 
      BCL_vldu(Vc22, Ci+VLEN*22); 
      BCL_vldu(Vc23, Ci+VLEN*23); 
      BCL_vldu(Vc24, Ci+VLEN*24); 
      BCL_vldu(Vc25, Ci+VLEN*25); 
      BCL_vldu(Vc26, Ci+VLEN*26); 
      BCL_vldu(Vc27, Ci+VLEN*27); 
      BCL_vldu(Vc28, Ci+VLEN*28); 
      BCL_vldu(Vc29, Ci+VLEN*29); 
      BCL_vldu(Vc30, Ci+VLEN*30); 
      BCL_vldu(Vc31, Ci+VLEN*31); 
      BCL_vldu(Vc32, Ci+VLEN*32); 
      BCL_vldu(Vc33, Ci+VLEN*33); 
      BCL_vldu(Vc34, Ci+VLEN*34); 
      BCL_vldu(Vc35, Ci+VLEN*35); 
      BCL_vldu(Vc36, Ci+VLEN*36); 
      BCL_vldu(Vc37, Ci+VLEN*37); 
      BCL_vldu(Vc38, Ci+VLEN*38); 
      BCL_vldu(Vc39, Ci+VLEN*39); 
      BCL_vldu(Vc40, Ci+VLEN*40); 
      BCL_vldu(Vc41, Ci+VLEN*41); 
      BCL_vldu(Vc42, Ci+VLEN*42); 
      BCL_vldu(Vc43, Ci+VLEN*43); 
      BCL_vldu(Vc44, Ci+VLEN*44); 
      BCL_vldu(Vc45, Ci+VLEN*45); 
      BCL_vldu(Vc46, Ci+VLEN*46); 
      BCL_vldu(Vc47, Ci+VLEN*47); 
      BCL_vldu(Vc48, Ci+VLEN*48); 
      BCL_vldu(Vc49, Ci+VLEN*49); 
      BCL_vldu(Vc50, Ci+VLEN*50); 
      BCL_vldu(Vc51, Ci+VLEN*51); 
      BCL_vldu(Vc52, Ci+VLEN*52); 
      BCL_vldu(Vc53, Ci+VLEN*53); 
      BCL_vldu(Vc54, Ci+VLEN*54); 
      BCL_vldu(Vc55, Ci+VLEN*55); 
      BCL_vldu(Vc56, Ci+VLEN*56); 
      BCL_vldu(Vc57, Ci+VLEN*57); 
      BCL_vldu(Vc58, Ci+VLEN*58); 
#endif
      for (INDEXTYPE j = pntrb[i]; j < pntre[i]; j++)
      {
         VTYPE Vb0, Vb1, Vb2, Vb3, Vb4, Vb5, Vb6, Vb7, Vb8, Vb9, Vb10, Vb11,
               Vb12, Vb13, Vb14, Vb15, Vb16, Vb17, Vb18, Vb19, Vb20, Vb21, Vb22,
               Vb23, Vb24, Vb25, Vb26, Vb27, Vb28, Vb29, Vb30, Vb31, Vb32, Vb33,
               Vb34, Vb35, Vb36, Vb37, Vb38, Vb39, Vb40, Vb41, Vb42, Vb43, Vb44,
               Vb45, Vb46, Vb47, Vb48, Vb49, Vb50, Vb51, Vb52, Vb53, Vb54, Vb55,
               Vb56, Vb57, Vb58;
         VTYPE Va0; 
         float a0 = val[j];
         INDEXTYPE colidj = indx[j];
         INDEXTYPE jindex = colidj*472;
         const float *Bj = b + jindex; 
         // load Vxj 
         BCL_vldu(Vb0, Bj+VLEN*0); 
         BCL_vldu(Vb1, Bj+VLEN*1); 
         BCL_vldu(Vb2, Bj+VLEN*2); 
         BCL_vldu(Vb3, Bj+VLEN*3); 
         BCL_vldu(Vb4, Bj+VLEN*4); 
         BCL_vldu(Vb5, Bj+VLEN*5); 
         BCL_vldu(Vb6, Bj+VLEN*6); 
         BCL_vldu(Vb7, Bj+VLEN*7); 
         BCL_vldu(Vb8, Bj+VLEN*8); 
         BCL_vldu(Vb9, Bj+VLEN*9); 
         BCL_vldu(Vb10, Bj+VLEN*10); 
         BCL_vldu(Vb11, Bj+VLEN*11); 
         BCL_vldu(Vb12, Bj+VLEN*12); 
         BCL_vldu(Vb13, Bj+VLEN*13); 
         BCL_vldu(Vb14, Bj+VLEN*14); 
         BCL_vldu(Vb15, Bj+VLEN*15); 
         BCL_vldu(Vb16, Bj+VLEN*16); 
         BCL_vldu(Vb17, Bj+VLEN*17); 
         BCL_vldu(Vb18, Bj+VLEN*18); 
         BCL_vldu(Vb19, Bj+VLEN*19); 
         BCL_vldu(Vb20, Bj+VLEN*20); 
         BCL_vldu(Vb21, Bj+VLEN*21); 
         BCL_vldu(Vb22, Bj+VLEN*22); 
         BCL_vldu(Vb23, Bj+VLEN*23); 
         BCL_vldu(Vb24, Bj+VLEN*24); 
         BCL_vldu(Vb25, Bj+VLEN*25); 
         BCL_vldu(Vb26, Bj+VLEN*26); 
         BCL_vldu(Vb27, Bj+VLEN*27); 
         BCL_vldu(Vb28, Bj+VLEN*28); 
         BCL_vldu(Vb29, Bj+VLEN*29); 
         BCL_vldu(Vb30, Bj+VLEN*30); 
         BCL_vldu(Vb31, Bj+VLEN*31); 
         BCL_vldu(Vb32, Bj+VLEN*32); 
         BCL_vldu(Vb33, Bj+VLEN*33); 
         BCL_vldu(Vb34, Bj+VLEN*34); 
         BCL_vldu(Vb35, Bj+VLEN*35); 
         BCL_vldu(Vb36, Bj+VLEN*36); 
         BCL_vldu(Vb37, Bj+VLEN*37); 
         BCL_vldu(Vb38, Bj+VLEN*38); 
         BCL_vldu(Vb39, Bj+VLEN*39); 
         BCL_vldu(Vb40, Bj+VLEN*40); 
         BCL_vldu(Vb41, Bj+VLEN*41); 
         BCL_vldu(Vb42, Bj+VLEN*42); 
         BCL_vldu(Vb43, Bj+VLEN*43); 
         BCL_vldu(Vb44, Bj+VLEN*44); 
         BCL_vldu(Vb45, Bj+VLEN*45); 
         BCL_vldu(Vb46, Bj+VLEN*46); 
         BCL_vldu(Vb47, Bj+VLEN*47); 
         BCL_vldu(Vb48, Bj+VLEN*48); 
         BCL_vldu(Vb49, Bj+VLEN*49); 
         BCL_vldu(Vb50, Bj+VLEN*50); 
         BCL_vldu(Vb51, Bj+VLEN*51); 
         BCL_vldu(Vb52, Bj+VLEN*52); 
         BCL_vldu(Vb53, Bj+VLEN*53); 
         BCL_vldu(Vb54, Bj+VLEN*54); 
         BCL_vldu(Vb55, Bj+VLEN*55); 
         BCL_vldu(Vb56, Bj+VLEN*56); 
         BCL_vldu(Vb57, Bj+VLEN*57); 
         BCL_vldu(Vb58, Bj+VLEN*58); 
         BCL_vset1(Va0, a0);
         // spmm vmac 
         BCL_vmac(Vc0, Va0, Vb0);
         BCL_vmac(Vc1, Va0, Vb1);
         BCL_vmac(Vc2, Va0, Vb2);
         BCL_vmac(Vc3, Va0, Vb3);
         BCL_vmac(Vc4, Va0, Vb4);
         BCL_vmac(Vc5, Va0, Vb5);
         BCL_vmac(Vc6, Va0, Vb6);
         BCL_vmac(Vc7, Va0, Vb7);
         BCL_vmac(Vc8, Va0, Vb8);
         BCL_vmac(Vc9, Va0, Vb9);
         BCL_vmac(Vc10, Va0, Vb10);
         BCL_vmac(Vc11, Va0, Vb11);
         BCL_vmac(Vc12, Va0, Vb12);
         BCL_vmac(Vc13, Va0, Vb13);
         BCL_vmac(Vc14, Va0, Vb14);
         BCL_vmac(Vc15, Va0, Vb15);
         BCL_vmac(Vc16, Va0, Vb16);
         BCL_vmac(Vc17, Va0, Vb17);
         BCL_vmac(Vc18, Va0, Vb18);
         BCL_vmac(Vc19, Va0, Vb19);
         BCL_vmac(Vc20, Va0, Vb20);
         BCL_vmac(Vc21, Va0, Vb21);
         BCL_vmac(Vc22, Va0, Vb22);
         BCL_vmac(Vc23, Va0, Vb23);
         BCL_vmac(Vc24, Va0, Vb24);
         BCL_vmac(Vc25, Va0, Vb25);
         BCL_vmac(Vc26, Va0, Vb26);
         BCL_vmac(Vc27, Va0, Vb27);
         BCL_vmac(Vc28, Va0, Vb28);
         BCL_vmac(Vc29, Va0, Vb29);
         BCL_vmac(Vc30, Va0, Vb30);
         BCL_vmac(Vc31, Va0, Vb31);
         BCL_vmac(Vc32, Va0, Vb32);
         BCL_vmac(Vc33, Va0, Vb33);
         BCL_vmac(Vc34, Va0, Vb34);
         BCL_vmac(Vc35, Va0, Vb35);
         BCL_vmac(Vc36, Va0, Vb36);
         BCL_vmac(Vc37, Va0, Vb37);
         BCL_vmac(Vc38, Va0, Vb38);
         BCL_vmac(Vc39, Va0, Vb39);
         BCL_vmac(Vc40, Va0, Vb40);
         BCL_vmac(Vc41, Va0, Vb41);
         BCL_vmac(Vc42, Va0, Vb42);
         BCL_vmac(Vc43, Va0, Vb43);
         BCL_vmac(Vc44, Va0, Vb44);
         BCL_vmac(Vc45, Va0, Vb45);
         BCL_vmac(Vc46, Va0, Vb46);
         BCL_vmac(Vc47, Va0, Vb47);
         BCL_vmac(Vc48, Va0, Vb48);
         BCL_vmac(Vc49, Va0, Vb49);
         BCL_vmac(Vc50, Va0, Vb50);
         BCL_vmac(Vc51, Va0, Vb51);
         BCL_vmac(Vc52, Va0, Vb52);
         BCL_vmac(Vc53, Va0, Vb53);
         BCL_vmac(Vc54, Va0, Vb54);
         BCL_vmac(Vc55, Va0, Vb55);
         BCL_vmac(Vc56, Va0, Vb56);
         BCL_vmac(Vc57, Va0, Vb57);
         BCL_vmac(Vc58, Va0, Vb58);
      }
      BCL_vstu(Ci + VLEN*0, Vc0); 
      BCL_vstu(Ci + VLEN*1, Vc1); 
      BCL_vstu(Ci + VLEN*2, Vc2); 
      BCL_vstu(Ci + VLEN*3, Vc3); 
      BCL_vstu(Ci + VLEN*4, Vc4); 
      BCL_vstu(Ci + VLEN*5, Vc5); 
      BCL_vstu(Ci + VLEN*6, Vc6); 
      BCL_vstu(Ci + VLEN*7, Vc7); 
      BCL_vstu(Ci + VLEN*8, Vc8); 
      BCL_vstu(Ci + VLEN*9, Vc9); 
      BCL_vstu(Ci + VLEN*10, Vc10); 
      BCL_vstu(Ci + VLEN*11, Vc11); 
      BCL_vstu(Ci + VLEN*12, Vc12); 
      BCL_vstu(Ci + VLEN*13, Vc13); 
      BCL_vstu(Ci + VLEN*14, Vc14); 
      BCL_vstu(Ci + VLEN*15, Vc15); 
      BCL_vstu(Ci + VLEN*16, Vc16); 
      BCL_vstu(Ci + VLEN*17, Vc17); 
      BCL_vstu(Ci + VLEN*18, Vc18); 
      BCL_vstu(Ci + VLEN*19, Vc19); 
      BCL_vstu(Ci + VLEN*20, Vc20); 
      BCL_vstu(Ci + VLEN*21, Vc21); 
      BCL_vstu(Ci + VLEN*22, Vc22); 
      BCL_vstu(Ci + VLEN*23, Vc23); 
      BCL_vstu(Ci + VLEN*24, Vc24); 
      BCL_vstu(Ci + VLEN*25, Vc25); 
      BCL_vstu(Ci + VLEN*26, Vc26); 
      BCL_vstu(Ci + VLEN*27, Vc27); 
      BCL_vstu(Ci + VLEN*28, Vc28); 
      BCL_vstu(Ci + VLEN*29, Vc29); 
      BCL_vstu(Ci + VLEN*30, Vc30); 
      BCL_vstu(Ci + VLEN*31, Vc31); 
      BCL_vstu(Ci + VLEN*32, Vc32); 
      BCL_vstu(Ci + VLEN*33, Vc33); 
      BCL_vstu(Ci + VLEN*34, Vc34); 
      BCL_vstu(Ci + VLEN*35, Vc35); 
      BCL_vstu(Ci + VLEN*36, Vc36); 
      BCL_vstu(Ci + VLEN*37, Vc37); 
      BCL_vstu(Ci + VLEN*38, Vc38); 
      BCL_vstu(Ci + VLEN*39, Vc39); 
      BCL_vstu(Ci + VLEN*40, Vc40); 
      BCL_vstu(Ci + VLEN*41, Vc41); 
      BCL_vstu(Ci + VLEN*42, Vc42); 
      BCL_vstu(Ci + VLEN*43, Vc43); 
      BCL_vstu(Ci + VLEN*44, Vc44); 
      BCL_vstu(Ci + VLEN*45, Vc45); 
      BCL_vstu(Ci + VLEN*46, Vc46); 
      BCL_vstu(Ci + VLEN*47, Vc47); 
      BCL_vstu(Ci + VLEN*48, Vc48); 
      BCL_vstu(Ci + VLEN*49, Vc49); 
      BCL_vstu(Ci + VLEN*50, Vc50); 
      BCL_vstu(Ci + VLEN*51, Vc51); 
      BCL_vstu(Ci + VLEN*52, Vc52); 
      BCL_vstu(Ci + VLEN*53, Vc53); 
      BCL_vstu(Ci + VLEN*54, Vc54); 
      BCL_vstu(Ci + VLEN*55, Vc55); 
      BCL_vstu(Ci + VLEN*56, Vc56); 
      BCL_vstu(Ci + VLEN*57, Vc57); 
      BCL_vstu(Ci + VLEN*58, Vc58); 
   }
#if defined(PTTIME) && defined(LDB)
   }
#endif
}
