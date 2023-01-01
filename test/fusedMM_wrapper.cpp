#include <stdio.h>
#include <stdlib.h>


#define VALUETYPE float

#define ATL_Cachelen 64
   #define ATL_MulByCachelen(N_) ( (N_) << 6 )
   #define ATL_DivByCachelen(N_) ( (N_) >> 6 )

#define ATL_AlignPtr(vp) (void*) \
        ATL_MulByCachelen(ATL_DivByCachelen((((size_t)(vp))+ATL_Cachelen-1)))



int main()
{
      int lda, ldb, ldc;
      lda = ldb = ldc = 2; // both row major, K multiple of VLEN 
      int M = 3;
      int N = 2;
      int szA, szB, szC, szAligned;
      float * pa, *pb, *pc, *pc0, *a, *b, *c, *c0;

      /*
      * NOTE: not sure about system's VLEN from this user code. So, make it cacheline
      * size aligned ....
      */
      szAligned = ATL_Cachelen / sizeof(VALUETYPE);
      szA = ((M*ldb+szAligned-1)/szAligned)*szAligned;  // szB in element
      szB = ((N*ldb+szAligned-1)/szAligned)*szAligned;  // szB in element
      szC = ((M*ldc+szAligned-1)/szAligned)*szAligned;  // szC in element 

      int lenA, lenB, lenC;
      
      lenA = szA*sizeof(VALUETYPE)+2*ATL_Cachelen;
      pa = (VALUETYPE*)malloc(lenA);
      
      a = (VALUETYPE*) ATL_AlignPtr(pa);

      lenB = szB*sizeof(VALUETYPE)+2*ATL_Cachelen;
      pb = (VALUETYPE*)malloc(lenB);
      b = (VALUETYPE*) ATL_AlignPtr(pb);


      pc0 = (VALUETYPE*)malloc(szC*sizeof(VALUETYPE)+2*ATL_Cachelen);
      
      c0 = (VALUETYPE*) ATL_AlignPtr(pc0); 

      lenC = szC*sizeof(VALUETYPE)+2*ATL_Cachelen;
      pc = (VALUETYPE*)malloc(lenC);
      
      c = (VALUETYPE*) ATL_AlignPtr(pc); 
      printf("%p\t%p\t%p\n", pa, pb, pc);
      printf("%p-%p\t%p-%p\t%p-%p\n", a, a+lenA, b, b+lenB, c, c+lenC);
      printf("%d\t%d\t%d\n", lenA, lenB, lenC);
}
