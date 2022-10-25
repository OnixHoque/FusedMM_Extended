#include <stdio.h>
#include "fusedMMtime.cpp"
int main(int narg, char **argv)
{
   INDEXTYPE M, K;
   VALUETYPE alpha, beta;
   int option, csKB, nrep, isTest, skHd, nrblk;
   char tkern;
   string inputfile; 
   GetFlags(narg, argv, inputfile, option, M, K, csKB, 
   	nrep, isTest, skHd, alpha, beta, tkern);
   printf("%s\n", inputfile);
   // GetSpeedup(inputfile, option, M, K, csKB, nrep, isTest, skHd, alpha, beta, 
         // tkern);
   return 0;
}
