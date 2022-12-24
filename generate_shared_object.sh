mkdir -p ./shared_object/C

cp ./bin/sOptFusedMM_pt.o ./shared_object/C/sOptFusedMM_pt.o
cp ./bin/sFusedMMtime_fr_pt.o ./shared_object/C/sFusedMMtime_fr_pt.o
cp ./kernels/lib/slibgfusedMM_pt.a ./shared_object/C/slibgfusedmm_pt.a
cd shared_object/C


mkdir -p ./objects
ar x slibgfusedmm_pt.a --output ./objects
cp *.o ./objects
ar rcs libmyfusedmm.a ./objects/*.o
gcc -O3 -march=native -std=c++11 -shared -fPIC -o libmyfusedmm_shared.so -Wl,--whole-archive libmyfusedmm.a -Wl,--no-whole-archive -Wl,-rpath='./' -lm -fopenmp 
