rm -fdr ./shared_object
make
mkdir -p ./shared_object

cp ./bin/sOptFusedMM_pt.o ./shared_object/sOptFusedMM_pt.o
cp ./bin/sFusedMMtime_spmm_pt.o ./shared_object/sFusedMMtime_spmm_pt.o
# cp ./bin/sFusedMMtime_fr_pt.o ./shared_object/sFusedMMtime_fr_pt.o
cp ./kernels/lib/slibgfusedMM_pt.a ./shared_object/slibgfusedmm_pt.a
cd shared_object/

mkdir -p ./objects
mv slibgfusedmm_pt.a ./objects
cp *.o ./objects
cd ./objects
ar x slibgfusedmm_pt.a
cd ..
ar rcs libmyfusedmm.a ./objects/*.o
cp libmyfusedmm.a ../pybind/
gcc -O3 -march=native -std=c++11 -shared -fPIC -o libmyfusedmm_shared.so -Wl,--whole-archive libmyfusedmm.a -Wl,--no-whole-archive -Wl,-rpath='./' -lm -fopenmp 

echo "Shared object generated in ./shared_object folder."
echo "You can also use pybind. Go to './pybind' and run './pybind_make.sh' or 'pip install .'"