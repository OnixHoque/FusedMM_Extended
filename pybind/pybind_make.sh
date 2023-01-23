c++ -O3 -march=native -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) fusedmm_pybind.cpp *.a -lm -fopenmp -o fusedmm$(python3-config --extension-suffix)
