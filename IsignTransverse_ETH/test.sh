

module load HDF5/1.12.0-gompi-2021a
module load OpenBLAS/0.3.18-GCC-11.2.0
module load intel/2022.00

icc main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising_test.o\
 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lblas -lopenblas -lhdf5\
 -no-multibyte-chars -pthread -lstdc++fs -qopenmp -qmkl=parallel -std=c++17 -O3

module load OpenBLAS/0.3.18-GCC-11.2.0
module load HDF5/1.12.0-gompi-2021a
module load intel/2022.00

./Ising_test.o -f input.txt >& run_test.log