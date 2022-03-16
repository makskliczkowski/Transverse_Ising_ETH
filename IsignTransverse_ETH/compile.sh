#!/bin/bash

module purge
#unset MKL_SERIAL
#export CPATH="/usr/include/hdf5/serial"
#export LD_LIBRARY_PATH="/home/rswietek/LIBRARIES_CPP/armadillo-10.8.2/lib64" #:"${MKLROOT}/lib/intel64"


module load OpenBLAS/0.3.18-GCC-11.2.0
module load intel/2022.00

icc main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising.o\
 -I../LIBRARIES_CPP/armadillo-10.8.2/include/ -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lblas -lopenblas\
 -no-multibyte-chars -pthread -lstdc++fs -qopenmp -qmkl=parallel -std=c++17 -O3


#g++ main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp\
# user_interface.cpp -o Ising.o -pthread -I../LIBRARIES_CPP/armadillo-10.8.2/include/ -larmadillo\
# -fopenmp -lpthread -lm -lstdc++fs -fcx-fortran-rules\
# -fomit-frame-pointer -std=c++20 -O3

#g++ -std=c++2a main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising.o\
# -DARMA_DONT_USE_WRAPPER -llapack -lopenblas -fopenmp -lpthread -lm -lstdc++fs -fomit-frame-pointer -O3 >& log_compile.txt

#clang++ -std=c++20 -I../LIBRARIES_CPP/armadillo-10.8.2/include/ main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising.o -pthread\
#   -L${MKLROOT}/lib/intel64 -fopenmp -lstdc++fs -llapack -fomit-frame-pointer -lopenblas -Ofast

#icpc main.cpp IsingModel.cpp IsingModel_disorder.cpp IsingModel_sym.cpp tools.cpp user_interface.cpp -o Ising.o\
# -I../LIBRARIES_CPP/armadillo-10.8.2/include/ -L${MKLROOT}/lib/intel64 -laramdillo\
# -lmkl_rt -lpthread -lm -ldl -lstdc++fs -qopenmp -mkl=parallel -std=c++20 -O3
