# snapper_core
C functions for core machinery of Snapper

You will need FFTW3 and Cblas in order to compile; using below command

gcc -c -fPIC snapper_core.c -lclblas -lblas -lm -lfftw3 -o snapper_core.o
