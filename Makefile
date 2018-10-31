OP = -Ofast -march=native
#OP = -O3 -march=native

WR = -Wall -Wextra

all:	y15_subroutine.o y07_subroutine.o y15c_subroutine.o dgchbv.o lapack.o blas.o simpleExample simplePerformance

dgchbv.o:	dgchbv.f
	gfortran -c $(OP) dgchbv.f

lapack.o:	lapack.f
	gfortran -c $(OP) lapack.f

blas.o:	blas.f
	gfortran -c $(OP) -std=legacy blas.f

y15_subroutine.o:	y15_subroutine.f90
	gfortran -c $(OP) y15_subroutine.f90

y07_subroutine.o:	y07_subroutine.f90
	gfortran -c $(OP) y07_subroutine.f90

y15c_subroutine.o:	y15c_subroutine.h y15c_subroutine.cc
	g++ -c $(OP) $(WR) -std=c++17 y15c_subroutine.cc

simpleExample:	simpleExample.cc y15_subroutine.o y07_subroutine.o y15c_subroutine.o dgchbv.o lapack.o blas.o
	g++ $(OP) $(WR) -std=c++17 simpleExample.cc -osimpleExample y15_subroutine.o y07_subroutine.o y15c_subroutine.o -lgfortran -L/usr/local/lib64/ -Wl,-rpath -Wl,/usr/local/lib64/ dgchbv.o lapack.o blas.o
	strip simpleExample

simplePerformance: simplePerformance.cc y15_subroutine.o y15c_subroutine.o dgchbv.o lapack.o blas.o
	g++ $(OP) $(WR) -std=c++17 simplePerformance.cc -osimplePerformance y15_subroutine.o y15c_subroutine.o -lgfortran -L/usr/local/lib64/ -Wl,-rpath -Wl,/usr/local/lib64/ dgchbv.o lapack.o blas.o -Wno-unused-variable
	strip simpleExample

clean:
	rm y15_subroutine.o y07_subroutine.o y15c_subroutine.o dgchbv.o lapack.o blas.o simpleExample simplePerformance
