FC=ifort
FC=gfortran

all:
	$(FC) -fopenmp -O3 module.f90 md.f90

clean:
	rm -f a.out properties epotke meansqdisp gr nr xvconf *.mod *.xyz *.out
