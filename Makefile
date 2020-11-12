## Compilateur

FC= gfortran 
## Options de compilation
FFLAGS = -O3 -g -fbacktrace -fbounds-check -fopenmp

## Editeur de liens
src= transport.f90 rusanov.f90 flux.f90 numerics.f90 musc.f90
# Fichiers objets
objs= transport.o rusanov.o flux.o numerics.o musc.o

EXEC = exe

## Programme

all: $(EXE)

transport.o: $(src) rusanov.o flux.o numerics.o musc.o
	$(FC) $(FFLAGS) -c $(src) 

numerics.o: numerics.f90
	$(FC) $(FFLAGS) -c numerics.f90

flux.o: flux.f90 numerics.o
	$(FC) $(FFLAGS) -c flux.f90

rusanov.o: rusanov.f90 numerics.o flux.o
	$(FC) $(FFLAGS) -c rusanov.f90

musc.o: musc.f90
	$(FC) $(FFLAGS) -c musc.f90

exe: $(objs) 
	$(FC) $(FFLAGS)   $(objs) -o  exe	

clean:
	rm *.o exe *.mod *.txt
