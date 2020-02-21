
FC=mpif90 -mkl -fopenmp 
#FCD=mpif90 -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan -lblas -llapack -fopenmp
FCD=mpif90 -Og -g -fbounds-check  -Wall -Wextra -ftrapv -fcheck=all  -Wuninitialized  -lblas -llapack -fopenmp
FCV = mpif90 -mkl=parallel -g -check -heap-arrays 1024 -fltconsistency


XLF= ${FCD}  ${INCLUDE}

prg.exe: ccm_library.o ccm_modules2.o main.o
	${XLF} ccm_library.o ccm_modules2.o main.o -o prg.exe ${LIBS}

ccm_modules2.o: ccm_modules2.f90
	${XLF} -c ccm_modules2.f90
main.o:	main.f90
	${XLF} -c main.f90
ccm_library.o: ccm_library.f90
	${XLF} -c ccm_library.f90

clean:
	rm *.o *.mod *.exe
