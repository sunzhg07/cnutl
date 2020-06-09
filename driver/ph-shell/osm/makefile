LIBS = -lblas -llapack -lm
FC=mpif90  -fbounds-check -Wall
FV = ${FC} -O3
FVD=${FC} -Og -g -fbounds-check  -Wall -Wextra -ftrapv -fcheck=all  -Wuninitialized
all = osm_module.o osm_io.o osm_int.o osm_lanczos.o osm_basis.o osm_matrix.o osm_main.o

sm.exe :	${all}
	   ${FV} -o sm.exe ${all} ${LIBS}

%.o: %.f90
	${FV} -c $< ${LIBS}

clean:
	rm *.mod *.exe *.o *__genmod.f90
		
