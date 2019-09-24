XLF = gfortran -Og -g -fbounds-check  -Wall -Wextra -ftrapv -fcheck=all -finit-local-zero


all_objects = j2m_module.o j2m_main.o ccm_library.o

prog_ccm_ex.exe : ${all_objects}
	${XLF} -o j2m.exe ${all_objects} 


%.o: %.f90
	${XLF} -c $<


clean:
	rm *.mod *.exe *.o *__genmod.f90


