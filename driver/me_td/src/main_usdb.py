#!/ usr / bin / env python3
#- * - coding : utf - 8 - * -

import numpy as np
import time
import sys
import coupling
from get_ME_usdb import *
from tensor_decompose import *
from sympy.physics.wigner import clebsch_gordan
from sympy.physics.wigner import wigner_9j
start = time.process_time()

file_usdb = "../Input/usdb.int"
file_out = "../Output/usdb_K2.int"


Coup_cpp = Coupling_F_cpp()

#sys.exit()
GM = Get_ME_usdb_int()
GM.read_data(file_usdb)
GM.print_sp()
#GM.print_me()

TD = Tensor_Decompose(GM, Coup_cpp)

K = 2
#sys.exit()
#val = TD.cal_print(2,2,1,3,1,0,K)

#print("Val_test = ",val)
print('Start to build')
TD.build(K)
#sys.exit()
print('Start to print')
#TD.print_me()
TD.print_me_f(file_out)

elapsed = (time.process_time() - start)
print("Time used:",elapsed)
