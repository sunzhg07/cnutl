#!/ usr / bin / env python3
#- * - coding : utf - 8 - * -

import numpy as np
import time
import sys
import coupling
from get_ME_gsm import *
from tensor_decompose import *

start = time.process_time()

K = 2
file_in_prefix = "../Input/Input_gsm_formate/"
file_sp_name  = "vshell_bare_2_2N-22_hw14_add_P_C_2p_1p_ct_sp.dat"
file_int_name = "vshell_bare_2_2N-22_hw14_add_P_C_2p_1p_ct.dat"

file_in_sp  = file_in_prefix + file_sp_name
file_in_int = file_in_prefix + file_int_name

file_out_prefix = "../Output/"
file_out_name = "vshell_bare_2_2N-22_hw14_add_P_C_2p_1p_ct"
file_tensor = str(K)
file_out_name = file_out_name + "_K_" + file_tensor
file_out_sp  = file_out_prefix + file_out_name + "_sp.dat"
file_out_int = file_out_prefix + file_out_name + ".dat"


Coup_cpp = Coupling_F_cpp()

#sys.exit()
GM = Get_ME_gsm()
GM.read_sp(file_in_sp)
GM.read_int(file_in_int)
GM.print_sp()
#sys.exit()
#GM.print_me()

TD = Tensor_Decompose(GM, Coup_cpp)
#
#
# #sys.exit()
# #val = TD.cal_print(2,2,1,3,1,0,K)
#
# #print("Val_test = ",val)
print('Start to build')
TD.build(K)
# #sys.exit()
print('Start to print')
# #TD.print_me()
TD.print_me_f(file_out_int)

elapsed = (time.process_time() - start)
print("Time used:",elapsed)
