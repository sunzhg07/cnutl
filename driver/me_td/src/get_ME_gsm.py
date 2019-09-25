#!/ usr / bin / env python3
#- * - coding : utf - 8 - * -

import numpy as np
import sys

class SP :
    'single particle state'
#sp_num = 0;
    def __init__(self, n, l, j2, t2z, E, index) :
        self.index = index
        self.n = n
        self.l = l
        self.j2 = j2
        self.t2z = t2z
        self.E = E
#Basis_SP.sp_num += 1

class Get_ME_gsm :
    def __init__(self):
        self.sp_state =[]
        self.me_dic = {}

    def phase(self,N):

        if(N%2 == 0):
            return 1
        else:
            return -1
    def read_sp(self, filename):
        f = open(filename)
        while True:
            line = f.readline()
            if not line:      #等价于if line == "":
                break
            A = line.splitlines()[0]
            B = A.split()
            index = int(B[0])
            n=int(B[4])
            l=int(B[5])
            j2=int(B[6])
            t2z=int(B[7])
            E_real = float(B[8])
            E_imag = float(B[9])
            E = complex(E_real, E_imag)
            sp_t = SP(n,l,j2,t2z,E,index)
            self.sp_state.append(sp_t)
            #print(index, n, l, j2)
        self.size_sp = len(self.sp_state)

    def read_int(self,filename):
        f = open(filename)
#-- - -- - get ME -- - -- -
        num = 0
        conf = "null"
        ME = "null"
        while True:
            line = f.readline()
            if not line:      #等价于if line == "":
                print(conf, ME)
                break
            #print(line)
            num += 1
            A = line.splitlines()[0]
            B = A.split()

            conf = B[0] +"-" + B[1] +"-"+ B[2] +"-"+ B[3]+"-"+ B[4]+"-"+ B[5]
            ME_real = float(B[6])
            ME_imag = float(B[7])
            ME = complex(ME_real, ME_imag)
            self.me_dic[conf] = ME
            if(num == 1):
                print(conf, ME)


        #
        self.size_ME = num
        print(" ME number : ", self.size_ME)

    def get_ME(self, a,b,c,d,J,T):

        val = 0.0
        a_s = str(a)
        b_s = str(b)
        c_s = str(c)
        d_s = str(d)
        J_s = str(J)
        T_s = str(T)

        ab_s = a_s+"-"+b_s
        ba_s = b_s+"-"+a_s
        cd_s = c_s+"-"+d_s
        dc_s = d_s+"-"+c_s

        conf_1 = ab_s+"-"+cd_s + "-" + J_s + "-" + T_s
        conf_2 = ab_s+"-"+dc_s + "-" + J_s + "-" + T_s
        conf_3 = ba_s+"-"+cd_s + "-" + J_s + "-" + T_s
        conf_4 = ba_s+"-"+dc_s + "-" + J_s + "-" + T_s

        conf_11 = cd_s +"-"+ab_s+ "-" + J_s + "-" + T_s
        conf_22 = dc_s +"-"+ab_s+ "-" + J_s + "-" + T_s
        conf_33 = cd_s +"-"+ba_s+ "-" + J_s + "-" + T_s
        conf_44 = dc_s +"-"+ba_s+ "-" + J_s + "-" + T_s

        val_t = self.me_dic.get(conf_1, "None")
        if(val_t != "None"):
            val = val_t
            return val

        val_t = self.me_dic.get(conf_11, "None")
        if(val_t != "None"):
            val = val_t
            return val

        j2c = self.sp_state[c-1].j2
        j2d = self.sp_state[d-1].j2
        phase_t = (j2c+j2d)/2 - J + 1 - T
        phase_2 = -1*self.phase(phase_t)

        val_t = self.me_dic.get(conf_2, "None")
        if(val_t != "None"):
            val = val_t*phase_2
            return val

        val_t = self.me_dic.get(conf_22, "None")
        if(val_t != "None"):
            val = val_t*phase_2
            return val

        j2a = self.sp_state[a-1].j2
        j2b = self.sp_state[b-1].j2
        phase_t = (j2a+j2b)/2 - J + 1 - T
        phase_3 = -1*self.phase(phase_t)

        val_t = self.me_dic.get(conf_3, "None")
        if(val_t != "None"):
            val = val_t*phase_3
            return val

        val_t = self.me_dic.get(conf_33, "None")
        if(val_t != "None"):
            val = val_t*phase_3
            return val

        phase_t1 = (j2a+j2b)/2 - J + 1 - T
        phase_t2 = (j2c+j2d)/2 - J + 1 - T
        phase_4 = self.phase(phase_t1+phase_t2)

        val_t = self.me_dic.get(conf_4, "None")
        if(val_t != "None"):
            val = val_t*phase_4
            return val

        val_t = self.me_dic.get(conf_44, "None")
        if(val_t != "None"):
            val = val_t*phase_4
            return val

        return val



    def print_me(self):
        for key,value in self.me_dic.items():
            print('{key}  :  {value}'.format(key = key, value = value))


    def print_sp(self):
        sp_state = self.sp_state
        for i in range(self.size_sp):
            print(sp_state[i].index,"\t ",sp_state[i].n,"\t ",sp_state[i].l,"\t ",sp_state[i].j2,"\t ",sp_state[i].t2z,"\t ",sp_state[i].E)




#



#
