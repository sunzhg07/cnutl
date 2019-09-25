#!/ usr / bin / env python3
#- * - coding : utf - 8 - * -

import numpy as np
import sys

class SP :
    'single particle state'
#sp_num = 0;
    def __init__(self, n, l, j2, index) :
        self.index = index
        self.n = n
        self.l = l
        self.j2 = j2
#Basis_SP.sp_num += 1

class Get_ME_usdb_int :
    def __init__(self):
        self.sp_state =[]
        self.me_dic = {}

    def phase(self,N):

        if(N%2 == 0):
            return 1
        else:
            return -1

    def read_data(self,filename):
        f = open(filename)
        line = f.readline()
        line = f.readline()
#-- - -- - get SP -- - -- -
        line = f.readline()
        size_sp = line.splitlines()[0].split()[0]
        size_sp = int(size_sp)
        print(size_sp)
        for i in range(size_sp):
            line = f.readline()
            A = line.splitlines()[0]
            B = A.split()
            index = int(B[0])
            n=int(B[1])
            l=int(B[2])
            j2=int(B[3])
            sp_t = SP(n,l,j2,index)
            self.sp_state.append(sp_t)
            print(index, n, l, j2)
#-- - -- - get ME -- - -- -
        line = f.readline()
        line = f.readline()
        line = f.readline()
        line = f.readline()
        size_ME = line.splitlines()[0].split()[0]
        size_ME = int(size_ME)
        print(size_ME)
        num = 0
        line = f.readline()
        while line:
            #print(line)
            A = line.splitlines()[0]
            B = A.split()

            conf = B[0] +"-" + B[1] +"-"+ B[2] +"-"+ B[3]+"-"+ B[4]+"-"+ B[5]
            self.me_dic[conf] = float(B[6])

            num += 1
            line = f.readline()

        if (num != size_ME):
            print("Wrong @ get_data(self,filename): ")
            sys.exit()
        self.size_sp = size_sp
        self.size_ME = size_ME

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
            print(sp_state[i].index,"\t ",sp_state[i].n,"\t ",sp_state[i].l,"\t ",sp_state[i].j2)

class Get_ME_usdb_snt :
    def __init__(self):
        self.sp_state =[]
        self.me_dic = {}
        self.phase = Get_ME_usdb_int().phase(N)
        self.print_me = Get_ME_usdb_int().print_me()
        self.print_sp = Get_ME_usdb_int().print_sp()

    def read_data(self,filename):
        f = open(filename)
        line = f.readline()
        line = f.readline()
        line = f.readline()
        line = f.readline()
#-- - -- - get SP -- - -- -
        line = f.readline()
        size_sp = line.splitlines()[0].split()[0]
        size_sp = int(size_sp)
        print(size_sp)
        for i in range(size_sp):
            line = f.readline()
            A = line.splitlines()[0]
            B = A.split()
            index = int(B[0])
            n=int(B[1])
            l=int(B[2])
            j2=int(B[3])
            sp_t = SP(n,l,j2,index)
            self.sp_state.append(sp_t)
            print(index, n, l, j2)
#-- - -- - get ME -- - -- -
        line = f.readline()
        line = f.readline()
        line = f.readline()
        line = f.readline()
        size_ME = line.splitlines()[0].split()[0]
        size_ME = int(size_ME)
        print(size_ME)
        num = 0
        line = f.readline()
        while line:
            #print(line)
            A = line.splitlines()[0]
            B = A.split()

            conf = B[0] +"-" + B[1] +"-"+ B[2] +"-"+ B[3]+"-"+ B[4]+"-"+ B[5]
            self.me_dic[conf] = float(B[6])

            num += 1
            line = f.readline()

        if (num != size_ME):
            print("Wrong @ get_data(self,filename): ")
            sys.exit()
        self.size_sp = size_sp
        self.size_ME = size_ME

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


#



#
