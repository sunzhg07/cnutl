#!/ usr / bin / env python3
#- * - coding : utf - 8 - * -

import numpy as np
import math
import sys
from numba import jit
from get_ME_gsm import *
from get_ME_usdb import *
from coupling import *
from sympy.physics.wigner import wigner_6j
from sympy.physics.wigner import wigner_9j
from tqdm import tqdm # prograss bar



class Tensor_Decompose :
# Beyond Hartree-Fock in Nuclear Physics; a tale of correlations
#  M. Hjorth-Jensen
    def __init__(self, me, C_cpp):
        self.sp_state = me.sp_state
        self.ME = me
        self.me_K_dic = {}
        self.Coupling = C_cpp
        self.sqrt_2 = math.sqrt(2)
        sqrt_num_list = []
        sqrt_num_list.append(1)
        for i in range(1,1000):
            val = math.sqrt(i)
            sqrt_num_list.append(val)
        self.sqrt_num_list = sqrt_num_list

    def phase(self,N):
        if(N%2 == 0):
            return 1
        else:
            return -1

    def Norm_AB(self,a,b):
        val = 1
        #delta = 1
        if(a==b):
            #delta = 0
            val = self.sqrt_2
        norm = 1.0/val
        return norm
    #@jit
    def NinJ_cal(self,ja,jb,jab,jc,jd,jcd,jac,jbd,J):
        val = 0.0
        # val = wigner_9j(ja,jb,jab,jc,jd,jcd,jac,jbd,J,prec=10)
        val = self.Coupling.NinJ(ja,jb,jab,jc,jd,jcd,jac,jbd,J)
        #print(ja,jb,jab,jc,jd,jcd,jac,jbd,J)
        return val
    def SixJ_cal(self,ja,jb,jc,jd,je,jf):
        val = 0.0
        # val = wigner_9j(ja,jb,jab,jc,jd,jcd,jac,jbd,J,prec=10)
        val = self.Coupling.SixJ(ja,jb,jc,jd,je,jf)
        #print(ja,jb,jab,jc,jd,jcd,jac,jbd,J)
        return val
    #@jit
    def ablsj(self,a,b,L,S,J):
        #return 1.0
        j2a = self.sp_state[a-1].j2
        j2b = self.sp_state[b-1].j2

        l_a = self.sp_state[a-1].l
        l_b = self.sp_state[b-1].l
        j_a = j2a/2.0
        j_b = j2b/2.0

        pf = (j2a+1)*(j2b+1)*(2*L+1)*(2*S+1)
        pf = math.sqrt(pf)
        #pf = self.sqrt_num_list[pf]

        #ninJ = wigner_9j(l_a,0.5,j_a,l_b,0.5,j_b,L,S,J,prec=16)
        ninJ = self.NinJ_cal(l_a,0.5,j_a,l_b,0.5,j_b,L,S,J)
        val = pf*ninJ

        return val

    #@jit
    def Tri(self, ja, jb, jc):
        # // I+J>=K, I+K>=J, J+K>=I,
        # // I/2+J/2+K/2 = INTEGER.
        # // TRI=1, WHEN TRIADIC CONDITION IS FULFILLED, TRI=-1 OTHERWISE
        j2a = int(2*ja)
        j2b = int(2*jb)
        j2c = int(2*jc)
        L2 = j2a + j2b + j2c;
        if (int(L2) != (L2 / 2) * 2):
            return -1
        L2 = L2 / 2
        if (j2a * j2b * j2c < 0):
            return -1
        if ((L2 - j2a) * (L2 - j2b) * (L2 - j2c) < 0):
            return -1
        return 1
        # val = 1
        # if(J12 < abs(j1-j2)):
        #     val = -1
        # if(J12 > abs(j1+j2)):
        #     val = -1
        # return val

    def ab_group_nl(self,a,b):
        val = 1
        if(self.sp_state[a-1].n != self.sp_state[b-1].n):
            val = -1
            return val
        if(self.sp_state[a-1].l != self.sp_state[b-1].l):
            val = -1
            return val
        return val
    #@jit
    def build(self, K):
        for conf,val in tqdm(self.ME.me_dic.items()):
            #print(conf,val)
            s = conf.split("-")
            #print(s)
            a = int(s[0])
            b = int(s[1])
            c = int(s[2])
            d = int(s[3])
            J = int(s[4])
            T = int(s[5])

            val_k_me = self.cal(a,b,c,d,J,T,K)
            conf_k = conf + "-" + str(K)
            self.me_K_dic[conf_k] = val_k_me


    def print_me(self):
        for conf,val in self.me_K_dic.items():
            s = conf.split("-")
            a = int(s[0])
            b = int(s[1])
            c = int(s[2])
            d = int(s[3])
            J = int(s[4])
            T = int(s[5])
            K = int(s[6])

            print(a,b,c,d,"\t ",J,"\t ",T,"\t ",K,"\t " ,val)

    def print_me_f(self,filename):
        with open(filename, 'r+') as f:
            #f.write(mobile+'\n')
            for conf,val in self.me_K_dic.items():
                s = conf.split("-")
                a = int(s[0])
                b = int(s[1])
                c = int(s[2])
                d = int(s[3])
                J = int(s[4])
                T = int(s[5])
                K = int(s[6])

                f.writelines([str(a)," ",str(b)," ",str(c)," ",str(d),"\t ",str(J),"\t ",str(T),"\t ",str(K),"\t " ,str(val), "\n"])
        f.close()

    #@jit
    def cal_print(self,a,b,c,d,J,T,K):

        val = 0.0
        size_sp = self.ME.size_sp

        phase_pf = self.phase(J)
        pf = (2*K + 1)

        l_a = self.sp_state[a-1].l
        l_b = self.sp_state[b-1].l
        l_c = self.sp_state[c-1].l
        l_d = self.sp_state[d-1].l

        L_ab_min = abs(l_a - l_b)
        L_ab_max = abs(l_a + l_b)

        L_cd_min = abs(l_c - l_d)
        L_cd_max = abs(l_c + l_d)

        S_ab_min = 0
        S_ab_max = 1
        S_cd_min = 0
        S_cd_max = 1

        part_1 = 0.0
        for S_ab in range(S_ab_min,S_ab_max+1):
            for S_cd in range(S_cd_min,S_cd_max+1):
                if(self.Tri(S_ab,S_cd,K)<0):
                    continue
                for L_ab in range(L_ab_min,L_ab_max+1):
                    if(self.Tri(L_ab,S_ab,J)<0):
                        continue
                    for L_cd in range(L_cd_min,L_cd_max+1):
                        if(self.Tri(L_cd,S_cd,J)<0):
                            continue
                        if(self.Tri(L_ab,L_cd,K)<0):
                            continue
                        #part_1_t =  wigner_6j(L_ab, S_ab, J, S_cd, L_cd, K,prec=16)
                        part_1_t = self.SixJ_cal(L_ab, S_ab, J, S_cd, L_cd, K)
                        if(part_1_t == 0.0):
                            continue
                        part_1_t = part_1_t * self.ablsj(a,b,L_ab,S_ab,J)*self.ablsj(c,d,L_cd,S_cd,J)
                        if(part_1_t == 0.0):
                            continue

                        J_ab_min = abs( L_ab - S_ab )
                        J_ab_max = abs( L_ab + S_ab )
                        J_cd_min = abs( L_cd - S_cd )
                        J_cd_max = abs( L_cd + S_cd )

                        J_p_min = max(J_ab_min, J_cd_min)
                        J_p_max = min(J_ab_max, J_cd_max)

                        part_2 = 0.0
                        for J_p in range(J_p_min,J_p_max+1):
                            pf_2 = self.phase(J_p)*(2*J_p+1)
                            part_2_t = pf_2 * self.SixJ_cal(L_ab, S_ab, J_p, S_cd, L_cd, K)
                            if(part_2_t == 0.0):
                                continue
                            part_3 = 0.0
                            for i in range(size_sp):
                                i = i+1
                                if (self.ab_group_nl(i,a)<0):
                                    continue
                                for j in range(size_sp):
                                    j = j+1
                                    if (self.ab_group_nl(j,b)<0):
                                        continue
                                    pf_ijLSJp = self.ablsj(i,j,L_ab, S_ab, J_p)
                                    if(pf_ijLSJp == 0.0):
                                        continue
                                    for k in range(size_sp):
                                        k = k+1
                                        if (self.ab_group_nl(k,c)<0):
                                            continue
                                        for l in range(size_sp):
                                            l = l+1
                                            if (self.ab_group_nl(l,d)<0):
                                                continue
                                            pf_klLSJp = self.ablsj(k,l,L_cd, S_cd, J_p)
                                            if(pf_klLSJp == 0.0):
                                                continue

                                            val_me = self.ME.get_ME(i,j,k,l,J_p,T)
                                            #print("--- ",i,j,k,l,"\t ",J_p,"\t ",T)
                                            #print("\t\t ",val_me, pf_ijLSJp, pf_klLSJp)
                                            part_3_t = val_me*pf_ijLSJp*pf_klLSJp

                                            part_3 = part_3 + part_3_t
                                        #l
                                    #k
                                #j
                            #i
                            part_2 = part_2 + part_2_t*part_3
                        #J_p
                        part_1 = part_1 + part_1_t * part_2
                    #L_cd
                #L_ab
            #S_cd
        #S_ab
        val = phase_pf * pf * part_1
        return val

    #@jit
    def cal(self,a,b,c,d,J,T,K):

        val = 0.0
        size_sp = self.ME.size_sp

        phase_pf = self.phase(J)
        pf = (2*K + 1)

        l_a = self.sp_state[a-1].l
        l_b = self.sp_state[b-1].l
        l_c = self.sp_state[c-1].l
        l_d = self.sp_state[d-1].l

        L_ab_min = abs(l_a - l_b)
        L_ab_max = abs(l_a + l_b)

        L_cd_min = abs(l_c - l_d)
        L_cd_max = abs(l_c + l_d)

        norm_ab = self.Norm_AB(a,b)
        norm_cd = self.Norm_AB(c,d)

        S_ab_min = 0
        S_ab_max = 1
        S_cd_min = 0
        S_cd_max = 1

        part_1 = 0.0
        for S_ab in range(S_ab_min,S_ab_max+1):
            for S_cd in range(S_cd_min,S_cd_max+1):
                if(self.Tri(S_ab,S_cd,K)<0):
                    continue
                for L_ab in range(L_ab_min,L_ab_max+1):
                    if(self.Tri(L_ab,S_ab,J)<0):
                        continue
                    for L_cd in range(L_cd_min,L_cd_max+1):
                        if(self.Tri(L_cd,S_cd,J)<0):
                            continue
                        if(self.Tri(L_ab,L_cd,K)<0):
                            continue
                        # part_1_t =  wigner_6j(L_ab, S_ab, J, S_cd, L_cd, K,prec=16)
                        part_1_t =  self.SixJ_cal(L_ab, S_ab, J, S_cd, L_cd, K)
                        if(part_1_t == 0.0):
                            continue
                        part_1_t = part_1_t * self.ablsj(a,b,L_ab,S_ab,J)*self.ablsj(c,d,L_cd,S_cd,J)
                        if(part_1_t == 0.0):
                            continue

                        J_ab_min = abs( L_ab - S_ab )
                        J_ab_max = abs( L_ab + S_ab )
                        J_cd_min = abs( L_cd - S_cd )
                        J_cd_max = abs( L_cd + S_cd )

                        J_p_min = max(J_ab_min, J_cd_min)
                        J_p_max = min(J_ab_max, J_cd_max)

                        part_2 = 0.0
                        for J_p in range(J_p_min,J_p_max+1):
                            pf_2 = self.phase(J_p)*(2*J_p+1)
                            part_2_t = pf_2 * self.SixJ_cal(L_ab, S_ab, J_p, S_cd, L_cd, K)
                            if(part_2_t == 0.0):
                                continue
                            part_3 = 0.0
                            for i in range(size_sp):
                                i = i+1
                                if (self.ab_group_nl(i,a)<0):
                                    continue
                                for j in range(size_sp):
                                    j = j+1
                                    if (self.ab_group_nl(j,b)<0):
                                        continue
                                    pf_ijLSJp = self.ablsj(i,j,L_ab, S_ab, J_p)
                                    if(pf_ijLSJp == 0.0):
                                        continue
                                    norm_ij = self.Norm_AB(i,j)
                                    for k in range(size_sp):
                                        k = k+1
                                        if (self.ab_group_nl(k,c)<0):
                                            continue
                                        for l in range(size_sp):
                                            l = l+1
                                            if (self.ab_group_nl(l,d)<0):
                                                continue
                                            pf_klLSJp = self.ablsj(k,l,L_cd, S_cd, J_p)
                                            if(pf_klLSJp == 0.0):
                                                continue
                                            norm_kl = self.Norm_AB(k,l)
                                            val_me = self.ME.get_ME(i,j,k,l,J_p,T)
                                            part_3_t = 1.0/(norm_ij*norm_kl)*val_me*pf_ijLSJp*pf_klLSJp

                                            part_3 = part_3 + part_3_t
                                        #l
                                    #k
                                #j
                            #i
                            part_2 = part_2 + part_2_t*part_3
                        #J_p
                        part_1 = part_1 + part_1_t * part_2
                    #L_cd
                #L_ab
            #S_cd
        #S_ab
        val = norm_ab*norm_cd*phase_pf * pf * part_1
        return val
