#!/usr/bin/env python
import ctypes
#import coupling
import Coupling_Module
from sympy.physics.wigner import clebsch_gordan
from sympy.physics.wigner import wigner_9j
# ll = ctypes.cdll.LoadLibrary
# lib = ll('./coupling.so')
#cpp = ctypes.CDLL('./libcoupling.so')
#
# class coupling(object):
#     def __init__(self):
#         self.add = lib.add()
#
# #coupling c

#B = cpp.add_c(1,2)
# 设置sum()函数传入参数的类型
# lib.add_c.argtypes = [ctypes.c_double, ctypes.c_double]
# # 这是sum()函数返回参数的类型
# lib.add_c.restype = ctypes.c_double
#
# a = ctypes.c_double(2.0)
# b = ctypes.c_double(2.0)
a = 2.1
b = 2.3
# hw = helloworld.helloworld("lanyulei", 18)
# hw.printinfo()
CP = Coupling_Module.Coupling()
c = CP.add(a,b)
ja = 1.5
jam = 1.5
jb = 0.5
jbm = -0.5
jc = 1
jcm = 1

cg_0 = clebsch_gordan(ja,jb,jc,jam,jbm,jcm)
cg_1 = CP.CG(int(2*ja),int(2*jb),int(2*jc),int(2*jam),int(2*jbm),int(2*jcm))

print("cg_0 = ",cg_0, "\t cg_1 = ",cg_1)


print(c)
