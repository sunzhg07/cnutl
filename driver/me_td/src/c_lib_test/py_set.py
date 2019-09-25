#!/usr/bin/env python
#import ctypes
#import coupling
import helloworld
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
a = 2
b = 2
hw = helloworld.helloworld("lanyulei", 18)
hw.printinfo()
#c = oupling.add(a,b)
#print(c)
