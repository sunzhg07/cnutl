#!/ usr / bin / env python3
#- * - coding : utf - 8 - * -
import Coupling_Module

C_cpp = Coupling_Module.Coupling()

class Coupling_F_cpp :
    def __init__(self):
        self.A = 0

    def frac(self,n):
        val = C_cpp.frac(n)
        return val

    def CG(self, ja, jb, jc, jam, jbm, jcm):
        j2a = int(2*ja)
        j2b = int(2*jb)
        j2c = int(2*jc)
        j2am = int(2*jam)
        j2bm = int(2*jbm)
        j2cm = int(2*jcm)
        val =  C_cpp.CG(j2a, j2b, j2c, j2am, j2bm, j2cm)
        return val

    def SixJ(self, ja, jb, jc,jd, je, jf):
        j2a = int(2*ja)
        j2b = int(2*jb)
        j2c = int(2*jc)
        j2d = int(2*jd)
        j2e = int(2*je)
        j2f = int(2*jf)
        #ll = ctypes.cdll.LoadLibrary
        #print(self.lib_file)
        #C_cpp = ll(self.lib_file)
        #val = self.C_cpp.add_c(1,2)#SixJ(j2a, j2b, j2c, j2d, j2e, j2f)
        val = C_cpp.SixJ(j2a, j2b, j2c, j2d, j2e, j2f)
        return val

    def NinJ(self, ja, jb, jc,jd, je, jf, jg, jh, ji):
        j2a = int(2*ja)
        j2b = int(2*jb)
        j2c = int(2*jc)
        j2d = int(2*jd)
        j2e = int(2*je)
        j2f = int(2*jf)
        j2g = int(2*jg)
        j2h = int(2*jh)
        j2i = int(2*ji)
        #ll = ctypes.cdll.LoadLibrary
        #print(self.lib_file)
        #C_cpp = ll(self.lib_file)
        #val = self.C_cpp.add_c(1,2)#SixJ(j2a, j2b, j2c, j2d, j2e, j2f)
        val = C_cpp.NinJ(j2a, j2b, j2c, j2d, j2e, j2f, j2g, j2h, j2i)

        return val
