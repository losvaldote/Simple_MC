## This is CDM cosmology with Ok.
## Cold Dark Matter and Cosmological Constant.


import math as N
from LCDMCosmology import *

class curvatureLCDM(LCDMCosmology):
    def __init__(self, varyOk=True):
        ## two parameters: Om and h

        self.varyOk=varyOk

        self.Ok=Ok_par.value   
        
        LCDMCosmology.__init__(self)

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        if (self.varyOk): l.append(Ok_par)
        return l


    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        for p in pars:
              if p.name=="Ok":
                self.Ok=p.value
                self.setCurvature(self.Ok)
                if (abs(self.Ok)>1.0):
                   return False
        return True

    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        NuContrib=self.NuDensity.rho(a)/self.h**2
        return (self.Ocb/a**3+self.Ok/a**2+self.Omrad/a**4+NuContrib+(1.0-self.Om-self.Ok))

