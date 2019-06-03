## This is the model for the metric g_3 for the group SU(2).
## We have to add two parameters b and \dot{b}.
## As a first approcimation we consider that these parameters are numberlike
## and they are not a scale factor.
## We consider a general curvature.


import math as N
from LCDMCosmology import *

class fibradosCosmology(LCDMCosmology):
    def __init__(self, varyb=True, varybdot=True, varybddot=True, varyadot=True, varyOk=True):
        ## two parameters: Om and h

        self.varyb=varyb
        self.varybdot=varybdot
        self.varyOk=varyOk
	self.varybddot=varybddot
	self.varyadot=varyadot

        self.Ok=Ok_par.value   
        self.b=b_par.value
        self.bdot=bdot_par.value
	self.bddot=bddot_par.value
	self.adot=adot_par.value
        LCDMCosmology.__init__(self)

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        if (self.varyb): l.append(b_par)
        if (self.varybdot): l.append(bdot_par)
	if (self.varybddot): l.append(bddot_par)
	if (self.varyadot): l.append(adot_par)
        if (self.varyOk): l.append(Ok_par)
        return l


    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        for p in pars:
            if p.name=="b":
                self.b=p.value
            elif p.name=="bdot":
                self.bdot=p.value
	    elif p.name=="bddot":
		self.bddot=p.value
	    elif p.name=="adot":
		self.adot=p.value
            elif p.name=="Ok":
                self.Ok=p.value
                self.setCurvature(self.Ok)
                if (abs(self.Ok)>1.0):
                   return False
        return True

    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        NuContrib=self.NuDensity.rho(a)/self.h**2
	#Pe=(3.*(self.bddot/self.b) + 6.*((self.adot*self.bdot)/a*self.b) + 3.*(self.bdot/self.b)**2. + 3./(8.*self.b**2))
	#rhoe=-3.*(self.bdot/self.b)**2. - 9.*((self.adot*self.bdot)/a*self.b) - 3./(8.*self.b**2)
	Pe = a + 16.*self.b*self.adot*self.bdot + 8.*a*self.bdot**2 + 8.*a*self.b*self.bddot
	rhoe = a + 8.*self.bdot*(3.*self.b*self.adot + a*self.bdot)
## we es la ecuacion de estado de la energia oscura para este grupo y metrica.
        we= -Pe/rhoe
        return ((self.Ocb/a**3+self.Ok/a**2+self.Omrad/a**4+NuContrib+(1.0-self.Om-self.Ok)*a**(-3.0*(1.0+we)))/self.b**3)

