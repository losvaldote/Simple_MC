## This is CDM cosmology with w, wa and w for DM
##No se considera curvatura.
##Vamos a introducir un nuevo parametro asociado con la ecuacion de estado de DM
##Esta es una combinacion de los modelos wDMCosmology.py y owa0CDMCosmology.py
##Para materia oscura el parametro es wMO.
##Se llama EKCDMCosmology por mis sobrinas E y K.


import math as N
from LCDMCosmology import *

class EKCDMCosmology(LCDMCosmology):
    def __init__(self, varyw=True, varywa=True, varyOk=False, varywMO=True):
        ## two parameters: Om and h

        self.varyw=varyw
        self.varywa=varywa
        self.varyOk=varyOk
	self.varywMO=varywMO

        self.Ok=Ok_par.value   
        self.w0=w_par.value
        self.wa=wa_par.value
	self.wMO=wMO_par.value
        LCDMCosmology.__init__(self)

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        if (self.varyw): l.append(w_par)
        if (self.varywa): l.append(wa_par)
        if (self.varyOk): l.append(Ok_par)
	if (self.varywMO): l.append(wMO_par)
        return l


    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        for p in pars:
            if p.name=="w":
                self.w0=p.value
            elif p.name=="wa":
                self.wa=p.value
	    elif p.name=="wMO":
		self.wMO=p.value
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
	self.Odm = self.Ocb -self.Obh2/(self.h**2)
	self.Ob = self.Ocb -self.Odm 
        rhow= a**(-3*(1.0+self.w0+self.wa))*N.exp(-3*self.wa*(1-a))
        return (self.Ob/a**3+self.Ok/a**2+ self.Odm*a**(-3*(1.0+self.wMO)) +self.Omrad/a**4+NuContrib+(1.0-self.Om-self.Ok)*rhow)

