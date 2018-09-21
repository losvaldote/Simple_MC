#
# Early Dark Energy cosmology. 
#
#Se considera un espacio plano, es decir, curvatura cero.
#Vamos a introducir un nuevo parametro asociado con la ecuacion de estado de DM.
#Esta es una combinacion de los modelos wDMCosmology.py y EarlyDECosmology.py
#Para materia oscura el parametro es wMO.

import math as N
from LCDMCosmology import *

class Early_DE_DMCosmology(LCDMCosmology):
    def __init__(self, varyw=True, varyOde=True, userd_DE=True, varywMO=True):
        ## two parameters: Om and h

        self.userd_DE = userd_DE
	print 'userd', userd_DE	

        self.varyw= varyw
        self.varyOde=varyOde
	self.varywMO=varywMO
  
        self.w0=  w_par.value
        self.Ode=Ode_par.value
	self.wMO=wMO_par.value

        self.oC=LCDMCosmology()
        LCDMCosmology.__init__(self)

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        if (self.varyw): l.append(w_par)
        if (self.varyOde): l.append(Ode_par)
	if (self.varywMO): l.append(wMO_par)
        return l


    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        self.oC.updateParams(pars)
        if not ok:
            return False
        for p in pars:
            if p.name=="w":
               self.w0=p.value
            if p.name=="Ode":
               self.Ode=p.value  
	    elif p.name=="wMO":
		self.wMO=p.value
        return True

    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        NuContrib=self.NuDensity.rho(a)/self.h**2
	self.Odm = self.Ocb -self.Obh2/(self.h**2) #Se agrega por DM
	self.Ob = self.Ocb -self.Odm               #Se agrega por DM
        Omega_d0 = 1.0-self.Ocb-self.Omrad-self.NuDensity.rho(1.0)/self.h**2
        factor= self.Ob/a**3+self.Odm*a**(-3*(1.0+self.wMO))+self.Omrad/a**4+NuContrib
        Omega_d =(Omega_d0-self.Ode*(1.0-a**(-3*self.w0)))/(Omega_d0+factor*a**(3*(1.0+self.w0)))+ self.Ode*(1.0-a**(-3*self.w0)) 
        return factor/(1.0-Omega_d)


    def prefactor(self):
	if self.userd_DE:
		self.rd = self.oC.rd*((1.-self.Ode)**(0.5))
	else:
		self.rd = self.oC.rd
	return 	self.c_/(self.rd*self.h*100)



#for printing purposes only
    def Omega_de(self, a):
        NuContrib=self.NuDensity.rho(a)/self.h**2
        Omega_d0 = 1.0-self.Ocb-self.Omrad-self.NuDensity.rho(1.0)/self.h**2 
        factor= self.Ob/a**3+self.Odm*a**(-3*(1.0+self.wMO))+self.Omrad/a**4+NuContrib
        return (Omega_d0-self.Ode*(1.0-a**(-3*self.w0)))/(Omega_d0+factor*a**(3*(1.0+self.w0)))+ self.Ode*(1.0-a**(-3*self.w0))




