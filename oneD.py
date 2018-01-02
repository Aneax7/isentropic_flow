# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 11:37:20 2017

@author: xajsc
"""  
p=1890 # Pascal
T=450  # Kelvin
M=1.5          # Mach number

k=1.4  # Specific heat ratio
R=1716 

import speed as sp
import isentropicflow as iflow
stag=iflow.IsentropicFlow(k,T=[0.25])
#
#class stagnation1D:
#    def __init__(self):
#        pass    

def density(p,T,R):
    density=p/(R*T)
    return density

## Initial density ##
rho=density(p,T,R)
# Stagnation stags

# Mach Number
M=stag.Mach
print('Mach:',M)

#Pressure Ratio
p_p0=stag.p_p0
print('P_P0:',p_p0)

# Temperature Ratio
T_T0=stag.T_T0
print('T_T0:',T_T0)

# Density Ratio
rho_rho0=stag.rho_rho0
print('rho_rho0:',rho_rho0)

# Area Ratio
A_Astar=stag.A_Astar
print('A_Astar:',A_Astar)

p0=1/p_p0*p
T0=1/T_T0*T
rho0=1/rho_rho0*rho
#rho0=1/rho_rho0*rho
print('Total Pressure: {} \nTotal Temperature: {}\n'.format(p0,T0)) 
print('rho0:',rho0)

#%% Critical stags
crit=iflow.IsentropicFlow(k,1)

M_crit=crit.Mach

# Pressure Ratio
pstar_p0=crit.p_p0
print('p_pcrit:',pstar_p0)

# Temperature Ratio
Tstar_T0=crit.T_T0
print('T_Tcrit:',Tstar_T0)

# Density Ratio
rho_rhocrit=crit.rho_rho0
print('rho_rhocrit',rho_rhocrit)

# Area Ratio
A_Astar=crit.A_Astar
print('A_Astar',A_Astar)

T_crit=Tstar_T0*T0
p_crit=(pstar_p0)*p0

#rho0=1/rho_rho0*rho
print('Pcritical:{}\nTcritical:{}'.format(p_crit,T_crit)) 
   

#%% Flow Peoperties 
a,V=sp.speed(k,T,R,M) ## Speed of Sound and Flow Velocity
print('Speed of Sound: {}\nFlow Velocity: {}'.format(a,V))
    

