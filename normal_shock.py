"""
Created on Fri Dec 15 17:08:02 2017

@author: xajsc

Calculates the flow parameters after the shock for given input flow parameters.

"""
import numpy as np
#%%
def shock_parameters(p_1,T_1,k=1.4,R=287,**kwargs):
    print(kwargs)
    key = list(dict.keys(kwargs))
    key=key[0]
    value=list(dict.values(kwargs))
    value=np.array(value)
    if key=='M_1':
        M_1=value
    elif key in ['M1_star','M1*','M_1*']:
        M_1=2/((k+1)/value**2-(k-1))
    return shock_relations(p_1,T_1,M_1,k,R)

#%%
def shock_relations(p_1,T_1,M_1,k,R) :   
    ## Normal Shock Relations##
    rho_1=p_1/(R*T_1) #Inital Density kg/m3
    
    ## Final Mach Number
    M_2=((1+((k-1)/2)*M_1**2)/(k*M_1**2-(k-1)/2))**0.5
#    print('\nM2=',M_2)
    
    ## Final Parameters
    rho_2=rho_1*((k+1)*M_1**2)/(2+(k-1)*M_1**2)
    p_2=p_1*(1+(2*k/(k+1))*(M_1**2-1))
    T_2=(p_2/p_1)*(rho_1/rho_2)*T_1
     
    return M_1,M_2,p_2,T_2,rho_2

#def round_function(*args):
#    arg2=[]
#    for i in args:
#        i=round(i,3)
#        arg2.append(i)
#        return arg2
#%%
p_1=0.5*101325  # Initial Pressure (Pascal)
T_1=200         # Initial Temperature (Kelvin)
k=1.4           # Specific heat ratio
R=287           #J/kg-K

M_1,M_2,p_2,T_2,rho_2=shock_parameters(p_1,T_1,M1_star=5)
#ppp=round_function(M_2,p_2,T_2,rho_2)
import speed as sp
a_2,u_2=sp.speed(k,T_2,R,M_2)

print('Final Pressure=',p_2/101325,'atm\nFinal Velocity=',u_2,'m/s')
print('Final Temperature=',T_2,'K\nFinal Speed of Sound=',a_2,'m/s')
import isentropicflow as ifl
emp1=ifl.IsentropicFlow(k,M_1)
emp2=ifl.IsentropicFlow(k,M_2)
p01=emp1.p0_p*p_1
p02=emp2.p0_p*p_2
del_s=-R*np.log(p01/p02)