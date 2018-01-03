"""
Created on Fri Dec 15 17:08:02 2017

@author: xajsc

Calculates the flow parameters after the shock for given input flow parameters.

"""
import numpy as np
import matplotlib.pyplot as plt
import math
#%%
def shock_parameters(k,*args,**kwargs):
    print(kwargs)
    if len(kwargs)==0:
        M_1=args
    else:
        key = list(dict.keys(kwargs))
        key=key[0]
        value=list(dict.values(kwargs))
        value=np.array(value)
        if key=='M_1':
            M_1=value
        elif key in ['M1_star','M1*','M_1*']:
            M_1=2/((k+1)/value**2-(k-1))
    return shock_relations(M_1,k)

#%%
def shock_relations(M_1,k) :   
    ## Normal Shock Relations#
    ## Final Mach Number
    M_2=((1+((k-1)/2)*M_1**2)/(k*M_1**2-(k-1)/2))**0.5
#    print('\nM2=',M_2)
    
    ## Final Parameters
    rho=((k+1)*M_1**2)/(2+(k-1)*M_1**2)
    p=(1+(2*k/(k+1))*(M_1**2-1))
    T=(p)/(rho)

     
    return M_1,M_2,p,T,rho

#%%
def entropy_calc(M_1,M_2,p_1,p_2):
    import isentropicflow as ifl
    emp1=ifl.IsentropicFlow(k,M_1)
    emp2=ifl.IsentropicFlow(k,M_2)
    p01=emp1.p0_p*p_1
    p02=emp2.p0_p*p_2
    return -R*np.log(p02/p01)
#%%
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
rho_1=p_1/(R*T_1) #Inital Density kg/m3

M_1,M_2,p,T,rho=shock_parameters(k,M_1=3)
rho_2=rho_1*rho
p_2=p_1*p
T_2=T*T_1
#ppp=round_function(M_2,p_2,T_2,rho_2)
import speed as sp
a_2,u_2=sp.speed(k,T_2,R,M_2)

#print('Final Pressure=',p_2/101325,'atm\nFinal Velocity=',u_2,'m/s')
#print('Final Temperature=',T_2,'K\nFinal Speed of Sound=',a_2,'m/s')

del_s=entropy_calc(M_1,M_2,p_1,p_2)
#print("Entropy Change",del_s)

#%% Plot
M1,M2,P2,T2,rho2,P02=[],[],[],[],[],[]
M=1
while M<=4.2:
    M_1,M_2,p_2,T_2,rho_2=shock_relations(M,k)
    M1.append(M)
    M2.append(M_2)
    P2.append(p_2)
    T2.append(T_2)
    rho2.append(rho_2)
    del_s=entropy_calc(M_1,M_2,p_1,p_2)
    P0=math.exp(-del_s/R)
    P02.append(P0)
    M=M+0.2
    
plt.figure(1)
plt.plot(1)
plt.plot(M1,P2,label='P2/P1')
plt.plot(M1,T2,label='T2/T1')
plt.plot(M1,rho2,label='rho2/rho1')
plt.legend(bbox_to_anchor=(0, 1), loc=2, borderaxespad=0.) 
plt.axis([1, 4, 0, 20])

plt.figure(2)
plt.plot(M1,M2,label='Mach')
plt.plot(M1,P02,label='Total P2/P1')
plt.axis([1, 4, 0, 1])
plt.legend(bbox_to_anchor=(0.5, 1), loc=2, borderaxespad=0.) 
