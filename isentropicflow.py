# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 13:35:56 2017

@author: xajsc
"""
#    def isentropicflow(k,*args,**kwargs):
'''Calculates the isentropic flow properties for inputs k and v2 'any one of [M,P,T,A]'
    
        k=Gamma
        M=Mach Number
        P=PressureRatio (P/P0)
        T=TemperatureRatio (T/T0)
        R=DensityRatio (rho/rho0)
        A=AreaRatio (A/A*)
            One AreaRatio gives two solutions of Mach number.So, Asub for subsonic 
            solution and Asup for supersonic solution
        
        
        Input Pattern: isentropic(k,v2)
        For input : Gamma=1.4 and Mach=2 enter input as (1.4,2) or (1.4,M=2) 
        For TemperatureRatio: T=2.5 or T=[2.5,3] 
        Same as TemperatureRatio for PressureRatio (P), DensityRatio(R,Rho)
        For AreaRatio as input enter 
                : Asub for subsonic solution and 
                  Asup for supersonic solution
        Anish Acharya
'''        
import numpy as np
import scipy

def return_value(param):
    if  __name__=='__main__':
        print("Parameter:",param)
    else:
        pass

    #%% Final Values           
class IsentropicFlow(object):  
    
     def __init__ (self,k,*args,**kwargs): 
        print ("Gamma :", k)

        if type(k) in [int,float]:
            k=[float(k)]
    #%% ## Check if gamma (k) is real value. ##
        if any(np.iscomplex(k)):
            raise TypeError('Error: NonReal Number input in Gamma')
    
    #%% ## Check if gamma (k) is numeric value or not. ##
        if all(isinstance(x, (int, float)) for x in k):
            pass
        else:
            raise TypeError("Error: String input in Gamma")
        
    #%% ## Check if k>1 ##
        if any(k)<1:
            raise ValueError('Error: Gamma less than 1')
    #%% ## Changing k to numpy array ##  
    
        
        k=np.array(k)
        self.k = k
        
        a,b=len(args),len(kwargs)
    #%% ##Check the number of input argument. Required=2.  ##
    
        if a+b >1 or  a>1 or b>1:
            raise ValueError('Too many inputs')
        elif a+b<1:
            raise ValueError('Error: Not enough input arguments')
        elif a==1:
            param='Mach'
            return_value(param)
            M_calc=args[0]
            self.M=M_calc
            
        else:
            key = list(dict.keys(kwargs))  
            value=list(dict.values(kwargs))
#            if type(value) in [int,float]:
#                value=[float(value)]
            value=np.array(value)
            key=[key.lower() for key in kwargs]
            check=key[0]
            self.M = self.mach_calculator(value,check)
     @property       
     def Mach(self):
         return self.M
    
     @property       
     def T0_T(self):
#        M=mach()
        T0_T = (1 + (self.k - 1)/2 * np.power(self.M,2))
        return T0_T    

     @property
     def p0_p(self,*args,**kwargs):     
        p0_p=np.power(self.T0_T,self.k/(self.k-1))
        return p0_p    
     @property
     def rho0_rho(self):
        rho0_rho=np.power(self.T0_T,1/(self.k-1))
        return rho0_rho
     @property
     def A_Astar(self):    
        b=np.divide((self.k+1),2*(self.k-1))
        A_Astar=1/self.M*np.power((2/((self.k+1))*(1+(self.k-1)/2*np.power(self.M,2))),b)
        return A_Astar   
    
     @property
     def T_T0(self):
         return 1/self.T0_T
     
     @property   
     def p_p0(self):
         return 1/self.p0_p
     @property   
     def rho_rho0(self):
         return 1/self.rho0_rho
     
     
     #%%Mach Calculator
     def mach_calculator(self,value,check):
            value=value[0]
#            value=np.array(value)
            if check in ['mach','m','machnumber']:
                param='Mach'
                if any(i<0 for i in value):
                    raise ValueError('Error: Value MissMatch')
                return_value(param)
                M_calc=value
                       
            elif check in ['temp','tempratio','temperatureratio','t','tratio','t_t0']:
               param='TempRatio'
               return_value(param)
               if any(value<0):
                   raise ValueError('Error : "TemperatureRatio" Value less than 0')
               elif any(value>1):
                   raise ValueError('Error : "TemperatureRatio" Value greater than 1')
               else:
                   M_calc=np.power(2*(1/value-1)/(self.k-1),0.5)
               
            elif check in ['pressureratio','p','pratio','p_p0']:
               param='PressureRatio'
               return_value(param)
               if any(value<0):
                   raise ValueError('Error : "PressureRatio" Value less than 0')
               elif any(value>1):
                   raise ValueError('Error : "PressureRatio" Value greater than 1')
               else:
                   b=(self.k-1)/self.k
                   M_calc=np.power((2/(self.k-1)*((1/np.power(value,b))-1)),0.5)
               
            elif check in ['rhoratio','rratio','dratio','densityratio','r','d','rho_rho0','r_r0']:
               param='DensityRatio'
               return_value(param)
               if any(value<0):
                   raise ValueError('Error : "DensityRatio" Value less than 0')
               elif any(value>1):
                   raise ValueError('Error : "DensityRatio" Value greater than 1')
               else:
                   b=(self.k-1)
                   M_calc=np.power((2/(self.k-1)*((1/np.power(value,b))-1)),0.5)
            
            elif check in ['areasubratio','asubratio','asub']:
                param='AreaRatio Subsoinc'
                return_value(param)
                if any(value<1):
                    raise ValueError('Error : "AreaRatio" Value less than 1')
                afunc = lambda x : value-(1/x)*np.power((2/((self.k+1))*(1+(self.k-1)/2*np.power(x,2))),(self.k+1)/(2*(self.k-1)))
                x = np.ones(len(value))*0.0001
                M_calc=scipy.optimize.fsolve(afunc,x)
            elif check in ['areasupratio','asupratio','asup']:
                param='AreaRatio Supersonic'
                return_value(param)
                if any(value<1):
                    raise ValueError('Error : "AreaRatio" Value less than 1')
                afunc = lambda x : value-(1/x)*np.power((2/((self.k+1))*(1+(self.k-1)/2*np.power(x,2))),(self.k+1)/(2*(self.k-1)))
                x = np.ones(len(value))*20
                M_calc=scipy.optimize.fsolve(afunc,x)
            else:
                print()
                raise ValueError('Error:ParameterError \nFor input : Gamma=1.4 and Mach=2 enter input as (1.4,2) or (1.4,M=2)\nFor TemperatureRatio: T=2.5 or T=[2.5,3]\nSame as TemperatureRatio for PressureRatio (P), DensityRatio(R)\nFor AreaRatio as input enter Asub for subsonic solution and Asup for supersonic solution')
            return M_calc  
 #%%      
def mach_angle(M):
    if M<1:
        raise('Mach number must be supersonic')
    else:
        return np.arcsin(1/M)
        
#%%  
emp1=IsentropicFlow(1.4,[6,3,6])    

     
