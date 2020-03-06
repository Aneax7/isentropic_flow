# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 18:50:21 2017

@author: xajsc
Calculates the speed of sound and flow velocity

Inputs:
    k           =specific heat ratio
    R           =universal gas constant
    T_initial   = Initial Temperature
    M           =Flow Mach Number
"""
def speed(k,T_initial,R,M):
     ## Speed of Sound##
    a=(k*R*T_initial)**(1/2)
    V=M*a
    return a,V