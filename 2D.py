# -*- coding: utf-8 -*-
"""
Created on Sep 9 
@author: Darly Castro
"""
#Librerias
import numpy as np
from scipy.integrate import simps
from scipy import interpolate 
from filon import cos_integral


def Vpot(rd,re)
    Vpot = 0
    xd = rd[0,:]
    yd = rd[
    Vpot = - 2.0/np.sqrt()
    return Vpot


def int_vpot(Lx,Ly,Nelectron,RD_grid,RE_grid): 
        
    #----------- Donor position--------------   
    Ndonor= len(RD_grid)

    #---------# Electron position-------------
    
    xe = np.linspace(0, Lx, Nelectron)
    ye = np.linspace(0, Ly, Nelectron)
    Re_x, Re_y = np.meshgrid(xe, xy)

    xE_grid = np.linspace(0, sf, Nelectron) 
    phi_e=np.zeros(Nelectron) # Electron position
    #SD_SF= np.zeros(Ndonor)
    #-------------------   # Potential#  --------------------------------
    VS=np.zeros(shape=(Nelectron,Nelectron,Ndonor,Ndonor))
    kindex=np.zeros(Nelectron,Nelectron)
    for i in range(1,  Ndonor+1):
        
    for j  in range(1, Ndonor+1):
        sd = SD_grid[j-1]
    #    phi_D[j-1]= fiss(sd,Long_arc,phi)
        for k in range(1, Ndonor+1):
                se =SE_grid[k]
             #   phi_e[k]= fiss(se,Long_arc,phi)
                kindex[k]= (k*np.pi/Lx)
                VS[f,j-1,k] = Vpot(se,sd,CAMP_F[f],teta,a,b,Long_arc,phi)       
    return VS,kindex,SE_grid,SD_SF,Ndonor,phi_e,nf


def coef(VS,kindex,sf,nf,Nelectron,Ndonor):
    ds= sf/(Nelectron)
    COEF=np.zeros((nf,Ndonor,Nelectron))
    for j in range(0,Ndonor):
        for i in range(0,nf):
            COEF[i,j,:]=(1/sf)*(cos_integral(VS[i,j,:],ds, kindex, x0=0.0, axis=0))
    return COEF

def schr(COEF,sf,nF,Nelectron,Ndonor):
    nd=int(Nelectron/2)
    AR=np.zeros(shape=(nd,nd)) # Empty Array "2Dimencional"
    eigenvalues = np.zeros(shape=(nF,Ndonor,nd))
    eigenvectors = np.zeros(shape=(nF,Ndonor,nd,nd))

    for j in range(0,nF):
        for s in range(0,Ndonor):
            for k in range(0,nd):
                for k1 in range(0,nd):
                    ka=abs(k-k1)
                    AR[k,k1]= 0.5*(COEF[j,s,ka] - COEF[j,s,k+k1])
                    if (k == k1):
                        AR[k,k1]=AR[k,k1]+(np.pi*k/sf)**2
            eigenvalues[j,s] ,eigenvectors[j,s] = np.linalg.eigh(AR)
    return eigenvalues,eigenvectors,nd

#  w : (…, M) ndarray  #The eigenvalues in ascending order, each repeated according to its multiplicity.

#  v : {(…, M, M) ndarray, (…, M, M) matrix}  #The column v[:, i] is the normalized eigenvector corresponding to the eigenvalue w[i]. Will return a matrix object if a is a matrix object



#--------- Wave Function: -----------
#VecR = Cn = Eigenvectors 
# 1 component : different electric fields  : nF: 20
# 2 component : different donor positions  : Ndonor 10 
# 3 component : nd  numbers of row egeinvectores  100
# 4 component : nd numbers of colums egeinvectores 100



