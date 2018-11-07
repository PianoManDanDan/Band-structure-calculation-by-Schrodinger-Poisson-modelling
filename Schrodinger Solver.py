# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 19:53:53 2015

@author: Chris
"""

#==============================================================================
# Imports
#==============================================================================
from pylab import *


close('all')


#~~~~~~~~~~~~~~~~~~~~~~~~~Constants/Set up~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

global q, eps
mev=1.60217656e-22; hbar=1.0545718e-34; q=1.60217656e-19; kb=1.38e-23; T=300.; eps0=8.85e-12; m_e=9.10938356e-31; Ang=1.0e-10

#==============================================================================
# material properties
#==============================================================================

me_ = 0.014*m_e
mh_ = 0.43*m_e
Eg = 0.17 #eV
Eg = 0.23 #eV



#==============================================================================
# Schrodinger paramaters
#==============================================================================


n = 5 #Number of states to find
nstep = 0.1 #Energy step for searching (mev)


# x
N = 3000 #Number x points
x = linspace(-200.,200.,N)*1e-9 # metres
dx = abs(x[0]-x[1])

#Potential - must be array like x, minimum of array at 0

#Square well, 100nm wide, 100meV high
Pot = ones_like(x)*200
Pot[abs(x)<50e-9] = 0
Pot*=mev


#Assymetric well
Pot = ones_like(x)*200
Pot[abs(x)<40e-9] = 0
Pot[(x <40e-9) & (x>0e-9)] = 20
Pot*=mev


#Well inside a well
Pot = ones_like(x)*200
Pot[abs(x)<40e-9] = 20
Pot[abs(x)<20e-9] = 0

Pot*=mev



#~~~~~~~~~~~~~~~~~~~~~~~~~Schrodinger Solver~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def Nicer_figures():
    label_size = 15.
    rcParams["font.size"] = label_size
    rcParams['xtick.labelsize'] = label_size
    rcParams['ytick.labelsize'] = label_size
    rcParams['axes.titlesize'] = label_size
    rcParams['axes.labelsize'] = label_size
#
    rcParams['axes.linewidth'] = 2.0
    rcParams['xtick.major.size'] = 7.0
    rcParams['xtick.minor.size'] = 4.0
    rcParams['xtick.major.width'] = 2.0
    rcParams['ytick.major.size'] = 7.0
    rcParams['ytick.major.width'] = 2.0
    rcParams['ytick.minor.size'] = 4.0
    rcParams['ytick.minor.width'] = 2.0
    rcParams['lines.linewidth'] = 2.0
    rcParams['legend.fontsize'] = label_size
    rcParams['legend.frameon'] = True

Nicer_figures()




#~~~~~~~~~~~~~~~~~~~~~~~~~Schrodinger Solver~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def Shooting_Method(Etemp,V,m,dz):                                              # Creates a phi profile using V array, constants E and m
                                                                                # Variable effective mass
    phi=zeros(len(V))                                                           # Empty Bulk phi array
    phi[0]=0.                                                                    # Standard Conditions
    phi[1]=1.
    alpha=(2*(dz/float(hbar))**2)
    mplus=(m[2:]+m[1:-1])/2.                                                # First forward average step of mass profile
    mminus=(m[:-2]+m[1:-1])/2.
    aa = alpha*(V[1:-1]-Etemp)                                                # First backward average step of mass profile
    for i in range(len(phi)-2):                                                 # For loop over the range of interior points                                                # First backward average step of mass profile
        phi[2+i]=(((aa[i])+(1./mplus[i])+(1./mminus[i]))*phi[i+1]-(phi[i]/mminus[i]))*mplus[i] # Varible mass shooting method
    return phi, phi[-1]

def Newton_Raphson_Method(phi0,phi1,phi2,En,Estep):

        difphi=(phi2-phi0)/(2*Estep) #0.05 dE
        Enew=En-(phi1/difphi)                   # Standard NR equation
        return Enew                             # Returns Confined Energy

def Schrodinger_Solver(n,nstep, x, Pot):                          # n is the numebr of states to find
    Phinew     =[]
    Phinorm    =[]
    Xstr       =[]
    dZstr      =[]
    Potstr     =[]
    E_tempstr  =[]
    Points = len(Pot)
    dZ = abs(x[1]-x[0])
    E_start=0.005*mev
    E_step =0.005*mev
    M_eff=ones_like(Pot)* me_

    phistates=[]
    Eftemp=0*mev
    _, phi_initial_End=Shooting_Method(0,Pot,M_eff,dZ)
    _  , phi_final_End=Shooting_Method(Eftemp,Pot,M_eff,dZ)
    Estate=[]
    while len(Estate)<n:
        phi_initial_End = phi_final_End
        Eftemp=Eftemp+nstep*mev
        _, phi_final_End=Shooting_Method(Eftemp,Pot,M_eff,dZ)
        if round(Eftemp/mev - nstep,1)%10 == 0:
            print ("still searching in well, reached: ", Eftemp/mev-nstep, 'mev')

        if phi_final_End*phi_initial_End <0.0:
            Estate.append(Eftemp-nstep*mev)
            print ("\033[31 ; 1 ; 43m finding states " + str(Eftemp/mev-nstep) + "meV\033[0m")
        if Eftemp > max(Pot):
            print ('nope nope nope')
            break
    print ("\033[31 ; 2m Energy states: \033[0m \033[32 ; 2m" + "%.2f meV, "*len(Estate)%(tuple(array(Estate)/mev)) + "\033[0m")

    for j in range(n):
#        E_temp=E_start
        E_temp=Estate[j]

        _, phi_initial_End=Shooting_Method(0    ,Pot,M_eff,dZ)             # starting shooting method only for last point
        _  , phi_final_End  =Shooting_Method(E_temp,Pot,M_eff,dZ)             # final shooting method only for last point

        while phi_final_End*phi_initial_End < 0.0:                                  # Checking to see if the confined energy has been past(finding rough confined energy)
            phi_initial_End = phi_final_End
            _  , phi_final_End=Shooting_Method(E_temp+E_step,Pot,M_eff,dZ)         # running that energy through shooting method                          # going up an energy step to run through shooting method

        phinew=zeros(Points)
        phinew_End=1.                                                 # Empty phi profile
        i=0

        while abs(phinew_End)>0.0005:                                                # Setting the accuracy for the final point
           E2=E_temp+E_step                                                          # A bit above guess
           E0=E_temp-E_step                                                          # A bit below guess
           _, phi1_End=Shooting_Method(E_temp,Pot,M_eff,dZ)
           _, phi2_End=Shooting_Method(E2,Pot,M_eff,dZ)
           _, phi0_End=Shooting_Method(E0,Pot,M_eff,dZ)
           Enew=Newton_Raphson_Method(phi0_End,phi1_End,phi2_End,E_temp,E_step)
           phinew, phinew_End=Shooting_Method(Enew,Pot,M_eff,dZ)
           E_temp=Enew
           if i ==50:
               break
           i+=1


        inti     =zeros(Points)
        normedphi=zeros(Points)
        global index1
        index1 = where(phinew[-len(phinew)//3:]>0.01*min(phinew))[0][0]
        index1 += len(phinew)-len(phinew)//3
        inti=(phinew[:index1]*phinew[:index1])*dZ
        const=(1./(sum(inti)))**0.5

        normedphi=(const*phinew)/1000. +E_temp/mev
        normedphi[index1:] = normedphi[index1]
        Phinorm  .append(normedphi)
        Phinew   .append(phinew)
        dZstr    .append(dZ)
        E_tempstr.append(E_temp)
        print ("\033[32 ; 2m Found energy state " +str(n) + " at %.2f mev"%(E_temp/mev) + "\033[0m")
    return Phinew,Phinorm,dZstr,E_tempstr





#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
# # # # Find solutions
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================


phinew1,Phinorm1,dz1,Etemp1=Schrodinger_Solver(n, nstep, x, Pot)

figure(0)
clf()
figManager = matplotlib.pyplot.get_current_fig_manager()
# figManager.window.showMaximized()

for p in range(n):
    plot(x*1e9,Phinorm1[p], label = '$E = %.2f\ meV$'%(Etemp1[p]/mev))

plot(x*1e9,Pot/mev, 'k')

leg = legend(loc = 0)
leg.set_draggable(True)
xlim(min(x)*1e9,max(x)*1e9)
ylim(min(Pot)*0.9/mev-5, max(Pot)*1.1/mev+5)

xlabel('$x\ (nm)$')
ylabel('$V\ (meV)$')
plt.show()
