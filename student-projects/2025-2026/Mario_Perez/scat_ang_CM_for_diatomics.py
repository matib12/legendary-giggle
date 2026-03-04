import numpy as np
import os 
from scipy.constants import c,e,N_A,m_p,h #speed of light, elementary charge, Avogadros const., mass of proton, Plancks const.
import matplotlib.pyplot as plt

os.makedirs('output', exist_ok=True)
# The functionality of this code is to generate 3 graphs one for thi
#
#
#
#
#
#

print('###################################################################################################')
print('Welcome to my program. This program is designed to evaluate the scattering angle for classical mechanics for diatomics. The output is a plot of the scattering angle as a function of the impact parameter as well as a plot of all the trajectories In the next step you will be asked to set the input parameters.')
print('###################################################################################################')
while True:
    input('to continue press enter')
    break
print('###################################################################################################')   

#============================= DEFINING INPUT PARAMETERS =================================================================================================

while True:
    m1 = 19 #float(input('set value for the mass of atom A in amu: '))
    break
while True:
    m2 = 19 # float(input('set value for the mass of atom B in amu: '))
    break
while True:
    ve = 916.93  #float(input('set value for the harmonic frequency in cm-1: '))
    break
while True:
    vexe = 11.32 #float(input('set value for the anharmonicity constant in cm-1: '))
    break
while True:
    re = 1.412 #float(input('set value for the equillibrium distance in Angstrom: '))
    break
while True:
    D0 = 1.644 #float(input('set value for the dissociation energy in eV: '))
    break
while True:
    init_d = float(input('set value for the initial displacement of the atoms in Angstrom: '))
    break
while True:
    Ekin_initial = float(input('set value for the initial kinetic energy in eV: '))
    break
 
print('###################################################################################################')
 
#================================ DEFINING THE FUNCTIONS ==================================================================================================
    
# reduced mass
def m_red(m1,m2):
    return (m1*m2/(m1+m2))*m_p # times mass of proton to convert u to kg
# r
def R(x,y):
    return np.sqrt(x**2+y**2)  

# Morse potential

def MP(r,De,re,beta):
    return (De*e*(1-np.exp(-beta*(1e-10*(r-re))))**2)-De*e #substract De times elementary charge e, because the input of De is given in eV
                                                           #put 1e-10 to convert from Angström to m 

# derivatives of the potential with respect to x and y, using the chain rule ---> dV/dx = dV/dr * dr/dx

def dVdx(r,De,re,beta,x,y):                                    
    return De*e*beta*(np.exp(-beta*(1e-10*(r-re)))-np.exp(-2*beta*(1e-10*(r-re))))*2*x/np.sqrt(x**2+y**2)    #put 1e-10 to convert from Angström to m 
def dVdy(r,De,re,beta,x,y):
    return De*e*beta*(np.exp(-beta*(1e-10*(r-re)))-np.exp(-2*beta*(1e-10*(r-re))))*2*y/np.sqrt(x**2+y**2)    #put 1e-10 to convert from Angström to m 
    
    
#================================ CALCULATING PARAMETERS ====================================================================================================  


# the well depth is calculated by adding the ZPE to the dissociation energy

De = D0 + (c*h*100/e*(0.5*ve - 0.25*vexe)) # c*h*100/e* for converting from cm-1 to eV

# the force constant

k_force = m_red(m1,m2)*(ve*2*np.pi*c*100)**2 # put times 100 to convert cm-1 to m-1

# the morse parameter

beta = np.sqrt(k_force/(2*De*e))

#================ FOR INTEGRATION OF dTHETA/dr ==============================================================================================================

def dthetadr(r,b,V,E):                          # this function needs to be integrated numerically
    return b/r**2 /np.sqrt(1 - b**2/r**2 - V/E)

theta_conv = 0.000001 # convergence value for integral

r_step = 0.0001e-9    # stepsize for integration


theta_val = []        # list for values of theta for specific impact paramter b

#=============================================================================================================================================================


b_val = np.arange(0.1,5.6,0.2) # this creates a list a b values ranging from 0.1 to 5.5 with steps of 0.1

X = []   # lists of x/y values of each cycle for plotting the trajectories at the end
Y = []   #
E = []   # list of energy values of each cycle in order to check if the total energy remains constant

# defining initial conditions 

mr = m_red(m1,m2) # reduced mass
t_step = 1e-17    # time step for integrating newton equations

#========================== THE LOOP OVER THE DIFFERENT VALUES OF b ================================================

for i in range(len(b_val)):
       
    px =  np.sqrt(Ekin_initial*2*mr*e) # initial momentum in x
    py =  0                            # initial momentum in y
                                       
    r_val = [] # list for r values
    xpos = []  # list for x values
    ypos = []  # list for y values
    energ = []
    
    b = b_val[i]                 #selecting the current b value
    x = init_d*1e-10             #converting to m    
    y = b*1e-10                  #converting to m
    r = R(x,y)*1e10              #converting to A
    
    E_prev = (px**2+py**2)/2/mr + MP(r,De,re,beta)  # initial total energy
    
#============================ THE LOOP FOR TRAJECTORY INTEGRATION ====================================================
    
    for k in range(100000): #int(t_lim/t_step)
    
        px = px + dVdx(r,De,re,beta,x,y)*t_step    # current momenta
        py = py + dVdy(r,De,re,beta,x,y)*t_step
        
        x = x - px*t_step/mr                       # current positions
        y = y - py*t_step/mr
         
        r = R(x,y)*1e10                            # current r value
        
        # calculate total energy
        
        H = (px**2+py**2)/2/mr + MP(r,De,re,beta)
        #E_diff = E_prev - H
    
        #E_prev = H
        
        r_val.append(r)  # add all the calculated values to the corresponding lists
        xpos.append(x)
        ypos.append(y)
        energ.append(H)
        
    X.append(xpos) # after completion of the cycles for ONE value of b, add the lists containing the relevant values to the corresponding lists
    Y.append(ypos)
    E.append(energ)
    
    a = min(r_val)            # defining turning point
    
    r_run = a*1e-10           # this radius will be the starting point for the integration of the scattering angle
    
#=============================== INTEGRATION OF SCATTERING ANGLE ==================================================
    
    Ij = []            # integrals of every iteration
    j = 0              # iteration variable 
    int_val = 0        # the initial integral value is set to zero
    
    # for the numerical integration, the trapezoidal rule is used
    
    while True:
        r_run += r_step                                             # the 1e10 / 1e-10 are for conversion between Angstrom and m
        int_val += (dthetadr(r_run-r_step,b*1e-10,MP((r_run-r_step)*1e10,De,re,beta),E_prev)+dthetadr(r_run,b*1e-10,MP(r_run*1e10,De,re,beta),E_prev))/2*(r_run-(r_run-r_step))
        Ij.append(int_val)
        #print(int_val)
        if j > 2:
            if Ij[j] - Ij[j-1] < theta_conv:                   # if convergence is reached, the value of the integral is appended to the list Ij
                theta_val.append((np.pi - 2*Ij[-1])/np.pi*180)
                break
        j += 1  
        
    print('cycles completed: ',i+1,'/',len(b_val))    
    
    
#================================= CREATING THE PLOTS ==============================================================

print('###################################################################################################')

theta_val = np.array(theta_val)  # transform lists to numpy arrays. this is more practical for plotting

X = np.array(X)
Y = np.array(Y)
E = np.array(E)
 
# plot of theta as a function of b

plt.scatter(b_val,theta_val/180,label='calc. values')
plt.title(f'E_kin initial: {Ekin_initial} eV')
plt.hlines(np.min(theta_val)/180, 0, b_val[-1], color='k',linestyle='dotted',label='rainbow angle')
plt.xlabel(r'$\mathrm{b\;/\;\AA}$')
plt.ylabel(r'$\mathrm{\theta}$ / 180°')
plt.legend()
plt.savefig('output/theta_function_b.png')
# plot of the trajectories of each cycle of b

fig,ax = plt.subplots()

for i in range(len(X)):
    ax.plot(X[i],Y[i])
ax.set_title(f'E_kin initial: {Ekin_initial} eV')
ax.set_xlim(-1e-9,1e-9)
ax.set_ylim(-1e-9,1e-9)
ax.set_xlabel('x values / m')
ax.set_ylabel('y values / m')
circle1 = plt.Circle((0, 0), 0.5*re*1e-10, color='b')
ax.add_patch(circle1)
ax.set_aspect('equal', adjustable='box')
plt.savefig('output/Trajectories.png')

# plot of the total energy as a function of time for each cycle of b

time = np.arange(0,len(xpos))

fig2, ax2 = plt.subplots()

for i in range(len(X)):
    ax2.plot(time,E[i]/e)
ax2.set_title(f'E_kin initial: {Ekin_initial} eV')
ax2.set_xlabel('time step')
ax2.set_ylabel('total energy / eV')

plt.savefig('output/energy.png')
