#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------------------------------------------------------------------------
''' THIS PROGRAM IS WRITTEN IN PYTHON3 TWO SOLVE THE PROBLEM OF DOUBLE WELL POTENTIAL WHICH IS HARMONIC
IN y-axis AND HAS A BRRIER OF 0.01au IN x-axis WITH THE TWO MINIMUM AT +1.5 AND -1.5
{U = 0.01*(x**2 - a**2)**2 - 0.5*k*y**2}
REQUIRED LIBRARY FOR THIS POROGRAM IS "numpy, math and matplotlib(IF YOU WANT TO PLOT USING PYTHON).
%% WRITTEN IN APRIL 3rd 2019 in SL302 by RAHUL VERMA %% '''
#----------------------------------------------------------------------------------------------------------------------------------------------
import os
import random
import math as m
import numpy as np
#import matplotlib.pyplot as plt
#filename = input("Enter filename:  ") #ask user to input filename
traj = open("traj.xyz","w+")  #open file / 'r+' read and write in the file / 'w+' create file and write 
ener = open("ener.txt","w+")  #open file / 'r+' read and write in the file / 'w+' create file and write 
pote = open("pot.txt", "w+")
#noise = open("noise.txt","w+")
#temp = open("Temp.txt","w+")
phase = open("ps.txt","w+")
#----------------------------------------------------------------------------------------------------------------------------------------------
# %% SCALING THE INITIAL VELOCITIES
#----------------------------------------------------------------------------------------------------------------------------------------------
''' THIS PART OF PROGRAM GENERATES THE RANDOM VELOCITIES FOR PARTICLE AND 
THEN SCALES ACCORDING TO THE BOLTZMANN DISTRIBUTION FOR A GIVEN TEMPERATURE
'''
def random_velocity():
    vel_ran1 = random.random() ; vel_ran2 = random.random()
    vel = np.zeros([2,2])
    vel[0,0] = m.sqrt(-2.0*m.log10(vel_ran1)) * 0.5*m.cos(2*m.pi*vel_ran2)
    vel[0,1] = m.sqrt(-2.0*m.log10(vel_ran1)) * 0.5*m.sin(2*m.pi*vel_ran2)
    #print('\n',vel_ran1,vel_ran2,vel[0,0],vel[0,1])
    return vel
#----------------------------------------------------------------------------------------------------------------------------------------------
def scale_velocity(m1, Kb, T):
    vel = random_velocity()
    vel[0,0] = vel[0,0]*(m.sqrt(m1/(2.0*m.pi*Kb*T)))
    vel[0,1] = vel[0,1]*(m.sqrt(m1/(2.0*m.pi*Kb*T)))

    bath_T = (0.5*m1*(vel[0,0]**2 + vel[0,1]**2))/Kb
        
    vel[0,0] = vel[0,0]*(m.sqrt((T/bath_T)))
    vel[0,1] = vel[0,1]*(m.sqrt((T/bath_T)))

    bath_T = (0.5 * m1*(vel[0,0]**2 + vel[0,1]**2))/Kb
    print('\n',vel[0,0],vel[0,1],bath_T,'\n')
    return vel
#----------------------------------------------------------------------------------------------------------------------------------------------
def plot(file):
    import matplotlib.pyplot as plt
    import numpy as np
    a = np.loadtxt(file, unpack=True, skiprows=2)
    plt.plot(a[0], a[7], label = "Langevin Thermostat (NVT)");
    #plt.plot(a[0],a[1], label = "Langevin Thermostat (NVT)");
    plt.xlabel('MD Steps')
    plt.ylabel('Temperature')
    plt.legend();
    plt.show();
#----------------------------------------------------------------------------------------------------------------------------------------------
def running_average(values,window):
    weights = np.repeat(1.0, window)/window
    mv_avg = np.convolve(values,weights,'valid')
    return mv_avg
#----------------------------------------------------------------------------------------------------------------------------------------------
# %% FUNCTIONS FOR THE POSITION, VELOCITY AND FORCES
#----------------------------------------------------------------------------------------------------------------------------------------------
def position (pos,vel,acc,delT):
    pos[0,0] = pos[0,0] + vel[0,0]*delT + (0.5*delT**2)*(acc[0,0])
    pos[0,1] = pos[0,1] + vel[0,1]*delT + (0.5*delT**2)*(acc[0,1])
    return pos
#----------------------------------------------------------------------------------------------------------------------------------------------
def force(pos,m1,U,k,a):
    acc = np.zeros([2,2])
    acc[0,0] = -(4*U*pos[0,0]*(pos[0,0]**2 - a**2))/m1
    acc[0,1] = -k*pos[0,1]/m1
    return acc
#----------------------------------------------------------------------------------------------------------------------------------------------
def velocity (acc,vel,delT):
    vel[0,0] = vel[0,0] + 0.5*acc[0,0]*delT
    vel[0,1] = vel[0,1] + 0.5*acc[0,1]*delT
    return vel
#----------------------------------------------------------------------------------------------------------------------------------------------
def langevin_force (pos,delT,gamma,acc,vel,U,m1,k,a,Kb,T):
   ''' LANGEVIN MODIFIED FORCES TO CONTROL THE SYSTEM TEMPERATURE WHERE,
   THE "gamma" IS THE FRICTION CO-EFFICIENT. THIS CONSISTS OF ADDING A RANDOM FORCE
   AND SUBTRACTING A FRICTION FORCE FROM EACH ATOM DURING EACH INTEGRATION STEP.
   THE RANDOM FORCE IS CALCULATED SUCH THAT THE AVERAGE FORCE IS ZERO.
   '''
   acc = np.zeros([2,2])
   
   r1 = random.random() ;   r2 = random.random()

   r11 = m.sqrt(-2.0*m.log10(r1))*m.cos(2.0*m.pi*r2)
   r22 = m.sqrt(-2.0*m.log10(r1))*m.sin(2.0*m.pi*r2)

   sigma = m.sqrt((2.31*Kb*T*gamma*m1)/delT)
   noise1 = sigma*r11 ; noise2 = sigma*r22
    
   fx = -(4.0*U*pos[0,0]*(pos[0,0]**2 - a**2)) - gamma*vel[0,0]*m1 + noise1
   fy = -k*pos[0,1] - gamma*vel[0,1]*m1 + noise2
    
   acc[0,0] = fx/m1 ; acc[0,1] = fy/m1
   #noise.write('{0}\t{1}\t{2}\n\t{3}\t{4}'.format(i,'%.8f'%acc[0,0],'%.8f'%acc[0,1],A1,A2))
   return acc
#----------------------------------------------------------------------------------------------------------------------------------------------
# %% VELOCITY VERLET INTEGRATION MODULE
#----------------------------------------------------------------------------------------------------------------------------------------------
'''
MAIN PART OF THE PROGRAM TO RUN VELOCITY VERLET INTEGRATION 
MODULE TO SOLVE THE NEWTON's EQUATION OF MOTION
Calculate velocity  ====>
Calculate position  ====>
Derive forces       ====>
Calculate velocity  ====>
'''
def main():
#    condition = input("\n CONSTANT TEMPERATURE LANGEVIN DYNAMICS ?? (Y/N):   ")
#    if (condition in ('Y','y','yes')):   gamma = 1e-2
#    elif (condition in ('N','n','no')): gamma = 0.0
#    else :
#        sys.exit(print("\n PLEASE CONSISTENT WITH THE OPTIONS... THANK YOU ! \n")) 

    steps = 1 ; n = 500000
    a = 1.5 ; U = 0.01 ; k = 1.0 ; m1 = 20000.0 ; T = 300.0 ; Kb = 3.15e-6 ; delT =  100.0 ; gamma = 0.002
    ''' INITIAL POSTION DEFINED AT ONE MINIMUM (1.5AU) BARRIER BETWEEN TWO MINIMUM
    IN x-axis IS 0.01 AU AND SPIRNG CONSTANT 1.0 AU MASS OF THE PARTICLE IS 20000amu 
    Time Step IS 50 AU AND FRICTION COEFFICIENT IS 0.001 WITH THE UNIT OF time^-1
    '''
    vel = np.zeros([steps+1, steps+1])
    pos = np.zeros([steps+1, steps+1])
    acc = np.zeros([steps+1, steps+1])

    pos[0,0] = -1.5 ; pos[0,1] = 0.0
    vel = scale_velocity(m1, Kb, T)
    acc = langevin_force(pos, delT, gamma, acc, vel, U, m1, k, a, Kb, T)
    for i in range (1,n):
        pos = position(pos, vel, acc, delT)
        vel = velocity(acc,vel,delT)
#        acc = force(pos,m1,U,k,a)
        acc = langevin_force(pos, delT, gamma, acc, vel, U, m1, k, a, Kb, T)
        vel = velocity(acc,vel,delT)
#=============================================================================================================================================
        t = delT*i
        KE = 0.5*m1*(vel[0,0]**2 + vel[0,1]**2)
        PE = U*((pos[0,0]**2 - a**2)**2) + 0.5*k*pos[0,1]**2
        Temp = 2*KE/Kb
        #count += i ; T += Temp
        #Te = count/T
        #print(count,T,Te)
        if (i%10 == 0 ):
            print('1','\n',file=traj)
            print("M \t",'%.6f'%pos[0,0,],'\t %.6f'%pos[0,1],'\t 0',file=traj)
            print(t,'\t %.6f'%pos[0,0],'\t %.6f'%pos[0,1],'\t %.6f'%PE,file=pote)
            print(t,'\t %.6f'%pos[0,0],'\t %.6f'%pos[0,1],'\t %.6f'%vel[0,0],'\t %.6f'%vel[0,1],file=phase)
            ener.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(i,'%.6f'%vel[0,0],'%.6f'%vel[0,1],'%.6f'%KE,'%.6f'%PE,'%.6f'%(KE+PE),'%.4f'%Temp))
#=============================================================================================================================================
    #dataset = np.loadtxt("ener.txt")[:,6]
    #T = running_average(dataset,200)
    #np.savetxt('temp.txt',T,fmt='%5.2f')
    #m = int(n/15.9)
    #print("\n JOB DONE \t Temp =",T[m:m+1])
    print("JOB DONE \n")
    os.system('./average_temp.x') 
    plot('ENERGIES.txt')
    traj.close(); ener.close(); pote.close(); phase.close()
    return None
#----------------------------------------------------------------------------------------------------------------------------------------------

if __name__=="__main__":
    main()
#----------------------------------------------------------------------------------------------------------------------------------------------
