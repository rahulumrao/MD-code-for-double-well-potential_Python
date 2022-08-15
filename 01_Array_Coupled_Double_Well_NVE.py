#program velocity verlet for 1D Harmonic oscillator
#filename = input("Enter filename:  ")
import random
import math
import numpy as np
traj = open("traj.xyz","w+")  #open file / 'r+' read and write in the file / 'w+' create file and write 
ener = open("ener.txt","w+")  #open file / 'r+' read and write in the file / 'w+' create file and write 
pote = open("pot.txt", "w+")
#print(filename)
steps = 2 ; n = 10000
k = 1.0 ; T = 10.0
m = 20.0; delT = 1.0 ; U = 0.01
Kb = 1.0 ; a = 1.5 #!Kb = 1.3806*10E-23

vel = np.zeros([steps+1, steps+1])
pos = np.zeros([steps+1, steps+1])
acc = np.zeros([steps+1, steps+1])

pos[0,0] = 1.5 ; pos[0,1] = 0.0
#vel[0,0] = 1.0 ; vel[0,1] = 0.05
#--------------------------------------------------------------------------------
def random_velocities(vel):
	vel_ran1 = random.random()
	vel_ran2 = random.random()
	vel[0,0] = math.sqrt(-2.0*math.log10(vel_ran1)) * 0.5*math.cos(2*math.pi*vel_ran2)
	vel[0,1] = math.sqrt(-2.0*math.log10(vel_ran1)) * 0.5*math.sin(2*math.pi*vel_ran2)
	print(vel_ran1,vel_ran2,vel[0,0],vel[0,1])
#-------------------------------------------------------------------------------------------------------------------
def scale_velocities(vel):

	random_velocities(vel)

	vel[0,0] = vel[0,0]*(math.sqrt(m/(2.0*math.pi*Kb*T)))
	vel[0,0] = vel[0,1]*(math.sqrt(m/(2.0*math.pi*Kb*T)))

	bath_T = (0.5*m*(vel[0,0]**2 + vel[0,1]**2))

	vel[0,0] = vel[0,0]*(math.sqrt((T/bath_T)))
	vel[0,1] = vel[0,1]*(math.sqrt((T/bath_T)))

	bath_T = (0.5 * m*(vel[0,0]**2 + vel[0,1]**2))
	print(vel[0,0],vel[0,1],bath_T)
#--------------------------------------------------------------------------------------
acc[0,0] = -(4*U*pos[0,0]*(pos[0,0]**2 - a**2))/m 
acc[0,1] = -k*pos[0,1]/m
#random_velocities(vel)
scale_velocities(vel)
for i in range (1,n):
	
	pos[1,0] = pos[0,0] + vel[0,0]*delT + (0.5*delT**2)*(acc[0,0])
	pos[1,1] = pos[0,1] + vel[0,1]*delT + (0.5*delT**2)*(acc[0,1])
    
	acc[1,0] = -(4*U*pos[1,0]*(pos[1,0]**2 - a**2))/m 
	acc[1,1] = -k*pos[1,1]/m
    
	vel[1,0] = vel[0,0] + (0.5*(acc[0,0] + acc[1,0])*delT)	
	vel[1,1] = vel[0,1] + (0.5*(acc[0,1] + acc[1,1])*delT)
    
	t = delT*i
	pos[0,0] = pos[1,0] ; pos[0,1] = pos[1,1]
	vel[0,0] = vel[1,0] ; vel[0,1] = vel[1,1]
	acc[0,0] = acc[1,0] ; acc[0,1] = acc[1,1]
	KE = 0.5*(m*(vel[0,0]**2 + vel[0,1]**2))
	PE = U*((pos[0,0]**2 - a**2)**2) + 0.5*k*pos[0,1]**2 
	Temp = (2.0/3.0)*(KE)/Kb
	
	print('1',"\n",file=traj)
	print("M","\t",'%.6f'%pos[0,0],"\t",'%.6f'%pos[0,1],"\t",'0',"\t",'%.6f'%vel[0,0],file=traj)
	print(i,"\t",'%.6f'%pos[0,0],"\t",'%.6f'%pos[0,1],"\t",'%.6f'%PE,file=pote)
	print('%.2f'%t,"\t",'%.6f'%KE,"\t",'%.6f'%PE,"\t",'%.6f'%(KE+PE),"\t",'%.4f'%Temp,file=ener)
print('\n',"Temp",Temp)
traj.close() ; ener.close()
