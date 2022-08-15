#program velocity verlet for 1d harmonic oscillator
#filename = input("Enter filename:  ")
import random
import math
traj = open("traj.xyz","w+")  #open file / 'r+' read and write in the file / 'w+' create file and write 
ener = open("ener.txt","w+")  #open file / 'r+' read and write in the file / 'w+' create file and write 
pote = open("pot.txt", "w+")
#print(filename)
#-------------------------------------------------------------------------------------------------------------------
k = 1.0 ; n = 50000
m = 20.0; delT = 0.01 ; U = 0.01
x = 1.5 ; y = 0.0 ; T =10.0
vx = 1.0 ; vy = 0.5
Kb = 1.0 ; a = 1.5 #!Kb = 1.3806*10E-23
#-------------------------------------------------------------------------------------------------------------------
A1 = -(4*U*x*(x**2 - a**2))/m 
A2 = -k*y/m
for i in range (1,n):
	
	x1 = x + vx*delT + (0.5*delT**2)*(A1)
	y1 = y + vy*delT + (0.5*delT**2)*(A2)
    
	A11 = -(4*U*x1*(x1**2 - a**2))/m 
	A22 = -k*y1/m
    
	v1 = vx + (0.5*(A11 + A1)*delT)	
	v2 = vy + (0.5*(A22 + A2)*delT)
    
	t = delT*i
	x = x1 ; y = y1
	vx = v1 ; vy = v2
	A1 = A11 ; A2 = A22
	KE = 0.5*(m*(vx**2 + vy**2))
	PE = U*((x**2 - a**2)**2)+ 0.5*k*y**2 
	Temp = (2.0/3.0)*(KE)/Kb
#-------------------------------------------------------------------------------------------------------------------
	print('1',"\n",file=traj)
	print("M","\t",'%.6f'%x,"\t",'%.6f'%y,"\t",'0',"\t",'%.6f'%vx,file=traj)
	print(i,"\t",'%.6f'%x,"\t",'%.6f'%y,"\t",'%.6f'%PE,file=pote)
	print('%.2f'%t,"\t",'%.6f'%KE,"\t",'%.6f'%PE,"\t",'%.6f'%(KE+PE),"\t",'%.4f'%Temp,file=ener)
print('\n',"Temp = ",Temp,'\n')
traj.close() ; ener.close()
