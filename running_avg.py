#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#temp = open("temp.txt",'w+')
import numpy as np
#a = np.zeros([1,1])
#a = np.loadtxt("ener.txt")[:,0]
dataset = np.loadtxt("ener.txt")[:,6]
#print(b)
#dataset = [1,2,3,4,5,6,7,8,4,3,2,4,5,6,3,7,9,8,1,5,3,6,4,4,5,2,4,7]

def runningaverage(values,window):
    weights = np.repeat(1.0, window)/window
    mv_avg = np.convolve(values,weights,'valid')
    return mv_avg

a = runningaverage(dataset,100)
np.savetxt('temp.txt',a,fmt='%8.2f')
#np.savetxt("temp.txt", a, newline="\n")
#print('\n',a,file=temp)
print(a)

