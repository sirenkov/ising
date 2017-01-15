#!/usr/bin/env python
import numpy as np
import math as m
import random
import matplotlib.pyplot as plt
import ising_functions as isf

#FOR 2-D SQUARE LATTICE
#code implementing Metropolis algorithm
#user inputs size (n) and temperature (T)
#code outputs T, M, E, X, C, S at equilibrium
#we have to run export PYTHONPATH="${PYTHONPATH}:/home/vlad/Documents/Python/modules" 
#on terminal before we can use ising_functions as intended

#FOR TRIANGLURAR LATTICE, REPLACE isf.energy(matrix, h, J) with isf.hex_energy(matrix, h, J)
#FOR NON-PERIODIC BC, REPLACE isf.energy(matrix, h, J) with isf.np_energy(matrix, h, J)
#FOR NON-RANDOM INITIAL CONFIGURATION, REPLACE isf.lattice(n) with isf.c_lattice(n,1) or chequer_lattice(n) or half_lattice(n)

#print "input n (size) and T (temperature), space separated"
n, T=raw_input().split()
n=float(n)
T=float(T)

h, J, k =1.0,-1.0,1.0


matrix=isf.lattice(n) #initializes 2d lattice
mc_steps=n**3 #code will run throught the entire lattice 'mc_steps' times

s=0
while s<(mc_steps):
	i=0
	while i<n:
		j=0
		while j<n:
			
			#now implement metropolis algorithm
			matrix_new=np.copy(matrix)
			matrix_new[i,j]=-matrix[i,j]
			E_new=isf.energy(matrix_new,h,J)
			E_old=isf.energy(matrix,h,J)
			dE=E_new-E_old
	
			if dE<0:
				matrix=matrix_new
			
			if dE>0:
				p_flip=m.exp((-1.0/(T*k))*dE)
				rand_num=random.random() #random numbers for flip probability
				if rand_num<=p_flip:
					matrix=matrix_new
			j=j+1
			
		i=i+1
	s=s+1


#now we are 'close enough' to equilibrium and we can start taking measurements

s=0
sum_mag=0
sum_mag_sq=0
sum_E=0
sum_E_sq=0
S_sum=0
collect=mc_steps*2 

while s<(collect):
	i=0
	while i<n:
		j=0
		while j<n:
			
			matrix_new=np.copy(matrix)
			matrix_new[i,j]=-matrix[i,j]
			E_new=isf.energy(matrix_new,h,J)
			E_old=isf.energy(matrix,h,J)
			dE=E_new-E_old
	
			if dE<0:
				matrix=matrix_new
				

			if dE>0:
				p_flip=m.exp((-1.0/(T*k))*dE)
				
				rand_num=random.random()
				if rand_num<=p_flip:
					matrix=matrix_new
					
			j=j+1
		i=i+1

	
	#collect values of M, S, etc. for averaging
	M=isf.total_mag(matrix)
	mag_eq=float(M)/(n*n)
	S=m.log(isf.choose(int(n*n), int(n*n+M)//2))
	S_sum=S_sum+S
	E_eq=isf.energy(matrix,h,J)/(n*n)
	sum_mag=sum_mag+mag_eq
	sum_mag_sq=sum_mag_sq+mag_eq*mag_eq
	sum_E=sum_E+E_eq
	sum_E_sq=sum_E_sq+E_eq*E_eq

	s=s+1

#now take averages
av_mag=sum_mag/(collect)
av_mag_sq=sum_mag_sq/(collect)
av_E=sum_E/(collect)
av_E_sq=sum_E_sq/(collect)
av_S=S_sum/(collect)
S_max=m.log(isf.choose(int(n*n), int(n*n//2)))

xi=(av_mag_sq-av_mag*av_mag)/(T)
C_v=(av_E_sq-av_E*av_E)/(T*T)
	
print T, " ", av_mag, " ", av_E , " ", xi, " ", C_v, " ", av_S/S_max
