#!/usr/bin/env python
import numpy as np
import math as m
import random
import matplotlib.pyplot as plt
import ising_functions as isf

#Same as ising.py but for 3d Simple Cubic lattice
#outputs T, M, E, C, X, S at equilibrium
#FOR BCC LATTICE, REPLACE isf.cube_energy(matrix, h, J) with isf.bcc_energy(matrix, h, J)
#FOR FCC LATTICE, REPLACE isf.cube_energy(matrix, h, J) with isf.fcc_energy(matrix, h, J)

#print "input n (size) and T (temperature), space separated"
n, T=raw_input().split()
n=float(n)
T=float(T)

h, J, k =0.0,1.0,1.0


matrix=isf.cubic(n)	#initialize cubic lattice

mc_steps=n**3	#code will run through lattice mc_steps//n times

s=0

while s<(mc_steps//n):
	kk=0
	while kk<n:
		i=0
		while i<n:
			j=0
			while j<n:
			
				matrix_new=np.copy(matrix)
				matrix_new[kk,i,j]=-matrix[kk,i,j]
				E_new=isf.cube_energy(matrix_new,h,J)
				E_old=isf.cube_energy(matrix,h,J)
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
		kk=kk+1
	s=s+1

#now we are close enough to equilibrium	
s=0
sum_mag=0
sum_mag_sq=0
sum_E=0
sum_E_sq=0
S_sum=0
collect=mc_steps

#start collecting data

while s<(collect):
	kk=0
	while kk<n:
		i=0
		while i<n:
			j=0
			while j<n:
			
				matrix_new=np.copy(matrix)
				matrix_new[kk,i,j]=-matrix[kk,i,j]
				E_new=isf.cube_energy(matrix_new,h,J)
				E_old=isf.cube_energy(matrix,h,J)
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
		kk=kk+1

	#data for averaging
	M=isf.total_mag(matrix)
	mag_eq=float(M)/(n*n*n)
	S=m.log(isf.choose(int(n*n*n), int(n*n*n+M)//2))	
	S_sum=S_sum+S
	E_eq=isf.cube_energy(matrix,h,J)/(n*n*n)
	sum_mag=sum_mag+mag_eq
	sum_mag_sq=sum_mag_sq+mag_eq*mag_eq
	sum_E=sum_E+E_eq
	sum_E_sq=sum_E_sq+E_eq*E_eq

	s=s+1

av_mag=sum_mag/(collect)
av_mag_sq=sum_mag_sq/(collect)
av_E=sum_E/(collect)
av_E_sq=sum_E_sq/(collect)
av_S=S_sum/(collect)
S_max=m.log(isf.choose(int(n*n*n), int(n*n*n//2)))

xi=(av_mag_sq-av_mag*av_mag)/(T)
C_v=(av_E_sq-av_E*av_E)/(T*T)
	
print T, " ", av_mag, " ", av_E , " ", xi, " ", C_v, " ", av_S/S_max
