#!/usr/bin/env python
import numpy as np
import math as m
import random
import matplotlib.pyplot as plt
import ising_functions as isf

#Same as ising.py but for 1d periodic lattice
#outputs T, M, E, C, X, S at equilibrium



def plot_line(matrix):	#for plotting 1d spin lattice as black (up) and white (down) dots
	i=0
	n=matrix.shape[0]
	l_up=[]
	l_down=[]
	l_zero1=[]
	l_zero2=[]
	while i<n:
		if matrix[i]==1:
			l_up.append(i)
			l_zero1.append(0)
		if matrix[i]==-1:
			l_down.append(i)
			l_zero2.append(0)
		i=i+1
	plt.plot(l_up,l_zero1, 'o', color='black', markersize=5)
	plt.plot(l_down, l_zero2, 'o', color='white', markersize=5)
	plt.title('1-d Spin Lattice, n=100, h=0')
	plt.show()		
		


n, T=raw_input().split()
n=float(n)
T=float(T)

h, J, k =0.0,1.0,1.0


matrix=isf.line(n) #initialize 1d lattice
mc_steps=n**2 #for n=100 n**3 is too many times



s=0
while s<(mc_steps):
	i=0
	while i<n:

			
		matrix_new=np.copy(matrix)
		matrix_new[i]=-matrix[i]
		E_new=isf.line_energy(matrix_new,h,J)
		E_old=isf.line_energy(matrix,h,J)
		dE=E_new-E_old


		if dE<0:
			matrix=matrix_new


		if dE>0:
			p_flip=m.exp((-1.0/(T*k))*dE)
			rand_num=random.random()

			if rand_num<=p_flip:
				matrix=matrix_new


		i=i+1
	s=s+1
#now close enough to equilibrium to take measurements

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

			
		matrix_new=np.copy(matrix)
		matrix_new[i]=-matrix[i]
		E_new=isf.line_energy(matrix_new,h,J)
		E_old=isf.line_energy(matrix,h,J)
		dE=E_new-E_old

		if dE<0:
			matrix=matrix_new
		

		if dE>0:
			p_flip=m.exp((-1.0/(T*k))*dE)

			rand_num=random.random()

			if rand_num<=p_flip:
				matrix=matrix_new
		

		i=i+1


	M=isf.total_mag(matrix)
	mag_eq=float(M)/(n)
	S=m.log(isf.choose(int(n), int(n+M)//2))
	S_sum=S_sum+S
	E_eq=isf.line_energy(matrix,h,J)/(n)
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
S_max=m.log(isf.choose(int(n), int(n//2)))

xi=(av_mag_sq-av_mag*av_mag)/(T)
C_v=(av_E_sq-av_E*av_E)/(T*T)
	
print T, " ", av_mag, " ", av_E , " ", xi, " ", C_v, " ", av_S/S_max
