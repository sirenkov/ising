#!/usr/bin/env python
import numpy as np
import math as m
import random
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.patches as mpatches

#module containing all important functions for Ising Model:

def choose(n, k): #choose function for shannon entropy calculation
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in xrange(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0


def line(n):	#initializes 1-d lattice of size n of +1 and -1's, randomly assignes
	matrix=np.random.randint(2, size=(n)) #make random array of 1's and 0's
	matrix[matrix<1]=-1 #replaces 0's with -1's
	return matrix

def lattice(n): #initializes n*n lattice of 1 and -1, randomly assigned
	matrix=np.random.randint(2, size=(n, n)) #make random array of 1's and 0's
	matrix[matrix<1]=-1 #replaces 0's with -1's
	return matrix

def cubic(n):	#initializes n*n*n lattice of 1 and -1, randomly assigned
	matrix=np.random.randint(2, size=(n, n, n)) #make random array of 1's and 0's
	matrix[matrix<1]=-1 #replaces 0's with -1's
	return matrix

def c_lattice(n,c): #initializes n*n lattice of numbers "c"
	matrix=np.zeros((n,n))
	matrix[matrix<1]=c #replaces 0's with -1's
	return matrix

def chequer_lattice(n): #initializes spin lattice with + and -1 in chequer-board pattern
	matrix=np.ones((n,n))
	i=0
	while i<=(n-1):
		j=0
		while j<=(n-1):
			if (j%2==0 and i%2==0) or (j%2!=0 and i%2!=0):
				matrix[i,j]=-matrix[i,j]
			j=j+1
		i=i+1
	return matrix

def half_lattice(n): #intializes lattice that is half +1, half -1
	matrix=np.ones((n,n))
	half=n//2
	i=0
	while i<=(n-1):
		j=0
		while j<=(n-1):
			if (j>=half):
				matrix[i,j]=-matrix[i,j]
			j=j+1
		i=i+1
	return matrix


def periodic(matrix): #adds extra rows & columns identical to first & last
	n=matrix.shape[0]
	
	matrix1=np.concatenate(((np.array([matrix[:,-1]])).T,matrix,(np.array([matrix[:,0]])).T),axis=1)
	matrix2=np.concatenate(((np.array([matrix1[-1,:]])),matrix1,(np.array([matrix1[0,:]]))),axis=0)
	matrix2[-1,-1]=matrix[0,0]
	matrix2[0,-1]=matrix[-1,0]
	matrix2[-1,0]=matrix[0,-1]
	matrix2[0,0]=matrix[-1,-1]

	return matrix2

def total_mag(matrix): #calculates total magnetisation in lattice
	return np.sum(matrix)

def energy(matrix,h,J): #calculate energy of specific 2d square matrix
	matrix=periodic(matrix)
	n=matrix.shape[0]
	h_term=0
	j_term=0
	i=1
	while i<(n-1):
		j=1
		while j<(n-1):
			#nearest neighbour interaction term:
			j_term=j_term+matrix[i,j]*matrix[i+1,j]+matrix[i,j]*matrix[i,j+1]

			#external field dependency term
			h_term=h_term+matrix[i,j]
			j=j+1
		i=i+1
	
	E=-h*h_term-J*j_term
	return E


def np_energy(matrix,h,J): #calculate energy of 2d matrix, NON PERIODIC
	
	n=matrix.shape[0]
	h_term=np.sum(matrix)
	j_term=0
	i=0
	while i<=(n-1):
		j=0
		while j<=(n-1):
			if i!=(n-1) and j!=(n-1):
				j_term=j_term+matrix[i,j]*matrix[i+1,j]+matrix[i,j]*matrix[i,j+1]
			if i==(n-1) and j!=(n-1):#we are in the last row
				j_term=j_term+matrix[i,j]*matrix[i,j+1]
			if i!=(n-1) and j==(n-1):#we are in the last column
				j_term=j_term+matrix[i,j]*matrix[i+1,j]
			if i==(n-1) and j==(n-1):#last entry in matrix has been counted
				j_term=j_term
			j=j+1
		i=i+1
	
	E=-h*h_term-J*j_term
	
	return E



def hex_energy(matrix,h,J): #calculate energy of triangular lattice
	
	matrix=periodic(matrix)
	n=matrix.shape[0]
	h_term=0
	j_term=0
	i=1
	while i<(n-1):
		j=1
		while j<(n-1):
			j_term=j_term+matrix[i,j]*matrix[i+1,j]+matrix[i,j]*matrix[i,j+1]+matrix[i,j]*matrix[i+1,j+1]
			h_term=h_term+matrix[i,j]
			j=j+1
		i=i+1
	E=-h*h_term-J*j_term
	return E


def line_energy(matrix,h,J):	#energy of 1d matrix
	n=matrix.shape[0]

	matrix1=np.append(matrix,matrix[0])

	j_term=0
	h_term=0
	i=0
	while i<n:
		h_term=h_term+matrix1[i]
		j_term=j_term+matrix1[i]*matrix1[i+1]
		i=i+1

	return (-h*h_term-J*j_term)


def cube_energy(cube,h,J):	#energy of simple cubic matrix


	cube1=np.concatenate(((np.array([cube[-1,:]])),cube,(np.array([cube[0,:]]))),axis=0) #for periodic boundary

	k=1
	n=cube1.shape[0]
	h_term=0
	j_term=0

	while k<(n-1):
		matrix=periodic(cube1[k,:])
		i=1
		while i<(n-1):
			j=1
			while j<(n-1):
				j_term=j_term+matrix[i,j]*matrix[i+1,j]+matrix[i,j]*matrix[i,j+1]+matrix[i,j]*cube1[k+1][i-1,j-1]
				h_term=h_term+matrix[i,j]
				j=j+1
			i=i+1
		k=k+1

	E=(-h*h_term-J*j_term)
	return E

def bcc_energy(cube,h,J):	#energy of bcc lattice

	cube1=np.concatenate(((np.array([cube[-1,:]])),cube,(np.array([cube[0,:]]))),axis=0)	#for periodic boundary

	k=1
	n=cube1.shape[0]
	h_term=0
	j_term=0

	while k<(n-1):
		matrix=periodic(cube1[k,:])
		
		i=1
		while i<(n-1):
			j=1
			while j<(n-1):
				j_term=j_term+matrix[i,j]*cube1[k+1][i-1,j-1]+matrix[i,j+1]*cube1[k+1][i-1,j-1]+matrix[i+1,j]*cube1[k+1][i-1,j-1]+matrix[i+1,j+1]*cube1[k+1][i-1,j-1]
				h_term=h_term+matrix[i,j]
				j=j+1
			i=i+1

		k=k+1
	E=(-h*h_term-J*j_term)
	return E


def fcc_energy(cube,h,J):	#energy of fcc lattice

	cube1=np.concatenate(((np.array([cube[-1,:]])),cube,(np.array([cube[0,:]]))),axis=0)	#for periodic boundary

	k=1
	n=cube1.shape[0]
	h_term=0
	j_term1=0
	j_term2=0
	j_term3=0
	while k<(n-1):
		matrix=periodic(cube1[k,:])
		matrixb=periodic(cube1[k+1,:])
		matrixc=periodic(cube1[k-1,:])
		i=1
		while i<(n-1):
			j=1
			while j<(n-1):
				j_term1=j_term1+matrix[i,j]*matrix[i+1,j]+matrix[i,j]*matrix[i-1,j+1]
				j_term2=j_term2+matrix[i,j]*matrixb[i,j]+matrix[i,j]*matrixb[i+1,j-1]
				j_term3=j_term3+matrix[i,j]*matrixc[i,j+1]+matrix[i,j]*matrixc[i+1,j]
				h_term=h_term+matrix[i,j]
				j=j+1
			i=i+1

		k=k+1
	E=-h*h_term-J*(j_term1+j_term2+j_term3)
	return E

def plot_lattice(matrix, n, T, h, J):	#plots matrix as black (spin up) and white (spin down) square lattice
	n=matrix.shape[0]
	bounds4=[0, 0.00001, n*n]
	cmap4 = colors.ListedColormap(['white','black'])
	norm4 = colors.BoundaryNorm(bounds4, cmap4.N)
	img = plt.imshow(matrix,interpolation='none', cmap = cmap4, norm=norm4)
	patch4 = mpatches.Patch(color='black', label='spin up')
	plt.title('Spin Lattice (n=%s, T=%s, h=%s, J=%s)'%((n), (T), (h), (J)))
	plt.show()


def plot_hex_lattice(matrix, T, h, J):	#plots matrix as triangular lattice
	n=matrix.shape[0]
	up_x=[]
	up_y=[]
	down_x=[]
	down_y=[]
	
	i=0
	while i<=(n-1):
		j=0
		while j<=(n-1):
			i=float(i)
			j=float(j)
			if matrix[i,j]==-1:
				down_x.append((j)+(n-1-i)*0.5)
				down_y.append((0.5/0.559)*(n-1-i))#mult by wierd factor makes points equidistant
			if matrix[i,j]==1:
				up_x.append((j)+(n-1-i)*0.5)
				up_y.append((0.5/0.559)*(n-1-i))

			#print i,j, j+0.5

			j=j+1
		i=i+1
	
	plt.plot(up_x,up_y,'o',color='black', markersize=10)
	plt.plot(down_x,down_y, 'o', color='white', markersize=10)
	plt.title('Hexagonal Spin Lattice (n=%s, T=%s, h=%s, J=%s)'%((n), (T), (h), (J)))
	plt.xlim(-0.5,n-1+n*0.5)
	plt.ylim(-0.5,n-0.5)
	plt.show()


