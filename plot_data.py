#!/usr/bin/env python

import matplotlib.pyplot as plt
import math as m

#code that reads in a file_name from command line
#provided data inside file_name is in columns of numbers
#this program will plot the first column versus the other columns on different graphs
#will output M, E, X, C, S versus T plots if columns are T, M, E, X, C, S values

print "enter file to read data from, data must be stored in 8 columns of numbers"
name_of_file=raw_input()

fo=open(name_of_file, "rw+")

lines_list = fo.readlines()

T_vals=[]
M_vals=[]
M_err=[]
E_err=[]
E_vals=[]
X_vals=[]
C_vals=[]
S_vals=[]

i=0
while i<len(lines_list):
	entry=lines_list[i].split()
	T_vals.append(float(entry[0]))
	M_vals.append(m.fabs(float(entry[1])))
#	M_err.append(m.sqrt(float(entry[3])*float(entry[0])))	#for standard deviation in M
	E_vals.append(float(entry[2]))
#	E_err.append(m.sqrt(float(entry[4]))*float(entry[0]))	#for standard deviation in E
	X_vals.append(float(entry[3]))
	C_vals.append(float(entry[4]))
	S_vals.append(float(entry[5]))
	i=i+1

plt.plot(T_vals,M_vals, 'o')
plt.plot(T_vals,M_vals)
plt.title("Magnetisation versus Temperature")
plt.xlabel("Temperature [J/k]")
plt.ylabel("Magnetisation per site, |<M>| [mu]")
#plt.errorbar(T_vals, M_vals, yerr=M_err)#plots errorbars
plt.ylim(0.0, 1.0)
plt.show()

plt.plot(T_vals,E_vals, 'o')
plt.plot(T_vals,E_vals)
plt.title("Energy versus Temperature")
plt.xlabel("Temperature [J/k]")
plt.ylabel("Energy per site, <E> [J]")
#plt.errorbar(T_vals, E_vals, yerr=E_err)
plt.show()


plt.plot(T_vals, X_vals, 'o')
plt.plot(T_vals, X_vals)
plt.title("Magnetic Susceptibility versus Temperature")
plt.xlabel("Temperature [J/k]")
plt.ylabel("Magnetic Susceptibility, Xi [mu/k]")

plt.show()


plt.plot(T_vals,C_vals, 'o')
plt.plot(T_vals,C_vals)
plt.title("Specific Heat versus Temperature")
plt.xlabel("Temperature [J/k]")
plt.ylabel("Specific Heat, Cv [J/k**2]")

plt.show()

plt.plot(T_vals,S_vals, 'o')
plt.plot(T_vals,S_vals)
plt.title("Entropy versus Temperature")
plt.xlabel("Temperature [J/k]")
plt.ylabel("Entropy, S/S_max")

plt.show()


fo.close()

