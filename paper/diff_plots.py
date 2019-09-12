#!/usr/bin/env python3

import argparse, os, glob, math
import matplotlib.pyplot as plt
import numpy as np

markers = ['o', 'v','*', 'd', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p' , 'h', 'H', '+', 'x', 'D',  '|', '_', 'P', 'X', ',', '.']

atom_crowder = [11890, 11825, 12906, 25592, 1009, 16284, 19976, 12814, 2377]
atom_s_b = sorted(atom_crowder)

order_crowder=[4, 7, 6, 3, 8, 1, 2, 9, 5] #biggest to smallest
order_s_b=[5, 9, 2, 1, 8, 3, 6, 7, 4] #smallest to biggest

avg_cyto, ste_cyto =[],[]
for i in range(0, len(order_s_b)):
	with open('msd_cyto.dat', 'r') as infile:
		for line in infile:
			if "crowder_n"+str(order_s_b[i]) in line:
				avg_cyto.append(float(line.split()[1]))
				ste_cyto.append(float(line.split()[2]))
plt.errorbar(atom_s_b, avg_cyto, ste_cyto, marker='o', markersize=10, linestyle=' ')
plt.xlabel('# atoms', fontsize=12)
plt.ylabel('$D_{cyto} (10^{-5} cm^2/s)$', fontsize=12)
plt.text(-5000,0.018,'B',fontsize=14)
plt.savefig('diff_cyto.pdf', bbox_inches="tight")
plt.close()


avg_sing, ste_sing = [],[]
for i in range(0, len(order_s_b)):
	with open('msd_singles.dat', 'r') as infile:
		for line in infile:
			if "crowder_n"+str(order_s_b[i]) in line:
				avg_sing.append(float(line.split()[1]))
				ste_sing.append(float(line.split()[2]))
plt.errorbar(atom_s_b, avg_sing, ste_sing, marker='o', markersize=10, linestyle=' ')
plt.xlabel('# atoms', fontsize=12)
plt.ylabel('$D_{0} (10^{-5} cm^2/s)$', fontsize=12)
plt.text(-5000,0.14, 'A',fontsize=14)
plt.savefig('diff_singles.pdf', bbox_inches="tight")
plt.close()

AVG, STD = [],[]

for i in range(0, len(avg_sing)):
	AVG.append(avg_cyto[i]/avg_sing[i])
	STD.append((avg_cyto[i]/avg_sing[i])* math.sqrt( (ste_cyto[i]/avg_cyto[i])**2 + (ste_sing[i]/avg_sing[i])**2 ) )


plt.errorbar(atom_s_b, AVG, STD, marker='o', markersize=10, linestyle=' ')
plt.xlabel('# atoms', fontsize=12)
plt.ylabel('$D_{cyto}/D_{0}$', fontsize=12)
# plt.text(-5000,0.2,'C',fontsize=14)
plt.savefig('diff_cyto_over_singles.pdf', bbox_inches="tight")
plt.close()

                    # std_vec.append((avg/avg_norm)*math.sqrt( (std/avg)**2 + (std_norm/avg_norm)**2 ) )



