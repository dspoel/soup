#!/usr/bin/env python3

import argparse, os, glob, math
import matplotlib.pyplot as plt
import numpy as np

glrt = 1.618
width  = 5
height = width/glrt
fs=16
y = [15000, 20000, 25000, 30000, 35000, 40000]
ylbl =["15", "20", "25", "30", "35", "40"]

# for prop in ["volume", "density"]:
for prop in ["volume"]:
	fig, ax = plt.subplots()
	fig.set_figheight(height)
	fig.set_figwidth(width)
	for rep in range(1,4):
		v = []
		with open('shrinking_%s.%s.xvg' % (prop, str(rep)), 'r') as infile: 
			for line in infile:
				if line.find('#') < 0 and line.find("@") < 0:
					data = line.split()
					v.append(float(data[1]))
		plt.plot(v, label="replica "+str(rep))
	plt.xlabel('Time(ps)', fontsize=fs)
	if prop == "volume":
		plt.ylabel('%s ($10^3 nm^3$)' % prop, fontsize=fs)
		plt.yticks(y, ylbl)
		plt.text(-120,43100,'A', fontsize=fs)
	else:
		plt.ylabel('%s ($kg/m^3$)' % prop, fontsize=fs)
		plt.text(-80,1150,'B', fontsize=fs)

	plt.legend(fontsize=fs)
	plt.tick_params(labelsize=fs-2);
	plt.savefig('shrinking_%s.pdf' % prop, bbox_inches="tight")
	plt.close()
	# plt.plot(v, label=("replica %s" %(str(rep)) )