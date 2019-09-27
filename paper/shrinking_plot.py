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
		plt.text(-120,43100,'B', fontsize=fs)
	else:
		plt.ylabel('%s ($kg/m^3$)' % prop, fontsize=fs)
		# plt.text(-80,1150,'B', fontsize=fs)

	plt.legend(fontsize=fs)
	plt.tick_params(labelsize=fs-2);
	plt.savefig('shrinking_%s.pdf' % prop, bbox_inches="tight")
	plt.close()
	# plt.plot(v, label=("replica %s" %(str(rep)) )


import matplotlib.image as mpimg
fig = plt.figure()
# print(fig.get_size_inches())
pic = ['box_big.png', 'box_small_.png']
x1=[0, 0]
x2=[0, 1000]
x3=[1000, 1000]
x4=[1000, 0]
image = mpimg.imread(pic[0])
ax1 = fig.add_subplot(121)
plt.imshow(image)
plt.axis('off')
plt.text(-100,-5, 'A', fontsize = fs)
plt.plot(x1,x2, 'k-')
plt.plot(x2,x3, 'k-')
plt.plot(x4,x1, 'k-')
plt.axis('off')

image = mpimg.imread(pic[1])
ax2 = fig.add_subplot(122)
plt.imshow(image)
# plt.title('After')
plt.axis('off')
# plt.plot(x1,x2, 'k-')
plt.plot(x2,x3, 'k-')
plt.plot(x3,x4, 'k-')
plt.plot(x4,x1, 'k-')
plt.axis('off')
fig.subplots_adjust(wspace=-.1, hspace=0)
# plt.tight_layout()
plt.savefig('box.pdf',bbox_inches="tight")
plt.close()


# plt.show()