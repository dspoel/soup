#!/usr/bin/env python3

import argparse, os, glob, math
import matplotlib.pyplot as plt
import numpy as np


width  = 7.84
height = 13/2.
crowder_name=["TufA", "MetE", "IcdA", "AhpC", "CspC", "Ppa ", "GapA", "Eno", "tRNA$^{Phe}$"]

fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8) = plt.subplots(8, sharey=False, sharex=True)



fig.set_figheight(height)
fig.set_figwidth(width)
fs=16
c=0
# for n in range(1,9):
#     c=1+c
#     for base_dir in ["../singles", "../cyto"]:
#         L = base_dir.split("../")[1]
#         if L == "singles":
#             lab="Dilute"
#         elif L == "cyto":
#             lab = "Cyto."
#         xvg = "%s/rmsf_ca*crowder_n%s*xvg" % (base_dir, str(n))
#         rmsf_xvg = glob.glob(xvg)
#         N = len(rmsf_xvg)
#         V=[]
#         for file in rmsf_xvg:
#             v = []
#             with open(file, 'r') as infile:
#                 for line in infile:
#                     if line.find("@") < 0 and line.find("#") < 0:
#                         data = line.split()
#                         v.append(float(data[1]))
#             if n == 4:
#                 v = v[0:163]
#             V.append(v)
#         fig_s = plt.subplot(4,2,c)

#         plt.errorbar(list(range(1,len(v)+1)), np.mean(V, axis=0), np.std(V, axis=0)/math.sqrt(N-1), label=lab)
#         plt.tick_params(labelsize=fs-2);
#         plt.ylabel(crowder_name[c-1],fontsize=fs)
#         if n == 1:
#             fig_s.legend(fontsize=fs-3)
#         if n == 2:
#             plt.yticks([0, .5, 1])
#         elif n == 7 or n == 8:
#             plt.xlabel('Residue', fontsize=fs)
# fig.text(-0.03,0.6,'RMSF ($nm$)',fontsize=fs+2, rotation=90)        
# fig.tight_layout()
# plt.savefig("rmsf.pdf", bbox_inches="tight")



for base_dir in ["../singles", "../cyto"]:
    L = base_dir.split("../")[1]
    with open('msf_%s.dat' % L , 'w') as outfile:
        outfile.write('%s sim.\tmean\tste\n' % (L))
        for n in range(1,9):            
            xvg = "%s/rmsf_ca*crowder_n%s*xvg" % (base_dir, str(n))
            rmsf_xvg = glob.glob(xvg)
            N = len(rmsf_xvg)
            # print(N)
            V=[]
            for file in rmsf_xvg:
                v = []
                with open(file, 'r') as infile:
                    for line in infile:
                        if line.find("@") < 0 and line.find("#") < 0:
                            data = line.split()
                            v.append(float(data[1])**2)
                if n == 4:
                    v = v[0:163]
                V.append(np.mean(v))
            outfile.write('crowder_n%d:\t%0.3f\t%0.3f\n' % (n, np.mean(V), np.std(V)/math.sqrt(len(V))) )




