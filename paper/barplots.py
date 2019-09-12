#!/usr/bin/env python3
import argparse
import matplotlib
import matplotlib.pyplot as plt
import numpy as np




def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-prop", "--prop",  type=str,  default="")
    parser.add_argument("-title", "--title",  type=str,  default="")
    parser.add_argument("-ylabel", "--ylabel",  type=str,  default="")
    parser.add_argument("-xlabel", "--xlabel",  type=str,  default="")
    
    args = parser.parse_args()
    return args


crowder_name=["TufA", "MetE", "IcdA", "AhpC", "CspC", "Ppa ", "GapA", "Eno", "RNA$^{Phe}$"]
def autolabel(rects, xpos='center'):
    """
    Attach a text label above each bar in *rects*, displaying its height.

    *xpos* indicates which side to place the text w.r.t. the center of
    the bar. It can be one of the following {'center', 'right', 'left'}.
    """

    ha = {'center': 'center', 'right': 'left', 'left': 'right'}
    offset = {'center': 0, 'right': 1, 'left': -1}

    # for rect in rects:
    #     height = rect.get_height()
    #     ax.annotate('{}'.format(height),
    #                 xy=(rect.get_x() + rect.get_width() / 2, height),
    #                 xytext=(offset[xpos]*3, 3),  # use 3 points offset
    #                 textcoords="offset points",  # in both directions
    #                 ha=ha[xpos], va='bottom')

if __name__ == '__main__':
    oligomers = (crowder_name[0],crowder_name[2],crowder_name[3],crowder_name[5],crowder_name[6],crowder_name[7])
    protein_name =crowder_name[0:8]

    args   = parseArguments()
    ylabel = args.ylabel
    if len(ylabel) == 0:
        if args.prop == "sasa":
            ylabel = "\$nm^{2}\$"
        elif args.prop == "rmsd" or args.prop == "rg":
            ylabel = "\$nm\$"
        elif args.prop == "msd":
            ylabel = '\$D (10^{-5} cm^2/s)\$'
        elif args.prop == "hb":
            ylabel = "#HB"
        else:
            print("Unknown property " + args.prop)
            exit(1)
        
    singles_mean, singles_std, cyto_mean, cyto_std = [], [], [], []
    for env in ["singles", "cyto"]:
        with open("%s_%s.dat" % (args.prop, env), 'r') as infile:
            for line in infile:
                if line.find("crowder_n") > -1 and line.find("-") <= 0:
                    line = line.split()
                    n = line[0].split("crowder_n")[1]
                    if env == "singles":
                        singles_mean.append(float(line[1]))
                        singles_std.append(float(line[2]))
                    if env == "cyto":
                        cyto_mean.append(float(line[1]))
                        cyto_std.append(float(line[2]))


    ind = np.arange(len(singles_mean))  # the x locations for the groups
    width = 0.35  # the width of the bars
    
    fig, ax = plt.subplots()
    rects1 = ax.bar(ind - width/2, singles_mean, width, yerr=singles_std,
                    label='Sing.')
    rects2 = ax.bar(ind + width/2, cyto_mean, width, yerr=cyto_std,
                    label='Cyto.')
    
    # # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel(args.ylabel)
    ax.set_xlabel(args.xlabel)
    ax.set_title(args.title)
    ax.set_xticks(ind)
    if args.prop == "hb":
        ax.set_xticklabels(oligomers)
    elif args.prop == "msd":
        ax.set_xticklabels(crowder_name)
    else:
        ax.set_xticklabels(protein_name)

    ax.legend()

    autolabel(rects1, "left")
    autolabel(rects2, "right")
   
    plt.savefig('%s.pdf' % (args.prop) , bbox_inches="tight")

