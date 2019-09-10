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
    parser.add_argument("-xlabel", "--xlabel",  type=str,  default="crowder")
    
    args = parser.parse_args()
    return args


men_means, men_std = (20, 35, 30, 35, 27), (2, 3, 4, 1, 2)
women_means, women_std = (25, 32, 34, 20, 25), (3, 5, 2, 3, 3)




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
    olygomers = ('1', '3', '4', '6', '7', '8')
    crowders =('1', '2', '3', '4', '5', '6', '7', '8')
    args  = parseArguments()
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
    if args.prop != "hb":
        ax.set_xticklabels(crowders)
    else:
        ax.set_xticklabels(olygomers)

    ax.legend()

    autolabel(rects1, "left")
    autolabel(rects2, "right")
   
    plt.savefig('%s.pdf' % (args.prop) , bbox_inches="tight")

