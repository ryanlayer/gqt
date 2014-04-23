#!/usr/bin/python
import sys
import numpy as np
import matplotlib
import pylab
import random
from optparse import OptionParser

#lumpy_color_1='#003366'
#lumpy_color_2='#0080ff'
#lumpy_color_2='#0059b3'
#lumpy_color_3='#99ccff'
#lumpy_color_3='#007fff'
#lumpy_color_4='#9accff'
#gasv_color='#99ff99'
#delly_color='#ff9999'
#pindel_color='#ffff99'
#pindel_color='#bababa'

#plot_colors=[ \
#        delly_color, \
#        gasv_color, \
#        lumpy_color_1, \
#        lumpy_color_2, \
#        lumpy_color_3, \
#        lumpy_color_4, \
#        pindel_color]
#
#plot_markers=[ 'o', \
#                'o', \
#                'o', \
#                '^', \
#                's', \
#                'D', \
#                'o']
#

parser = OptionParser()

parser.add_option("-d",
                  "--data_file",
                  dest="data_file",
                  help="Data file")

parser.add_option("-o",
                  "--out_prefix",
                  dest="out_prefix",
                  help="Output file prefix")

parser.add_option("--y_max",
                  dest="y_max",
                  type="int",
                  help="Output file prefix")

parser.add_option("--y_min",
                  dest="y_min",
                  type="int",
                  help="Output file prefix")

parser.add_option("--x_max",
                  dest="x_max",
                  type="int",
                  help="Output file prefix")

parser.add_option("--x_min",
                  dest="x_min",
                  type="int",
                  help="Output file prefix")

parser.add_option("--minor",
                  dest="minor",
                  help="Output file prefix")

(options, args) = parser.parse_args()

if not options.data_file:
    parser.error('Data file not given')

if not options.out_prefix:
    parser.error('Out prefix not given')

if not options.y_max:
    parser.error('y_max not given')

if not options.y_min:
    parser.error('y_min not given')

if not options.x_max:
    parser.error('x_max not given')

if not options.x_min:
    parser.error('x_min not given')

if not options.minor:
    parser.error('minor not given')

minor = [int(x) for x in options.minor.split(":")]


f = open(options.data_file,'r')

L1={}
L2={}

L1_name = ''
L2_name = ''

for l in f:
    a = l.rstrip().split()

    L1_name = a[1]
    L1[int(a[0])] = [float(x) for x in a[2].split(',')]

    L2_name = a[3]
    L2[int(a[0])] = [float(x) for x in a[4].split(',')]

    

f.close()

print L1
print L2

#width = 0.11 
#matplotlib.rcParams.update({'font.size': 12})
#fig = matplotlib.pyplot.figure(figsize=(4.5,2.25),dpi=300)
#fig.subplots_adjust(wspace=.05,left=.01,bottom=.01)
#
#c=1
##print sorted(R.keys())
###for coverage in sorted(R.keys()):
##for coverage in [40.0]:
#
#
#max_fp =  max([max(R[r][1]) for r in R])
#min_fp =  min([min(R[r][1]) for r in R])
#max_tp =  max([max(R[r][0]) for r in R])
#min_tp =  min([min(R[r][0]) for r in R])
#
#max_tp = options.y_max
#min_tp = options.y_min
#
#max_fp = options.x_max
#min_fp = options.x_min
#
#ax = fig.add_subplot(1,1,c)
#ax.set_ylim([min_tp-100,max_tp+100])
#ax.set_xlim([min_fp-100,max_fp+100])
#
#print min_fp,max_fp,min_tp,max_tp
#
#minorLocator   = matplotlib.ticker.FixedLocator(range(minor[0],minor[1]+1,minor[2]))
#ax.yaxis.set_minor_locator(minorLocator)
#ax.yaxis.grid(b=True, which='major', color='0.75', linestyle='-')
#ax.yaxis.grid(b=True, which='minor', color='0.75', linestyle='--')
#ax.set_axisbelow(True) 
#
#line_width = 0.5
#marker_size = 5
#black_size = 0.75
#
#ps=[]
#i=0
#print R.keys()
#for r in sorted(R.keys()):
#    print R[r][1]
#
#    #ps.append(ax.plot(R[r][1],R[r][0],'.',
#        #markeredgecolor='none',
#        #color='black', 
#        #markersize=marker_size,
#        #marker='v',
#        #fillstyle='none'))
#
#    ps.append(ax.plot(R[r][1],R[r][0],'.',
#        markeredgecolor='none',
#        marker=plot_markers[i],
#        color='black',
#        markersize=marker_size+black_size,
#        fillstyle='full'))
#
#
#
#    ps.append(ax.plot(R[r][1],R[r][0],'-',
#        color='black', 
#        linewidth=line_width + black_size,
#        fillstyle='none'))
#
#    ps.append(ax.plot(R[r][1],R[r][0],'.',
#        markeredgecolor='none',
#        marker=plot_markers[i],
#        color=plot_colors[i],
#        markersize=marker_size,
#        fillstyle='full'))
#
#    ps.append(ax.plot(R[r][1],R[r][0],'-',
#        color=plot_colors[i],
#        linewidth=line_width,
#        markersize=6.5,
#        fillstyle='none'))
#
#
#
#
#
#
#    #ps.append(ax.plot(R[r][1],R[r][0],'.-',color=plot_colors[i],
#    #ps.append(ax.plot(R[r][1],R[r][0],'-',color=plot_colors[i],
#        #color='black', 
#        #linewidth=1,
#        #markersize=7,
#        #fillstyle='none'))
#
##    ps.append(ax.plot(R[r][1],R[r][0],'.-',color=plot_colors[i],
##        linewidth=1,
##        markersize=7,
##        markeredgecolor='black',
##        fillstyle='none'))
#
##    ps.append(ax.plot(R[r][1],R[r][0],'-',color=plot_colors[i],
##        linewidth=1,
##        markersize=7,
##        markeredgecolor='black',fillstyle='none'))
#
#    i+=1
##
##
##    l_pe_p = ax.bar(ind, l_pe, width, color=lumpy_color_1)
##    l_sr_p = ax.bar(ind+width*1, l_sr, width, color=lumpy_color_2)
##    l_pesr_p = ax.bar(ind+width*2, l_pesr, width, color=lumpy_color_3)
##    g_p = ax.bar(ind+width*3, g, width, color=gasv_color_1)
##    gp_p = ax.bar(ind+width*4, gp, width, color=gasv_color_2)
##    d_pe_p = ax.bar(ind+width*5, d_pe, width, color=delly_color_1)
##    d_sr_p = ax.bar(ind+width*6, d_sr, width, color=delly_color_2)
##    d_pesr_p = ax.bar(ind+width*7, d_pesr, width, color=delly_color_3)
##
##
##    #matplotlib.pyplot.title(str(int(coverage)) + 'x')
#matplotlib.pyplot.tick_params(axis='x',length=5,top="off",direction="out")
#matplotlib.pyplot.tick_params(axis='y',length=0)
##ax.get_xaxis().set_ticks_position('none')
#ax.get_xaxis().tick_bottom()
#ax.get_yaxis().set_ticks_position('none')
#ax.tick_params(axis='both', direction='in', length=3)
##ax.tick_params(axis='both', direction='in', width=0)
#
#
##ax.set_xticks(ind+width*3)
##ax.set_xticklabels([])
##ax.set_yticklabels([])
#ax.spines['top'].set_visible(False)
#ax.spines['right'].set_visible(False)
#ax.spines['left'].set_visible(False)
#
#matplotlib.pyplot.savefig(options.out_prefix,bbox_inches='tight')
