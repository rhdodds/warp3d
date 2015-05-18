from numpy import *
from scipy.integrate import trapz, cumtrapz
import matplotlib
matplotlib.use("PDF")  # non-interactive plot making
import matplotlib.pyplot as plt
import os
#
#
#        main program follows plotting function
#
#


#
# ----------------------------------------------------------------------------
#
#    plot_set_up -- call before any plt.plot(...)
#
# ----------------------------------------------------------------------------
#
def plot_set_up( type ):
   global plot_colors
#
#      type
#       0   1 plot landscape
#       1   1 plot portrait
#
   paper_width = 11.0   # inches
   paper_height = 8.5
   params = {'legend.fontsize': 12,
      'legend.linewidth': 0.75,
      'legend.frameon': True,
      'legend.numpoints': 1,
      'figure.figsize': (paper_width,paper_height),
      'axes.linewidth': 1.125,
      'axes.titlesize': 20,       # plot title
      'axes.labelsize': 16,
      'axes.labelcolor': 'k',
      'xtick.major.size': 10,     # major tick size in points
      'xtick.minor.size': 5,      # minor tick size in points
      'xtick.major.pad': 6,       # distance to major tick label in points
      'xtick.minor.pad': 4,       # distance to the minor tick label in points
      'xtick.color': 'k',         # color of the tick labels
      'xtick.labelsize': 13,      # fontsize of the tick labels
      'ytick.major.size': 10,     # major tick size in points
      'ytick.minor.size': 5,      # minor tick size in points
      'ytick.major.pad': 6,       # distance to major tick label in points
      'ytick.minor.pad': 4,       # distance to the minor tick label in points
      'ytick.color': 'k',         # color of the tick labels
      'ytick.labelsize': 13 }     # fontsize of the tick labels
   plt.rcParams.update(params)
   plt.subplots_adjust(left=0.2, right=0.8,
                    bottom=0.2, top=0.8) # nice margins on page
   plot_colors = ["r", "g", "b", "m", "c", "k" ]
   return
#
# ----------------------------------------------------------------------------
#
#    plot_start
#
# ----------------------------------------------------------------------------
#
def plot_start():
   global plot_colors

#           1 plot on 11.0 x 8.5 (landscape)
#           has nice margins around, plot is
#           centered on page with rectangular shape
#           saved as pdf with transparent background
#
#           PPT on OS X. drag-drop .pdf onto slide.
#           scaling then changes uniformly the size of
#           all fonts, line thicknesses etc.
#
   plot_set_up( 0 )
   plt.figure()
   plot_set_up( 0 )
#   
   return

#
# ----------------------------------------------------------------------------
#
#    plot_finish
#
# ----------------------------------------------------------------------------
#
def plot_finish():
   global plot_colors
#
   plt.annotate(' ',xy=(0.02, 0.92), xycoords='axes fraction')
   plt.grid(False)
#
   if os.path.isfile( plot_file_name ):
      os.remove( plot_file_name )
   plt.savefig( plot_file_name, transparent=True)
   return
#
   

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
#
#    begin main
#
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
#
set_printoptions(threshold='nan')
#
directory = './'

plot_start()
#
plt.ylim( [0.0, 3.2] )
plt.xlim( [0.0,3.0] )
plt.xlabel(r'Step Number (time)')
plt.ylabel(r'Punch Displacement')
plt.title( " " )
#
#         read data and plot data set 1
#
times, displ= loadtxt( directory+'punch_displ_time_plot_data.txt',
   dtype='float64', comments='#', usecols=(0,1), unpack=True )
num_points = count_nonzero( times ) +1
print ".. num points: ", num_points

plt.plot( times, displ,  'yo-', color=plot_colors[0])
plot_file_name = 'punch_displ_vs_time.pdf'
plot_finish()
exit(0)

