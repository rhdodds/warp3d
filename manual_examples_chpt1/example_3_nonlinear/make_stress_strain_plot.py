from numpy import *
from scipy.integrate import trapz, cumtrapz
import matplotlib
import matplotlib.pyplot as plt
import os


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
#
#    begin main
#
# ----------------------------------------------------------------------------
#
set_printoptions(threshold='nan')

directory = './'
data = zeros([2,7],float64)
pls_eps = zeros( 2, float64 )
stresses = zeros( [2,6], float64)

data = loadtxt( directory+'input_plot_stress_strain_curves.txt',
   dtype='float64', comments='#', usecols=(0,1,2,3,4,5,6) )
num_points = 2
#
pls_eps[0] = data[0,0]
pls_eps[1] = data[1,0]
for j in range(1,7):
  stresses[0,j-1] = data[0,j]
  stresses[1,j-1] = data[1,j]

#for j in range(2):
print pls_eps[0], pls_eps[1]
print stresses
print data
plot_start()
#
plt.ylim( [0.0, 800.0] )
plt.xlim( [0.0,3.5] )
plt.xlabel('Plastic Strain')
plt.ylabel('Uniaxial Stress (MPa)')
plt.title( " " )
      
plt.plot( pls_eps,stresses[:,0], color=plot_colors[0], label="Temperature Dependent" )
plt.plot( pls_eps,stresses[:,1], color=plot_colors[0], label="Temperature Dependent" )
plt.plot( pls_eps,stresses[:,2], color=plot_colors[0], label="Temperature Dependent" )
plt.plot( pls_eps,stresses[:,3], color=plot_colors[1], label="Temperature Dependent" )
plt.plot( pls_eps,stresses[:,4], color=plot_colors[1], label="Temperature Dependent" )
plt.plot( pls_eps,stresses[:,5], color=plot_colors[1], label="Temperature Dependent" )


#plt.plot( times_matlab, matlab, 's',markersize=3, color=plot_colors[4],label="MatLab" )

#plt.legend(loc='upper left')
#   
plot_file_name = 'plot_stress_strain.pdf'
plot_finish()



exit(0)        


