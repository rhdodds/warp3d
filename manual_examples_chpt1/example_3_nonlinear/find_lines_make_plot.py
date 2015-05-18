from numpy import *
from scipy.integrate import trapz, cumtrapz
import matplotlib
import matplotlib.pyplot as plt
import os
import sys
sys.path.append("/Users/rdodds/bin")


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
   



def read_KJ_values( fname, column ):
#
 global KJ_values
 f = open(fname, "r") 
 count = 0

 while True:

  while True:
   line = f.readline()
   if not line:
     break
   if line.find( "domain   KI plane stress   ") > 0:
       line = f.readline()
       line = f.readline()
       words = line.split()
       count += 1
       KJ_values[count-1,column] = float( words[2] )
       KJ_values[count-1,column] = KJ_values[count-1,column]/sqrt(1000.0)
       print KJ_values[count-1,column] 
     
     
  if not line:
    break

 f.close()
 print "Normal EOF on:", filename  
 return count


# ----------------------------------------------------------------------------
#
#    begin main
#
# ----------------------------------------------------------------------------
#
set_printoptions(threshold='nan')



KJ_values = zeros([1000,10], float64 )
temps = zeros(1000,float64)


filename = "./output_temp_dependent"
num_values = read_KJ_values( filename, 0 )
print "number of values: ", num_values
filename = "./output_room_temp"
num_values = read_KJ_values( filename, 1 )
print "number of values: ", num_values
filename = "./output_steady_temp"
num_values = read_KJ_values( filename, 2 )
print "number of values: ", num_values
filename = "./output_CTE_clad_zero"
num_values = read_KJ_values( filename, 3 )
print "number of values: ", num_values

now_temp = 288.0
delta_temp = -4.0
temps = zeros(1000,float64)

for j in range(num_values):
      temps[j] = now_temp
      now_temp += delta_temp
      print temps[j]

plot_start()
#
plt.ylim( [0.0, 120.0] )
plt.xlim( [288.0, 90.0] )
plt.xlabel(r'$\mathrm{Inner\ Surface\ Temperature\ (C)}$')
plt.ylabel(r'$K_J\ \ \mathrm{(MPa\sqrt{m})}$')
plt.title( " " )
      
plt.plot( temps[0:num_values], KJ_values[0:num_values,0], color=plot_colors[0], label="Temperature Dependent" )
plt.plot( temps[0:num_values], KJ_values[0:num_values,1], color=plot_colors[1], label="20C" )
plt.plot( temps[0:num_values], KJ_values[0:num_values,2], color=plot_colors[2], label="250C" )
plt.plot( temps[0:num_values], KJ_values[0:num_values,3], color=plot_colors[3], label="Clad CTE=0" )
#plt.plot( times_matlab, matlab, 's',markersize=3, color=plot_colors[4],label="MatLab" )

plt.legend(loc='upper left')
#   
plot_file_name = 'plot.pdf'
plot_finish()



exit(0)        


