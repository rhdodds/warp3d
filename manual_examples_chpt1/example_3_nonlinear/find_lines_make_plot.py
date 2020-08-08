#
#
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
import os
#from scipy.interpolate import spline  Python 2
from scipy import interpolate
#
#
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
#
#    functions to make engineering plots with Dodds specific style
#    available for anyone to use
#
# ----------------------------------------------------------------------------

from plot_support import *



def read_KJ_values( fname, column ):
#
 global KJ_values
 f = open(fname, "r") 
 count = 0
 ldebug = 0

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
       if ldebug: print(KJ_values[count-1,column]) 
     
     
  if not line:
    break

 f.close()
 print("Normal EOF on:", filename)  
 return count


# ----------------------------------------------------------------------------
#
#    begin main
#
# ----------------------------------------------------------------------------
#
#set_printoptions(threshold='nan')

ldebug = 0

KJ_values = zeros([1000,10], float64 )
temps = zeros(1000,float64)


filename = "./output_temp_dependent"
num_values = read_KJ_values( filename, 0 )
print("number of values: ", num_values)
filename = "./output_20C_temp"
num_values = read_KJ_values( filename, 1 )
print("number of values: ", num_values)
filename = "./output_250C_temp"
num_values = read_KJ_values( filename, 2 )
print("number of values: ", num_values)
filename = "./output_CTE_clad_zero"
num_values = read_KJ_values( filename, 3 )
print("number of values: ", num_values)

now_temp = 288.0
delta_temp = -4.0
temps = zeros(1000,float64)

for j in range(num_values):
      temps[j] = now_temp
      now_temp += delta_temp
      if ldebug: print(temps[j])

orientation = 0  #  landscape
orientation = 1  #  portrait
orientation = 0 
plot_start( orientation )
#

#
plt.ylim( [0.0, 120.0] )
plt.xlim( [288.0, 90.0] )
plt.xlabel(r'$\mathrm{Inner\ Surface\ Temperature\ (C)}$')
plt.ylabel(r'$K_J\ \ \mathrm{(MPa\sqrt{m})}$')
plt.title( " " )
      
plt.plot( temps[0:num_values], KJ_values[0:num_values,0], color=colors(0), label="Temperature Dependent" )
plt.plot( temps[0:num_values], KJ_values[0:num_values,1], color=colors(1), label="20C" )
plt.plot( temps[0:num_values], KJ_values[0:num_values,2], color=colors(2), label="250C" )
plt.plot( temps[0:num_values], KJ_values[0:num_values,3], color=colors(3), label="Clad CTE=0",
          linestyle=lines(1))
#plt.plot( times_matlab, matlab, 's',markersize=3, color=plot_colors[4],label="MatLab" )

plt.legend(loc='upper left')
#   
plot_file_name = 'plot_KJ_vs_temp.pdf'
plot_finish( plot_file_name )
# 
print("\n\n... normal end of program")
#
exit(0)        


