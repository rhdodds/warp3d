from numpy import *
from scipy.integrate import trapz, cumtrapz
import matplotlib
import matplotlib.pyplot as plt
import os
#
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
#
#    functions to make engineering plots with Dodds specific style
#    available for anyone to use
#
# ----------------------------------------------------------------------------

from plot_support import *

# ----------------------------------------------------------------------------
#
#    begin main
#
# ----------------------------------------------------------------------------
#
print(" ")
fname = "./input_plot_stress_strain_curves.txt"
file = open( fname, 'r' )
print(">> data file opened: ", fname)
#
#   skip comment lines
#
line = file.readline()
words = line.split()
num_comments = int( words[0] )
for i in range(num_comments):     # skip lines
   line = file.readline()
#
data = zeros([2,7],float64)
pls_eps = zeros( 2, float64 )
stresses = zeros( [2,6], float64)
nrows = 2
num_cols = 7
#
data = fromfile( file, dtype=float64, sep=" ", 
                  count=nrows*num_cols).reshape(nrows,num_cols)
#
num_points = 2
#
pls_eps[0] = data[0,0]
pls_eps[1] = data[1,0]
for j in range(1,7):
  stresses[0,j-1] = data[0,j]
  stresses[1,j-1] = data[1,j]
#
#print(( pls_eps[0], pls_eps[1]))
#print( stresses)
#print( data)
#
orientation = 0  #  landscape
orientation = 1  #  portrait
orientation = 0 
plot_start( orientation )
#
plt.ylim( [0.0, 800.0] )
plt.xlim( [0.0,3.5] )
plt.xlabel('Plastic Strain')
plt.ylabel('Uniaxial Stress (MPa)')
plt.title( " " )
#      
plt.plot( pls_eps,stresses[:,0], color=colors(0), label="Temperature Dependent" )
plt.plot( pls_eps,stresses[:,1], color=colors(0), label="Temperature Dependent" )
plt.plot( pls_eps,stresses[:,2], color=colors(0), label="Temperature Dependent" )
plt.plot( pls_eps,stresses[:,3], color=colors(1), label="Temperature Dependent" )
plt.plot( pls_eps,stresses[:,4], color=colors(1), label="Temperature Dependent" )
plt.plot( pls_eps,stresses[:,5], color=colors(1), label="Temperature Dependent" )
#   
plot_file_name = 'plot_stress_strain.pdf'
plot_finish(plot_file_name)
print("\n>>Plot file: ",plot_file_name)
#
print(">> normal termination")
#
exit(0)        


