#!/usr/local/bin/python3
#
#
#
from numpy import *
import os
#from scipy.interpolate import spline  Python 2
from scipy import interpolate
import subprocess
#
#
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
#
#    functions to make engineering plots with Dodds specific style
#    available for anyone to use
#
# ----------------------------------------------------------------------------
#
#from Dodds_plot_support import *
#
#
#              function set_step_list
#              ----------------------
#
def set_step_list( case ):
#
  first = 1
  increment = 1
  last = 793
  local_lst = arange( first, last+increment, increment ) 
  return len(local_lst), local_lst
#
#
#              function find_line
#              ------------------
#
def find_line( string  ):
#
   global f, endLocation

   debug = False
   if( debug ):
      print("..inside find_line")
      print("search string: ", string )
#
   while True:
        if( f.tell() == endLocation ):
          print("\n\n....... EOF reached .......")
          exit(0)
        line2 = f.readline( ) 
        if( debug ): print("     line2:",line2)
        if( line2.find( string ) == -1 ): continue
        return line2
#
#              function find/process a step
#              ----------------------------
#

def process_step( step ):
#
 global f, crack_growth, j_values
#
 debug = False # True
#
 while True:
#
#
    now_step = step
    da = crack_growth[step]
    ymod = 88000.0
    nu = 0.34
    J = j_values[step-1]
    K_I = sqrt(ymod*J*0.001/(1.0-nu*nu))
    print("{0:5.0f}  {1:8.2f}    {2:9.4f}   {3:8.4f}".format(  
        now_step, K_I, J, da ) )
#
    return
#
#    main
#
# ----------------------------------------------------------------------------
#
   
def main():
 global f, endLocation, crack_growth, j_values

 debug = False
 max_steps = 20000
#
 directory = "./"
#
#
 fname = directory + "crack_growth_values.dat"
 n = len(open(fname).readlines(  )) # number of data lines
 last_step = 0
 f = open(fname, "r") 
 for i in range(n):
   line = f.readline( ) 
   if ">> step, max da" in line:
      last_step = last_step + 1
 f.close()
 if debug: print("\n... number of steps with da: ", last_step)
#
 first = 1
 increment = 1
 step_list = arange( first, last_step+increment, increment ) 
#
 if( debug ):
  print(step_list)
#
#         read crack extension at each load step
#
 crack_growth = zeros( [max_steps], float64 )
 fname = directory + "crack_growth_values.dat"
 f = open(fname, "r") 
 line = f.readline( ) # 2 comment lines
 line = f.readline( ) 
 for i in range(last_step):
   line = f.readline( ) 
   words = line.split()
   delta_a = float64( words[5] )
   if debug: print("\n... step, da: ", i, delta_a)
   crack_growth[i+1] = delta_a
 f.close()

#
#         read J-values
#
 j_values = zeros( [max_steps], float64 )
 fname = directory + "j_values.inp"
 f = open(fname, "r") 
 line = f.readline()
 words = line.split()
 number_comments = int( words[0] )
 for i in range(number_comments):     # skip lines
    line = f.readline()
 line = f.readline()
 words = line.split()
 number_steps_J = int( words[0] )
 if debug: print( ".... steps w/ J-values: ", number_steps_J )
 for i in range(number_steps_J):
   line     = f.readline( ) 
   words    = line.split()
   j_values[i] = float64( words[0] )
 if debug: print( "   last of J-values: ",  j_values[number_steps_J-1] )
 f.close( )
#
#         header lines for output table
#
 print("5\n#")
 print("#     SSY, plane-strain. growth by element deletion")
 print("#     Coarse mesh")
 print("#\n#  step  K_I [Mpa-sqrt(m)]   J (N/mm)  Da (mm)")
 print(last_step+1)
 print("    0           0           0         0")
#
 for i in range(last_step):
  step = int( step_list[i] )
  if( debug ): print("... doing step: ", step)
  process_step( step )
#  f.seek(0)  # not needed since results are sequential by step
#   
 print("Normal end of program")
#
 exit(0)

# ----------------------------------------------------------------------------
#
#    run main
#
# ----------------------------------------------------------------------------
#
if __name__ == "__main__": main()

