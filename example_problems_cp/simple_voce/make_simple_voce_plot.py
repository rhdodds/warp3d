from numpy import *
import matplotlib
import matplotlib.pyplot as plt
import os

# ----------------------------------------------------------------------------
#
#    function plot_set_up -- call before any plt.plot(...)
#
# ----------------------------------------------------------------------------
#
def plot_set_up( type ):
#
#      type:
#       0   1 plot landscape
#       1   1 plot portrait
#
   paper_width = 11.0   # inches
   paper_height = 8.5
   adjust_left = 0.2
   adjust_right = 0.8
   adjust_bottom = 0.2
   adjust_top = 0.8
   if type == 1 :
     print " set portrait"
     paper_width = 8.5  # inches
     paper_height = 11.0
     adjust_left = 0.2
     adjust_right = 0.9
     adjust_bottom = 0.2
     adjust_top = 0.8
     
   params = {'legend.fontsize': 12,
      'legend.linewidth': 0.75,
      'legend.frameon': False,
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
   plt.rcParams.update( params )
   plt.subplots_adjust( left=adjust_left, right=adjust_right,
        bottom=adjust_bottom, top=adjust_top ) # nice margins on page
   return
#
# ----------------------------------------------------------------------------
#
#    function colors
#
# ----------------------------------------------------------------------------
#
def colors( selection ):
#                   0    1    2    3    4   5 
   xplot_colors = ["r", "g", "b", "m", "c", "k" ]
   return xplot_colors[selection]

#
# ----------------------------------------------------------------------------
#
#    function lines
#
# ----------------------------------------------------------------------------
#
# 0   solid
# 1   dashed
# 2   dotted
# 3   dash-dotted
#
def lines( selection ):
   xplot_lines = ["-", "--", ":", "-.", "None" ]
   return xplot_lines[selection]


#
# ----------------------------------------------------------------------------
#
#    function symbols
#
# ----------------------------------------------------------------------------
 
# 0   "None"              
# 1   "."	point
# 2   ","	pixel
# 3   "o"	circle
# 4   "v"	triangle_down
# 5   "^"	triangle_up
# 6   "<"	triangle_left
# 7   ">"	triangle_right
# 8   "1"	tri_down
# 9   "2"	tri_up
# 10  "3"	tri_left
# 11  "4"	tri_right
# 12  "8"	octagon
# 13  "s"	square
# 14  "p"	pentagon
# 15  "*"	star
# 16  "h"	hexagon1
# 17  "H"	hexagon2
# 18  "+"	plus
# 19  "x"	x
# 20  "D"	diamond
# 21  "d"	thin_diamond
#
def symbols( selection ):
   xplot_symbols = ["None",  ".", ",", "o", "v", 
                    "^", "<", ">", "1",
                    "2", "3", "4", "8", "s", "p", 
                    "*", "h", "H", "+", "x", "D", "d" ]
   return xplot_symbols[selection]
 
 

#
# ----------------------------------------------------------------------------
#
#    function plot_start
#
# ----------------------------------------------------------------------------
#
def plot_start( type ):
#
#           1 plot on 11.0 x 8.5 (landscape, type = 0)
#           or 8.5 x 11 (portrait, type = 1)
#           has nice margins around, plot is
#           centered on page with rectangular shape
#           saved as pdf with transparent background
#
#           PDF has a transparent background
#
   plot_set_up( type )
   plt.figure()
   plot_set_up( type )
#   
   return

#
# ----------------------------------------------------------------------------
#
#    function plot_finish
#
# ----------------------------------------------------------------------------
#
def plot_finish( plot_file_name ):
#
   plt.annotate(' ',xy=(0.02, 0.92), xycoords='axes fraction')
   plt.grid(False)
#
   if os.path.isfile( plot_file_name ):
      os.remove( plot_file_name )
   plt.savefig( plot_file_name, transparent=True)
#
   return


#
# ----------------------------------------------------------------------------
#
#    main
#
# ----------------------------------------------------------------------------
#
   
def main():

# matplotlib.use("PDF")  # non-interactive plot making
 set_printoptions(threshold='nan')
#
 directory = './'

 nsteps = 145
 num_ele = 1

 max_num_rows = 100000  # max number model elements
 num_data_cols = 26  # stresses
 adata = zeros( [max_num_rows,num_data_cols], float64 )
 sig_xx = zeros( [max_num_rows], float64 )
 eps_pls = zeros( [max_num_rows], float64 )
#

 for i in xrange(nsteps):
  step = i + 1
  fname =  'wes' + str(step).zfill(5) + '_text'
  adata = loadtxt(fname, dtype=float64,
         comments='#').reshape(num_ele,num_data_cols)
  sig_xx[i] = adata[0,0]  
  eps_pls[i] = adata[0,10]  
  print  "\t\t Selected values last row: %15.7f  %15.7f"%(
     sig_xx[i], eps_pls[i] )
#
#
#     make a plot
#
#
 orientation = 0  #  landscape
 orientation = 1  #  portrait
 orientation = 0 
# 
 plot_start( orientation )
#
 plt.ylim( [0.0, 2000.] )
 plt.xlim( [0.0, 0.25] )
 plt.ylabel('True stress (MPa)')
 plt.xlabel('Log plastic strain')
 symbsize=8
# plt.title( "Effects of Cavity Options" )

# =======
 plt.plot( eps_pls[0:nsteps-1], sig_xx[0:nsteps-1],
      linestyle=lines(0), markersize=symbsize, marker=symbols(0),
      color=colors(1) )

#
#         legend, annotations, etc. finish plot
#
 plt.legend(loc='upper right')
#             
 plot_finish( directory+"z.pdf" )

 print "... normal end of program"


 exit(0)

# ----------------------------------------------------------------------------
#
#    run main
#
# ----------------------------------------------------------------------------
#
if __name__ == "__main__": main()

