#!/usr/local/bin/python3


#
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
import os
from scipy.interpolate import make_interp_spline, BSpline
from scipy import interpolate

#
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
#
#    functions to make engineering plots with Dodds specific style
#    available for anyone to use
#
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
#
#
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
     print(" set portrait")
     paper_width = 8.5  # inches
     paper_height = 11.0
     adjust_left = 0.2
     adjust_right = 0.9
     adjust_bottom = 0.2
     adjust_top = 0.8

   params = {'legend.fontsize': 12,
#      'legend.linewidth': 0.75,   deprecated
      'legend.frameon': False,
      'legend.numpoints': 1,
      'figure.figsize': (paper_width,paper_height),
      'lines.linewidth': 0.75,   # width of data curves
      'axes.linewidth': 1.125,
      'axes.titlesize': 20,       # plot title
      'axes.labelsize': 15,
      'axes.labelcolor': 'k',
      'xtick.top': True,         # ticks top of graph
      'xtick.bottom': True,
      'xtick.direction': 'in',    # ticks point into graph
      'xtick.major.size': 10,     # major tick size in points
      'xtick.minor.size': 5,      # minor tick size in points
      'xtick.minor.visible': True,
      'xtick.major.pad': 6,       # distance to major tick label in points
      'xtick.minor.pad': 4,       # distance to the minor tick label in points
      'xtick.color': 'k',         # color of the tick labels
      'xtick.labelsize': 13,      # fontsize of the tick labels
      'ytick.direction': 'in',    # ticks point inside graph
      'ytick.major.size': 10,     # major tick size in points
      'ytick.minor.size': 5,      # minor tick size in points
      'ytick.minor.visible': True,
      'ytick.major.pad': 6,       # distance to major tick label in points
      'ytick.minor.pad': 4,       # distance to the minor tick label in points
      'ytick.color': 'k',         # color of the tick labels
      'ytick.labelsize': 13,    # fontsize of the tick labels
      'ytick.left': True,         # ticks top of graph
      'ytick.right': True }

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
#
#           0    1    2    3    4    5          6        7
#  colors: "r", "g", "b", "m", "c", "k", "darkviolet", "teal"
#             8              9            10          11
#         "darkgreen", "saddlebrown", "darkred", "darkmagenta",
#            12        13
#         "indigo",  "mediumblue"
#
#                   0    1    2    3    4   5
   xplot_colors = ["r", "g", "b", "m", "c", "k",
     "darkviolet",  # 6
     "teal",        # 7
     "darkgreen",   # 8
     "saddlebrown", # 9
     "darkred",     #10
     "darkmagenta", #11
     "indigo",      #12
     "mediumblue",  #13
     ]
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
# 1   "."	small filled circle
# 2   ","	pixel    just a point
# 3   "o"	circle   larger filled circle
# 4   "v"	filled triangle_down
# 5   "^"	filled triangle_up
# 6   "<"	filled triangle_left
# 7   ">"	filled triangle_right
# 8   "1"	tri_down
# 9   "2"	tri_up
# 10  "3"	tri_left
# 11  "4"	tri_right
# 12  "8"	filled octagon
# 13  "s"	filled square
# 14  "p"	filled pentagon
# 15  "*"	filled star
# 16  "h"	filled hexagon1 w/ border
# 17  "H"	filled hexagon1 w/ border
# 18  "+"	plus
# 19  "x"	x
# 20  "D"	filled diamond
# 21  "d"	filled thin_diamond
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
#
# ----------------------------------------------------------------------------
#
#    function plot_spline
#
# ----------------------------------------------------------------------------
#
def plot_spline( x, y, npts_x_y, npts_plot_on_curve= 100,
                 curve_label= " ", curve_color= 5,
                 curve_marker= 1, marker_size= 6 ,
                 line_style= 0 ):
#
#       vector of new points on X to evaluate spline fit. fit
#       spline using input number of points
#
 xnew = linspace( x.min(), x.max(), num=npts_plot_on_curve )
 tck =  interpolate.splrep(x, y, s=0)
 power_smooth = interpolate.splev(xnew, tck, der=0)
#
#       plot spline curve w/o markers, include label. label line
#       has right style, color but no marker point
#
 plt.plot( xnew, power_smooth, 
      linestyle=lines(line_style), markersize=0,
      color=colors(curve_color) )
#
#       plot input data points, no line, no label
#
 plt.plot( x[0:npts_x_y], y[0:npts_x_y],
      linestyle='None', markersize=marker_size, 
      marker=symbols(curve_marker),
      color=colors(curve_color),label=curve_label )
#
 return
# ----------------------------------------------------------------------------

#
# ----------------------------------------------------------------------------
#
#    function set_color_names
#
# ----------------------------------------------------------------------------
#
def plot_set_color_names():
#
 global red, green, blue, magenta, cyan, black, darkviolet
 global teal, darkgreen, saddlebrown, darkred, darkmagenta
 global mediumblue 
#
 red = 0
 green = 1
 blue = 2
 magenta = 3
 cyan = 4
 black = 5
 darkviolet = 6
 teal = 7
 darkgreen = 8
 saddlebrown = 9
 darkred = 10
 darkmagenta = 11
 mediumblue = 13
#
 return
#
# ----------------------------------------------------------------------------
#
#    function make_smooth
#
# ----------------------------------------------------------------------------
#
#        cubic bspline fit to data vectors (xvals,yvals).
#        return vectors of numpoints for plotting (x_smooth, y_smooth)
#
def make_smooth( xvals, yvals, numpoints ):
#
 x_smooth = linspace(amin(xvals),amax(xvals), num=numpoints, endpoint=True)
 spline_function = make_interp_spline(xvals,yvals, k=3)  # type: BSpline
 y_smooth = spline_function(x_smooth)
#
 return x_smooth, y_smooth


