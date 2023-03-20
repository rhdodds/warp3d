#
import numpy
from scipy import interpolate
#
#    functions to make engineering plots with Dodds specific style
#    available for anyone to use
#
# ----------------------------------------------------------------------------
#
from Dodds_plot_support import *   #from Dodds_plot_support import *
#
# ----------------------------------------------------------------------------
#
#    function set color values - names
#
# ----------------------------------------------------------------------------
def set_color_symbol_line_values():
#
 global red, green, blue, \
        magenta, cyan, black, darkviolet, teal, \
        darkgreen, saddlebrown, \
        darkred, darkmagenta, mediumblue, \
        solid, dashed, dotted, dash_dotted, filled_square,    \
        sm_filled_circle, lg_filled_circle, filled_down_tri,  \
        filled_up_tri, filled_diamond

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
 solid       = 0
 dashed      = 1
 dotted      = 2
 dash_dotted = 3
#
 filled_square    = 13
 sm_filled_circle = 1
 lg_filled_circle = 3
 filled_down_tri  = 4
 filled_up_tri    = 5
 filled_diamond   = 20

#
 return
# ----------------------------------------------------------------------------
#
#    Cubic Hermite functions defined over 0 <= t <= 1.0. the functions
#    and 1st derivatives
#
# ----------------------------------------------------------------------------
#
def h00( t ):
 return (1.+2.*t)*(1.-t)**2
def h10( t ):
 return t*(1.-t)**2
def h01( t ):
 return t*t*(3.-2.*t)
def h11( t ):
 return t*t*(t-1.0)
def dh00( t ):
 return 6.0*t**2 - 6.0*t
def dh10( t ):
  return 3.0*t*t -4.0*t + 1.0
def dh01( t ):
  return -6.0*t*t + 6.0*t
def dh11( t ):
  return 3.0*t*t -2.0*t
#
# ----------------------------------------------------------------------------
#
#    User input, read curve points file
#
# ----------------------------------------------------------------------------
def get_input():
#
 global pts_table, num_curve_pts, ymod, fname, \
        sig_low_limit, sig_up_limit, eps_low_limit, eps_up_limit, \
        hprime_low, hprime_limit
#
 debug = False
 directory = "./"
#
 tfname = input("> file w/ curve points (use quit to stop, CR for same file): ")
 if tfname == "quit" : exit(0)
 if tfname == "" : tfname = fname
 string = input("> min,max stress; min,max plastic strain for plot (4 values, no commas): ")
 words = string.split()
#
 sig_low_limit = float64( words[0] )
 sig_up_limit  = float64( words[1] )
 eps_low_limit = float64( words[2] )
 eps_up_limit  = float64( words[3] )
 string = input("> min,max range for hprime in hprime vs. plastic strain plot (2 values): ")
 words = string.split()
#
 hprime_low = float64( words[0] )
 hprime_limit = float64( words[1] )
#
 fname = tfname
 file = open( fname, 'r' )
 line = file.readline()
 words = line.split()
 ymod = float64( words[0] )
 line = file.readline()
 words = line.split()
 num_curve_pts = int( words[0] )
 last_index = num_curve_pts- 1
 n = num_curve_pts
 pts_table = fromfile( file, dtype=float64, sep=" ", 
                  count=n*2).reshape(n,2)
 print( ".... last curve values from file: ", pts_table[n-1,0], pts_table[n-1,1] )
 file.close( )
 return
# ----------------------------------------------------------------------------
#
#    Convert total strain values to plastic strain values
#
# ----------------------------------------------------------------------------
#
def convert_to_plastic_strains():
#
 global pts_table
#
 i = 0
 while i < num_curve_pts :
   s =   pts_table[i,1]   # skip (0,0) point
   eps = pts_table[i,0]
   pts_table[i,0] = eps - s/ymod 
   pts_table[i,1] = s
   i += 1
# ----------------------------------------------------------------------------
#
#    Setup for WARP3D pchip interpolation
#
# ----------------------------------------------------------------------------
def setup_warp3d_pchip():
#
 global x, y, delta, m
#
 debug = False
#
 x = zeros( [num_curve_pts +1], float64 )  #    1 indexed
 y = zeros( [num_curve_pts +1], float64 )
 delta = zeros( [num_curve_pts +1], float64 )
 m = zeros( [num_curve_pts +1], float64 )
 k = 1
 while k <= num_curve_pts :
   x[k] = pts_table[k-1,0]  # 1 indexed
   y[k] = pts_table[k-1,1]   # 1 indexed
   k += 1
#
 if x[1] != 0.0 : x[1] = 0.0  # handle roundoff
#
 k = 1
 while k <= num_curve_pts  - 1:
   h =  x[k+1] - x[k]
   delta[k] = ( y[k+1] - y[k] ) / h
   k += 1
#
 k = 1
 while k <= num_curve_pts :
   if debug: print(".. k, x, y: ",k, x[k], y[k] )
   k += 1
#
 return
# ----------------------------------------------------------------------------
#
#     build SciPy pchip interpolator 
#
# ----------------------------------------------------------------------------
def build_SciPy_chip_interpolator():
#
 global pchip
#
 pchip = interpolate.PchipInterpolator(x[1:],y[1:] )
 return
# ----------------------------------------------------------------------------
#
#     build WARP3D pchip slopes at pts on user curve 
#
# ----------------------------------------------------------------------------
def build_warp3d_pchip_slopes():
#
 print(".... Building pchip slopes to use at pts on user curve")
#
 debug = False
#
 k = 2
 while k <= num_curve_pts  - 1:
   G = 0.0
   S1 = delta[k-1]
   S2 = delta[k]  
   if sign(S1) * sign(S2) > 0 :
      h1 = x[k] - x[k-1]
      h2 = x[k+1] - x[k]
      alpha1 = h1 + 2.*h2
      alpha2 = 3.0*(h1 + h2)
      alpha = alpha1 / alpha2
      denom = alpha * S2 + (1.-alpha)*S1
      G = S1*S2 / denom
   m[k] = G
   k += 1
#
 del1 = delta[1] #  one sided harmonic mean at first point. matches SciPy
 del2 = delta[2]
 h1 = x[2] - x[1]
 h2 = x[3] - x[2]
 d = ((2.0*h1+h2)*del1 - h1*del2)/(h1+h2)
 if sign(d) != sign(del1) :
   d = 0.0
 else:
   if (sign(del1) != sign(del2) ) and (abs(d)>abs(3.0*del1)):
     d = 3.0 * del1
 m[1] = d
 m[num_curve_pts] = delta[num_curve_pts -1]
#
 if debug:
   k = 1
   while k <= num_curve_pts :
    print(".. k, x, y, delta, m: ",k, x[k], y[k], delta[k], m[k] )
    k += 1
#
 return
# 
# ----------------------------------------------------------------------------
#
#     build plot pts using WARP3D pchip interpolator
#
# ----------------------------------------------------------------------------
def build_warp3d_fit_points():
#
 global num_strain_values, strain_vec, stress_vec_warp3d, \
        hprime_vec_warp3d
#
 debug = False
 print(".... Building WARP3D pchip interpolated points to plot")
#
 num_strain_values = 50000  # change as desired for plotting
 deps = (eps_up_limit - eps_low_limit)/num_strain_values
 now_eps = eps_low_limit
#
 strain_vec = zeros( [num_strain_values], float64 ) # zero indexed
 stress_vec_warp3d= zeros( [num_strain_values], float64 ) #    "
 hprime_vec_warp3d = zeros( [num_strain_values], float64 ) #    "    
#
 i = 0
 while i < num_strain_values:
   if debug: print("@4.1, i: ",i)
#
   k = 1   # find the segment
   while k <= num_curve_pts - 1 :
     if debug: print("@4.2. k, now_eps,x[k],x[k+1]: ",k, now_eps,x[k],x[k+1]  )
     if now_eps >= x[num_curve_pts]:
       strain_vec[i] = now_eps #linear extend last user segment
       stress_vec_warp3d[i] = y[num_curve_pts] 
       break       
     if now_eps >= x[k] and now_eps <= x[k+1] : # pchip interpolate
       if debug: print("@4.2")
       h2 = x[k+1] - x[k]
       t = ( now_eps - x[k] ) / h2 
       pchip_value = y[k]*h00(t) + h2*m[k]*h10(t) + y[k+1]*h01(t) + \
                     h2 * m[k+1]*h11(t)
       if debug: print("@4.3     now_eps, t, s: ", now_eps, t, s)
       strain_vec[i] = now_eps
       stress_vec_warp3d[i] = pchip_value
       hp = y[k]*dh00(t) + h2*m[k]*dh10(t) \
              + y[k+1]*dh01(t) + h2*m[k+1]*dh11(t)
       hprime_vec_warp3d[i] = hp / h2
       break
     k += 1
   i += 1
   now_eps += deps
#
 return
# ----------------------------------------------------------------------------
#
#     build plot pts using SciPy pchip interpolator
#
# ----------------------------------------------------------------------------
def build_SciPy_fit_points():
#
 global ypchip, hprime_scipy_vec
#
 print(".... Building SciPy pchip interpolated points to plot")
 ypchip = pchip(strain_vec)
 hprime_scipy_vec = zeros( [num_strain_values], float64 )
#
 k = 1     # use central difference form. skip 1st. last curve points
 while k <= num_strain_values - 2:
    deps = strain_vec[k+1] - strain_vec[k-1]
    dsig = ypchip[k+1] - ypchip[k-1]
    hprime_scipy_vec[k] = dsig/deps
    k += 1
#
 return

# ----------------------------------------------------------------------------
#
#     make the plot
#
# ----------------------------------------------------------------------------
def make_plot():
#
 debug = False
#
 print(".... Starting plot stress-strain")
 orientation = 0  #  landscape
 orientation = 1  #  portrait
 orientation = 0
#
 plot_start( orientation )
#
 plt.ylim( [sig_low_limit,sig_up_limit] )
 plt.xlim( [eps_low_limit,eps_up_limit] )
 plt.xlabel('Plastic Strain')
 plt.ylabel('Stress')
#
 symbsize = 8
#
 np = num_strain_values - 1
 if debug: print("... np: ", np)
 xvals = copy( strain_vec[0:np] )
 yvals = copy( stress_vec_warp3d[0:np] )
 plt.plot( xvals, yvals,
       linestyle=lines(solid), markersize=symbsize, marker=symbols(0),
       color=colors(black),label="WARP3D pchip fit" )
#
 plt.plot( strain_vec, ypchip,
       linestyle=lines(dashed), markersize=symbsize, marker=symbols(0),
       color=colors(green),label="SciPy pchip fit" )
#
 plt.plot( x[1:], y[1:],
       linestyle='none', markersize=symbsize, marker=symbols(1),
       color=colors(red),label="User Points" )
#
 plt.legend(loc='lower right')
#
 print(".... PDF for plot: stress_strain_curve.pdf")
 plot_finish( "stress_strain_curve.pdf" )
#
#
#
 print(".... Starting plot hprime")
 orientation = 0  #  landscape
 orientation = 1  #  portrait
 orientation = 0
#
 plot_start( orientation )
#
 plt.ylim( [hprime_low, hprime_limit] )
 plt.xlim( [eps_low_limit,eps_up_limit] )
 plt.xlabel('Plastic Strain')
 plt.ylabel('h-prime')
#
 symbsize = 5
#
 np = num_strain_values - 1
 xvals  = copy( strain_vec[0:np] )
 y1vals = copy( hprime_vec_warp3d[0:np] )
 plt.plot( xvals, y1vals,
       linestyle=lines(solid), markersize=symbsize, marker=symbols(0),
       color=colors(black),label="WARP3D pchip fit" )
#
 xvals = copy( strain_vec[1:np] ) # notice skip 1st point
 y2vals = copy( hprime_scipy_vec[1:np] )
 plt.plot( xvals, y2vals,
       linestyle=lines(solid), markersize=symbsize, marker=symbols(0),
       color=colors(green),label="SciPy fit" )

#
 xvals = copy( x[1:num_curve_pts] )
 yvals = copy( m[1:num_curve_pts] )
 plt.plot( xvals, yvals,
       linestyle='none', markersize=symbsize, marker=symbols(1),
       color=colors(red),label="PCHIP @ user pts" )
 plt.legend(loc='upper right')
#
 print(".... PDF for plot: hprime.pdf")
 plot_finish( "hprime.pdf" )
#
 return

#
# ----------------------------------------------------------------------------
#
#    main
#
# ----------------------------------------------------------------------------
#
#  
def main():
#
#
 set_color_symbol_line_values()
#
 while (True):
  get_input()
  convert_to_plastic_strains()
  setup_warp3d_pchip()
  build_SciPy_chip_interpolator()
  build_warp3d_pchip_slopes()
  build_warp3d_fit_points()
  build_SciPy_fit_points()
  make_plot()
#
 exit(0)

# ----------------------------------------------------------------------------
#
#    run main
#
# ----------------------------------------------------------------------------
#
if __name__ == "__main__": main()

