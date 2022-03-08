#!/usr/local/bin/python3
#
#from numpy import *
import datetime
import os
#
# ----------------------------------------------------------------------------
#
#    functions to make engineering plots with Dodds specific style
#    available for anyone to use
#
# ----------------------------------------------------------------------------
#
from Dodds_plot_support import *
#
# ----------------------------------------------------------------------------
#
#
# ----------------------------------------------------------------------------
#
#    function: get_converged_step
#
# ----------------------------------------------------------------------------
#
def get_converged_step( find_this_step ):
#
 global file
 for j in range(1000000):
   local_line = file.readline()
   if " >> solution for step: " in local_line:
      if local_line[30] == '.' : continue
      if " converged: " in local_line:
        words = local_line.split()
        now_step = int( words[4] )
        if now_step == find_this_step:
          return now_step
#
# ----------------------------------------------------------------------------
#
#    function find next J results, extract
#
# ----------------------------------------------------------------------------
#
def get_J_values( skip ):
#
 global file, count
#
#         find next instance of domain integral output.
#         skip prescribed number of domains, then line with value
#         for desired domain
#
# 

 for j in range(1000000):
   local_line = file.readline()
   if " domain        dm1 " in local_line:
      for k in range(skip):
       local_line = file.readline()
      local_line = file.readline()
      words = local_line.split()
      J_value = float( words[9] ) # column with total J
      return J_value
#
# ----------------------------------------------------------------------------
#
#    function next line in output with string
#
# ----------------------------------------------------------------------------
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
        line2 = f.readline( ) 
        if( debug ): print("     line2:",line2)
        if( line2.find( string ) == -1 ): continue
        return line2
#
# ----------------------------------------------------------------------------
#
#    function process_step -- extract output values for 1 step
#
# ----------------------------------------------------------------------------
#
def process_step( step ):
#
 global f, crack_growth, J_values, npts,  sig_flow, b_0_mm
#
 debug = False # True
#
 while True:
#
    force    = 0.0
    half_LLD = 0.0
#
    line = find_line( " >> total applied load pattern factors through step") 
#
    if( debug ): print("@0 line: ",line)
    words = line.split()
    if( debug ): print("@0 words:", words)
    fnow_step = float64( words[8] )
    if( fnow_step > float64( int(fnow_step) ) ): continue  # adaptive subincrement
    now_step = int( words[8] )
    if( debug ): print( ". now_step: ", now_step)
    if( now_step != step ): continue
    if( debug ): print(">>> found results step: ", now_step)
#
    line = find_line( "   > constraints  ")
    words = line.split()
    total_pin_displ = float( words[2] ) * 2.0
    if( debug ):print("\t (total) pin displacement: ", total_pin_displ)
#
    idir = 3 #  direction for reactions, LLD
    sign = -1.0
#
    line = find_line( "   Totals:  ")
    words = line.split()
    force_half = sign * float( words[idir] )  # reaction is Z
    force = 2.0 * force_half
    if( debug ):print("\t force: ", force)
#
#
    line = find_line( "  nodal displacements" )
    line = find_line( "  node               u   " )
    line = f.readline( )    
    words = line.split()
    half_LLD_center = sign * float64( words[idir] )
    line = f.readline( )    
    words = line.split()
    half_LLD_outside = -1.0 * float64( words[3] )
    half_LLD = 0.5 * ( half_LLD_center + half_LLD_outside )
    if( debug ):print("..half_LLD: ", half_LLD)
    total_LLD = 2.0 * half_LLD
#
    npts = npts + 1 # so that FE results for step 1 in position [1]
    J_value = J_values[npts]
    forces[npts] = force
    LLD[npts] = total_LLD
    pin_displ[npts] = total_pin_displ
    M_fe = b_0_mm * sig_flow / J_value
#
    print(
       "{0:5.0f}  {1:8.5f}    {2:8.5f}    {3:10.4f}   {4:10.2f}    {5:6.1f}".
        format( now_step, total_pin_displ, total_LLD, J_value, force, M_fe ) )
#
    return
#
# ----------------------------------------------------------------------------
#
#    function CT_limit_load -- closed form limit load estimat
#
# ----------------------------------------------------------------------------
#
def limit_load( b, a, W, B, sig_y ):
#
#
 bW = b/W
 t1 = 0.364 + 0.195*bW + 0.233*bW**2 - 0.0458*bW**3
 t2 = bW**2
 t3 = B*W*sig_y
 PL = t3 * t2 * t1
 return PL
 

#
# ----------------------------------------------------------------------------
#
#    main
#
# ----------------------------------------------------------------------------
#
  
def main():
#
 global f, file, count, J_values, npts
 global forces, LLD, pin_displ, b_0_mm, sig_flow
#
#         model output units: mm, N, MPa
#         during computations, they are converted to kips, inches
#          (yeah crazy USA units !)
#         printed tables and plots are back in kN, mm, MPa
#
 debug = False
 max_steps = 10000
 max_pts = max_steps
#
 directory = "./"  # results
 fname = directory + "woutput"
#
#             get  # completed steps with domain results
#
 file = open( fname, 'r' )
 last_step = 0
 for line in file:
    if "Completed domain integral" in line:
        last_step +=1
 file.close()
#
 file = open( fname, 'r' )
#
 MPa_to_ksi         = 0.145033
 mm_to_in           = 0.0393701
 kN_to_kips         = 0.224809 
 J_N_mm_to_J_kip_in = 0.00570125
#
 npts = 0
 step_incr = 1  #   results at these steps in output file
 now_step  = 0
#
 ymod_MPa = 88000.
 nu       = 0.4
#
 sig_y    = (450.0 + 500.0 ) / 2.0  # make sure consistent with inout file
 sig_flow = sig_y
#
 W_mm     = 50.0
 B_net_mm = W_mm / 2.0
 a_0_mm   = 0.5 * W_mm
 b_0_mm   = W_mm - a_0_mm
#
 ymod_ksi = MPa_to_ksi * ymod_MPa # -> ksi
 B_net    = mm_to_in * B_net_mm # -> in
 a_0      = mm_to_in * a_0_mm  # -> in
 W        = mm_to_in * W_mm # -> in
 b_0      = mm_to_in * b_0_mm
#
 print_limit_load = False
#
 mises = True   #  check with inout file: mises or deformation
#
 eta_pl_stable = 1.95   #    from these results
#
 PL = limit_load( W_mm-a_0_mm, a_0_mm, W_mm, B_net_mm, sig_y )
 if print_limit_load:
   print("\n  ..... Estimated limit load (N): {0:10.2f} " .format(PL) )
   print(" ")
#
 J_values  = zeros( [max_steps], float64 )   # for N, mm, MPa results
 forces    = zeros( [max_steps], float64 )
 LLD       = zeros( [max_steps], float64 )
 pin_displ = zeros( [max_steps], float64 )
#
 J_values_kip_in = zeros( [max_steps], float64 ) # for Kips, inches results
 forces_kip      = zeros( [max_steps], float64 )
 LLD_in          = zeros( [max_steps], float64 )
 pin_displ_in    = zeros( [max_steps], float64 )
#
#             process output file. get J-values
#
 for steps in range(50000):
   now_step = now_step + step_incr
   if now_step > last_step: break
   this_step = get_converged_step( now_step ) # positions read point in results file
   if this_step != now_step:
          print( ".. Error @ 1. this_step:",this_step)
          exit( 0 )
   J_value = get_J_values( 39 ) # domains to skip over
   if debug: print(".. steps, now_step, J",steps, now_step, J_value)
   npts = npts + 1 # 1st FE results in [1]
   J_values[npts] = J_value 
   if debug: print( ".. this_step, J: ", this_step, J_value)
#
 file.close()
 if debug: print("\n... J-values values extracted")
#
#             get load, LLD, ... values from output file
#
 print(17) # comment lines before data
 print("#")
 print("#")
 print("#   3D C(T), plane-sided, straight front, focused mesh")
 print("#")
 if mises: print("#     *** incremental plasticity ***")
 if not mises: print("#     *** deformation plasticity ***")
 print("#")
 print("#   nodes: 26754 elements: 23055")
 print("#   Unit: N, mm, N-mm") 
 print("#   1/4 symmetric specimen. no side-grooves")
 print("#   W = 50.0, a_0 = 25.0, b_0 = 25.0, a_0/W = 0.5, B = 25")
 print("#   14 layers across front (B/2)")
 print("#")
 print("#   Material: Zr-2.5Nb.Rm Temp. power-law hardening n=10")
 now = datetime.datetime.now()
 print("#  ",now.strftime("%Y-%m-%d %H:%M:%S"))
 print("#   ")
 print("#  step   Pin V       LLD          J_FE        force     M_fe")
 print(npts+1)
 print("    0       0           0            0            0      100000 ")
#
 f = open(fname, "r") 
 now_step  = 0
 npts = 0
#
 for i in range(50000):
   now_step = now_step + step_incr
   if now_step > last_step: break
   step = now_step
   if( debug ): print("... doing step: ", step)
   process_step( step ) #  <<<< call function to get step values
#
#                        compute eta-plastic value. calc in English
#                        units - yeah I know ...
#
 print("\n\n")
 print("..... Compute an eta-plastic value at each load step")
 print("         value should reach constant at high levels of")
 print("         plastic deformation")
 print(" ")
 print("         units: N, mm")
 print(" ")
 debug = False
#
 i = 0
 while i <= npts:
   J_values_kip_in[i] = J_N_mm_to_J_kip_in * J_values[i] 
   forces_kip[i]      = kN_to_kips * ( forces[i] * 0.001 )
   LLD_in[i]          = mm_to_in * LLD[i]
   pin_displ_in[i]    = mm_to_in * pin_displ[i]
   i = i + 1
#
#        K-factor and compliance from computed values over first few steps
#
 use_step = 4 
 J =  J_values_kip_in[use_step]  # [0] is the (0,0) data pair
 p = forces_kip[use_step]
 lld = LLD_in[use_step]
 elastic_compliance = lld / p
 KI = sqrt(ymod_ksi*J/(1.-nu**2))
 K_factor = KI/p
#
 print( "\t>> computed K_I factor (ksi sqrt(in)/kip): {0:10.5f}".format( K_factor ) )
 print( "\t>> computed elastic compliance (in/kip): {0:10.7f}".format(
      elastic_compliance ) )
#
 i = 1  # step 1 in FE analysis
 im1 = 0
 area = 0.0 
#
 print("\tstable eta_pl value set for comparing")
 print("\tJ_estimated vs. J_fe: {0:5.2f}".format( eta_pl_stable ) )
 print("\n")
#
 print("    i     i-1    force_i         J_fe        J_pl    eta_pl   J_stable_eta_pl   J_eta_stable/J_fe")
 while i <= npts: # num steps in FE analysis
#   
   force_i   = forces_kip[i]       # kips
   force_im1 = forces_kip[i-1]     # kips
   v_im1     = LLD_in[im1]         # inches
   v_i       = LLD_in[i]           # inches
#
   KI        = force_i * K_factor # ksi-sqrt(in)
   J_fe      = J_values_kip_in[i] 
   J_elastic = KI*KI*(1.0-nu*nu)/ymod_ksi # kip-in/in^2
   J_plastic = J_fe - J_elastic
   if J_plastic < 0: J_plastic = 0.0
#
   force_avg = 0.5 * ( force_i + force_im1 )
   dv = v_i - v_im1
   area = area + force_avg * dv
   v_elastic = force_i * elastic_compliance
#
   area_plastic = area - 0.5 * force_i * v_elastic
   if area_plastic < 0.0 : area_plastic = 0.0
   eta_pl = 1.0
   if area_plastic > 0.00001:
       eta_pl = J_plastic * B_net * b_0 / area_plastic
#
   J_from_eta_pl = J_elastic + eta_pl_stable * area_plastic / B_net / b_0
#   
   s1 =  "{0:5.0f}  {1:5.0f}    {2:8.2f}    {3:10.4f}   {4:10.4f}"
   s2 =  "  {5:5.2f}  {6:10.4f}            {7:6.3f} "
   s3 = s1+ s2
   print( s3.format(
      i, im1, forces[i], J_values[i], J_plastic/J_N_mm_to_J_kip_in, eta_pl, 
      J_from_eta_pl/J_N_mm_to_J_kip_in,
      J_from_eta_pl/J_N_mm_to_J_kip_in /J_values[i]) )  
   i = i + 1
   im1 = im1 + 1  # end of while
#
#-------------------------------------------------------------------------------------
#
#                        make plots
#
#-------------------------------------------------------------------------------------
#
#
#             Force vs. LLD  :   kN, mm
#
#
 orientation = 0  #  landscape
 orientation = 1  #  portrait
 orientation = 0
#
 plot_start( orientation )
#
 plt.ylim( [0.0, 100.0] )
 plt.xlim( [0.0, 6.0] )
 plt.xlabel(r'$\Delta_{\rm{LLD}}\ (\rm{mm})$')
 plt.ylabel('Force (kN)')
#
 symbsize=8
#
#  npts = number of FE load steps, result vectors have (0,0)
#
 xvals = zeros( [npts+1], float64 )
 yvals = zeros( [npts+1], float64 )
#
 i = 0
 while i <= npts:
     xvals[i] = LLD[i]
     yvals[i] = forces[i] * 0.001
     i = i + 1
#
 plt.plot( xvals, yvals,
       linestyle=lines(0), markersize=symbsize, marker=symbols(0),
       color=colors(2),label=r'$a_0=10.0$ mm, incr. plasticity' )
#
 plt.legend(loc='lower right')
 plot_finish( "plot_results_force_displ.pdf" )
#
#
#             J  vs.   LLD     N/mm  vs. mmm
#
#
 orientation = 0  #  landscape
 orientation = 1  #  portrait
 orientation = 0
#
 plot_start( orientation )
#
 plt.ylim( [0.0, 1200.] )
 plt.xlim( [0.0, 6.0] )
 plt.xlabel(r'$\Delta_{\rm{LLD}}\ (\rm{mm})}$')
 plt.ylabel(r'$J\ \rm{(N/mm)}$')
#
 symbsize=8
# plt.title( "current a = 10.05 mm" )
#
#
#  npts = number of FE load steps, result vectors have (0,0)
#
 xvals = zeros( [npts+1], float64 )
 yvals = zeros( [npts+1], float64 )
#
 i = 0
 while i <= npts:
     xvals[i] = LLD[i]
     yvals[i] = J_values[i]
     i = i + 1
#
 plt.plot( xvals, yvals,
       linestyle=lines(0), markersize=symbsize, marker=symbols(0),
       color=colors(5),label=r'$a_0=10.0$ mm, incr. plasticity' )
#
 plt.legend(loc='lower right')
#
 plot_finish( "plot_results_J_LLD.pdf" )
#
 print("\nNormal end of program\n\n")
#
 exit(0)


# ----------------------------------------------------------------------------
#
#    run main
#
# ----------------------------------------------------------------------------
#
if __name__ == "__main__": main()

