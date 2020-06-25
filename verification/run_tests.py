import subprocess
import os
import platform
import sys
import time
#
#
#
#
#              see driver.inp in each problem directory for data
#              to rive running test problems and checking values.
#
#
#              function check_results
#              ----------------------
#
def check_results( problem_dir, results_file, nsearchlines, 
                   search_lines, skip_lines, answer, anspos ):         
#
 global  count_diffs, failed_search, message_diff
#
#
#              -> open file
#              -> loop over search lines
#                  - find line in the output file
#              -> handle skip lines. which line has data ?
#                    > 0 skip this number of lines after last search line,
#                        then read line containing answer
#                    = 0 the last searched line is also the line with data
#                   = -1 read next line after last search line. it has data
#              -> extract data value from line and compare to
#                 known good value
#                  
#
 ldebug = False
 print("\n  ... Checking results: " + results_file)
#
 with open( problem_dir + '/' + results_file, 'rt' ) as outfile:
   for sline in range( nsearchlines ):
     if ldebug: 
       print("\t\t ldebug2:", sline, search_lines[sline])
     line, cnt, found = find_line( search_lines[sline], outfile  )
     if found == False : print( failed_search ); return
   if skip_lines == -1: line = outfile.readline()
   if skip_lines > 0 :
      for i in range( skip_lines ):
          line = outfile.readline()
      line = outfile.readline()
   words = line.split()
   print("\t    .... comparison value:\t", answer )
   message = " "
   if answer != words[anspos] : message = message_diff; count_diffs += 1
   print("\t    .... value from output:\t", words[anspos], message, "\n" )
#
 return
#
#              function run_all
#              ----------------
#
def run_all():
#
 time_start = time.time()
#
 test_generic( 'test14' )
 test_generic( 'test18' )
 test_generic( 'test24' )
 test_generic( 'test39' )
 test_generic( 'test41' )
 test44()
 test_generic( "test47" )  
 test_generic( "test86" )
 test_generic( "test50" )
 test_generic( "test51" )
 test_generic( "test54" )
 test_generic( "test57" )
 test_generic( "test60" )
 test_generic( "test61" )
 test_generic( "test63" )
 test_generic( "test67" )
 test_generic( "test69" )
 test_generic( "test70" )
 test_generic( "test71" )
 test_generic( "test72" )
 test_generic( "test73" )
 test_generic( "test74" )
 test_generic( "test75" )
 test_generic( "test76" )
 test_generic( "test77" )
 test_generic( "test78" )
 test_generic( "test80" )
 test81()
 test_generic( "test82/voche_model" )
 test_generic( "test82/mts_model" )
 test_generic( "test82/mts_model_multi" )
 test_generic( "test82/ornl_model" )
 test_generic( "test82/mrr_model" )
 test_generic( "test82/mrr_model_diff1B" )
 test_generic( "test82/mrr_model_diff2A" )
 test_generic( "test82/mrr_model_diff2B" )
 test_generic( "test82/mrr_model_diff3B" )
 test_generic( "test82/mrr_model_diff4B" )
 test_generic( "test82/djgm_model" )
 test_generic( "test82/djgm_overlap_taylor" )
 test_generic( "test82/djgm_taylor" )
 test_generic( "test82/djgm_hard_work")
 test_generic( "test83" )
 test_generic( "test84" )
 test_generic( "test85" )
 test_generic( "test86" )
 test_generic( "test87" )
#
 time_end = time.time()
 print("  ... Elaspsed walltime: %0.1f" % (time_end-time_start))
 print("\n >>> Number of tests w/ detected solution differences: ",
               count_diffs)
 print(" ")

 return
#
#              function get_new_line  (skip comments and blanks lines)
#              ---------------------
#

def get_new_line( infile ):
 
 lwords = []
 while True:
    line = infile.readline()
    if len(line.strip()) == 0 : continue # skip blank lines
    lwords = line.split()
    if lwords[0] == "*" : continue
    if lwords[0] == "!" : continue
    if lwords[0] == 'c' : continue
    if lwords[0] == "C" : continue
    if lwords[0] == "#" : continue
    if lwords[0] == "/" : continue
    break
 return line, lwords
#
#              function test_generic
#              ---------------------
#
def test_generic( problem_dir ):
#
 global debug, run_warp, count_diffs, failed_search, message_diff
#
#
#
#              problem_dir must have a driver.inp file with 
#              data to execute test problems and check results
#
 ldebug = False
 infile = open( problem_dir + "/" + "driver.inp", 'rt' )
 line, words = get_new_line( infile )
 title1 = line.rstrip()
 title2 = infile.readline().rstrip()
 print("\n>>> ",title1); print("    ",title2)
#
#
#              get number of problems to solve in this directory
#              get problem ids, output files, input files
#
 line, words = get_new_line( infile )
 num_problems = int(words[0])
#
 test_ids = [ ]
 outfiles = [ ]
 infiles  = [ ] 
 line = " "
 for i in range( num_problems ):
   line, words = get_new_line( infile )
   test_ids.append( line.rstrip() )
 for i in range( num_problems ):
   line, words = get_new_line( infile )
   outfiles.append( line.rstrip()  )
 for i in range( num_problems ):
   line, words = get_new_line( infile )
   infiles.append( line.rstrip()  )
#
 cpoutfiles = outfiles.copy()  # for deleting working files
#
#            run the problems
#
 print(" ")
 print("  ... Running jobs: ...")
 for i in range( num_problems ):
   run_str = run_warp + '<' + infiles[i] + '>' + outfiles[i]
   time_start = time.time()
   run_chk = subprocess.run( run_str, shell = True, cwd = problem_dir )
   time_end = time.time()
   print("  ... Ran: ", test_ids[i]," ... walltime: %0.1f" % (time_end-time_start))
 print("  ... WARP3D jobs completed...")
#
#            check results for each problem
#               -> number of search lines
#               -> string to find for each search line
#               -> number of lines to skip after last search line
#               -> string with known good answer
#               -> position on check line to extract solution value
#
 for i in range( num_problems ):
  if ldebug: print("<<< checking results problem: ", i )
  line, words = get_new_line( infile )
  nsearchlines = int( words[0] )
  search_lines = [ ]
  apostrophe = "'"
#
  for j in range( nsearchlines ): # store each search lone in list
     line, words = get_new_line( infile )
     line = line.rstrip();  line = line.replace(apostrophe, "",2)
     search_lines.append( line )
#
  line, words = get_new_line( infile )
  skip_lines = int( words[0] )
  line, words = get_new_line( infile )
  answer = words[0]
  line, words = get_new_line( infile )
  anspos = int( words[0] )
  if ldebug: 
     print("\t\t ldebug. # search lines: ", nsearchlines )
     for k in range( nsearchlines ):
        print("\t\t ldebug. search line: ", search_lines[k] )
     print("\t\t ldebug. skip_lines, anspos: ", skip_lines, anspos )
     print("\t\t ldebug. answer: ", answer )
#
  check_results( problem_dir, outfiles[i], nsearchlines, search_lines, 
                skip_lines, answer, anspos )
# |< end for over problems to check results
 infile.close()
#
#            delete all working files in problem directory
#
# for i in range( num_problems ):
#    cleanup( cpoutfiles[i], problem_dir )
#
 return

#
#              function test44
#              ---------------
#
def test44():
#
 global debug, run_warp, count_diffs, failed_search, message_diff
#
 print("\n>>> Test 44 (Cohesive-interface growth, restart)")
 print("    ============================================")
# 
 problem_dir = 'test44'
#
 num_problems = 4
 test_ids = ['test_44_ppr','test_44_ppr_restart','test_44_exp1','test_44_exp1_restart']
 outfiles = ['test_44a_out','test_44b_out','test_44c_out','test_44d_out']
 cpoutfiles = outfiles.copy()
 infiles = ['test_44_ppr','test_44_ppr_restart','test_44_exp1','test_44_exp1_restart']
#
 print(" ")
 print("  ... Running jobs: ...")
 for i in range( num_problems ):
   run_str = run_warp + '<' + infiles[i] + '>' + outfiles[i]
   time_start = time.time()
   run_chk = subprocess.run( run_str, shell = True, cwd = problem_dir )
   time_end = time.time()
   print("  ... Ran: ", test_ids[i]," ... walltime: %0.1f" % (time_end-time_start))
 print("  ... WARP3D jobs completed...")
#
 i = 1
 nsearchlines = 2
 search_lines = [ ' step:     200 converged' ]
 search_lines.append( 'Totals:' )
 skip_lines = 0
 answer = '0.33949465E+00'
 anspos = 2
 check_results( problem_dir, outfiles[i], nsearchlines, search_lines, 
                skip_lines, answer, anspos )
#
 i = 3
 nsearchlines = 2
 search_lines = [ ' step:     200 converged:  2 iters' ]
 search_lines.append( 'Totals: ' )
 skip_lines = 0
 answer = '0.35246131E+00'
 anspos = 2
 check_results( problem_dir, outfiles[i], nsearchlines, search_lines, 
                skip_lines, answer, anspos )
#
 for i in range( num_problems ):
    cleanup( cpoutfiles[i], problem_dir )
#
 return

#
#              function test81
#              ---------------
#
#         special function which examines contents of flat files
#
def test81():
#
 global debug, run_warp, count_diffs, failed_search, message_diff
#
 print("\n>>> Test 81 (Flat and Patran result output files)")
 print("    =============================================")
# 
 problem_dir = 'test81'
#
#         checks on results in flat & Patran result files.
#
 test_ids = ['test_81']
 outfiles = ['test_81_out']
 cpoutfiles = outfiles.copy()
 infiles = ['test_81']
#
 print(" ")
 print("  ... Running jobs: ...")
 run_str = run_warp + '<' + infiles[0] + '>' + outfiles[0]
 time_start = time.time()
 run_chk = subprocess.run( run_str, shell = True, cwd = problem_dir )
 time_end = time.time()
 print("  ... Ran: ", test_ids[0]," ... walltime: %0.1f" % (time_end-time_start))
 print("  ... WARP3D jobs completed...")
#
 i = 0
 answer = "-0.430241E-01"
 anspos = 0
 rfile = "wns0000050_text"
 print("\n  ... Checking results: " + rfile )
 with open( problem_dir + '/' + rfile ) as f:
    for line in f:
        pass
    words = line.split()
 print("\t    .... comparison value:\t" + " " + answer )
 message = " "
 if answer != words[anspos] : message = message_diff; count_diffs += 1
 print("\t    .... value from output:\t", words[anspos], message, "\n" )
#
 i = 1
 answer = "-0.202360E-03"
 anspos = 0
 rfile = "wne0000050_text"
 print("\n  ... Checking results: " + rfile )
 with open( problem_dir + '/' + rfile ) as f:
    for line in f:
        pass
    words = line.split()
 print("\t    .... comparison value:\t" + " " + answer )
 message = " "
 if answer != words[anspos] : message = message_diff; count_diffs += 1
 print("\t    .... value from output:\t", words[anspos], message, "\n" )
#
 i = 2
 answer = "-0.698211E-02"
 anspos = 0
 rfile = "wee0000050_text"
 print("\n  ... Checking results: " + rfile )
 with open( problem_dir + '/' + rfile ) as f:
    for line in f:
        pass
    words = line.split()
 print("\t    .... comparison value:\t" + " " + answer )
 message = " "
 if answer != words[anspos] : message = message_diff; count_diffs += 1
 print("\t    .... value from output:\t", words[anspos], message, "\n" )
#
 i = 3
 answer = "0.565348E+02"
 anspos = 0
 rfile = "wes0000050_text"
 print("\n  ... Checking results: " + rfile )
 with open( problem_dir + '/' + rfile ) as f:
    for line in f:
        pass
    words = line.split()
 print("\t    .... comparison value:\t" + " " + answer )
 message = " "
 if answer != words[anspos] : message = message_diff; count_diffs += 1
 print("\t    .... value from output:\t", words[anspos], message, "\n" )
#
 i = 4
 answer = "-0.130775E-02"
 anspos = 0
 rfile = "wnd0000050_text"
 print("\n  ... Checking results: " + rfile )
 with open( problem_dir + '/' + rfile ) as f:
    for line in f:
        pass
    words = line.split()
 print("\t    .... comparison value:\t" + " " + answer )
 message = " "
 if answer != words[anspos] : message = message_diff; count_diffs += 1
 print("\t    .... value from output:\t", words[anspos], message, "\n" )
#
 i = 5
 answer = "0.358056E-01"
 anspos = 0
 rfile = "wem0000050_text_creep"
 print("\n  ... Checking results: " + rfile )
 with open( problem_dir + '/' + rfile ) as f:
    for line in f:
        pass
    words = line.split()
 print("\t    .... comparison value:\t" + " " + answer )
 message = " "
 if answer != words[anspos] : message = message_diff; count_diffs += 1
 print("\t    .... value from output:\t", words[anspos], message, "\n" )

#
 i = 6
 answer = "0.228269E-04"
 anspos = 2
 rfile = "wnfd0000050"
 print("\n  ... Checking results: " + rfile )
 with open( problem_dir + '/' + rfile ) as f:
    for line in f:
        pass
    words = line.split()
 print("\t    .... comparison value:\t" + " " + answer )
 message = " "
 if answer != words[anspos] : message = message_diff; count_diffs += 1
 print("\t    .... value from output:\t", words[anspos], message, "\n" )

#
 i = 7
 answer = "-0.160559E-04"
 anspos = 0
 rfile = "wnfs0000050"
 print("\n  ... Checking results: " + rfile )
 with open( problem_dir + '/' + rfile ) as f:
    for line in f:
        pass
    words = line.split()
 print("\t    .... comparison value:\t" + " " + answer )
 message = " "
 if answer != words[anspos] : message = message_diff; count_diffs += 1
 print("\t    .... value from output:\t", words[anspos], message, "\n" )
#
 i = 7
 answer = "-0.100000E+01-0.930621E-05"
 anspos = 0
 rfile = "wnfe0000050"
 print("\n  ... Checking results: " + rfile )
 with open( problem_dir + '/' + rfile ) as f:
    for line in f:
        pass
    words = line.split()
 print("\t    .... comparison value:\t" + " " + answer )
 message = " "
 if answer != words[anspos] : message = message_diff; count_diffs += 1
 print("\t    .... value from output:\t", words[anspos], message, "\n" )
#
 cleanup( cpoutfiles[0], problem_dir )
#
 return
#
#              function initialize
#              -------------------
#
def initialize():
#
 global debug, run_warp, count_diffs, failed_search, message_diff, \
        windows, linux, osx, menu_displayed
#
 debug = 0
#
#          must have Python >= 3.5
#
 ver = sys.version_info
 print("\n\n>>> Running with Python version: " \
      + str(ver.major)+'.'+str(ver.minor)) 
 if ver.major == 3 and ver.minor <= 4:
     print('\n\n  >>>>>>>>>>>>> Please upgrade to Python 3.5 or newer')
     print('                Execution terminated\n\n')
     exit(0)
#
 os_name = platform.system()
 windows = os_name == 'Windows'
 linux   = os_name == 'Linux'
 osx     = os_name == 'Darwin'
#
 if windows :
   continuation = '&'
   warp_name = '%WARP3D_HOME%"/run_windows/warp3d.exe" '
   print(">>> Platform: Windows\n" )
#
 if linux :
   continuation = ';'
   warp_name = '$WARP3D_HOME"/run_linux/warp3d_Intel.omp" '
   print(">>> Platform: Linux\n" )
#
 if osx :
   continuation = ';'
   warp_name = '$WARP3D_HOME"/run_mac_os_x/warp3d.omp" '
   print(">>> Platform: Mac OSX\n" )
#
 str_threads = str(input(">>> Number of threads to use: " )  )
 threads = 'set OMP_NUM_THREADS='+str_threads + continuation + \
           'set MKL_NUM_THREADS='+str_threads
 run_warp = threads + continuation + warp_name
# 
 print("\n>>> Note: comparison and output values may be" \
               " different in the last 1 or 2")
 print("          signficant digits with various number of threads")
#
 count_diffs = 0
 failed_search ="\n\n\t\t ***** cannot find required output to check\n\n"
 message_diff =   "\t **** difference in solution"
 menu_displayed = False
#
 return


#
#              function display_menu
#              ---------------------
#
def display_menu():
#
  global debug, run_warp, count_diffs, failed_search, message_diff, \
         windows, linux, osx, menu_displayed
  local_debug = False
#
  m = ["All Problems"]
  m.append( "Test 14: (Linear elastic impact, sparse iterative solver)" )
  m.append( "Test 18: (Gurson model growth, impact loading, restart)" )
  m.append( "Test 24: (hollow sphere, skew constraints, int press, large displ) " )
  m.append( "Test 39: (CTOA and SMCS crack growth)" )
  m.append( "Test 41: (rigid contact, large displacements)" )
  m.append( "Test 44: (Cohesive-interface growth, restart)" )
  m.append( "Test 47: (Bilinear, cyclic, thermo-plasticity)" )
  m.append( "Test 86: (SEB, FGMs, nonlinear, J-values)" )
  m.append( "Test 50: (Fatigue crack growth using node release)" )
  m.append( "Test 51: (Nonlinear patch tests: tet elements)" )
  m.append( "Test 54: (FGM, linear, SE(T), non-global constraints)" )
  m.append( "Test 57: (SSY, mises + hydrogen effects)" )
  m.append( "Test 60: (Advanced cyclic plasticty model)" )
  m.append( "Test 61: (FGMs, interaction integrals)" )
  m.append( "Test 63: (360-degree SSY model, KI,KII,KIII)" )
  m.append( "Test 67: (Mesh tieing of two cylinders)" )
  m.append( "Test 69: (Interface-cohesive elements)" )
  m.append( "Test 70: (generalized_plasticity w/ temperature dependent properties)" )
  m.append( "Test 71: (piston loading for unsteady aerodynamic pressures)" )
  m.append( "Test 72: (user-defined integer lists of nodes-elements)" )
  m.append( "Test 73: (UMAT included in WARP3D - bilinear mises, kinematic hardening)" )
  m.append( "Test 74: (Crystal plasticity and restart)" )
  m.append( "Test 75: (Release constraints: absolute and relative)" )
  m.append( "Test 76: (User-defined multi-point constraints)" )
  m.append( "Test 77: (Finite strain transformations/output)" )
  m.append( "Test 78: (High-rate loading panel with hole. 20-node, rate plasticity)" )
  m.append( "Test 80: (Norton creep. Non-global constraints)" )
  m.append( "Test 81: (Flat and Patran result output files" )
  m.append( "Test 82: (Crystal Plasticity options)" )
  m.append( "Test 83: (bar2 element)" )
  m.append( "Test 84: (link2 element; Periodic Boundary Conditions)" )
  m.append( "Test 85: (user-defined initial-stresses)" )
  m.append( "Test 86: (FGMs compute J. explicit terms needed for path independence)" )
  m.append( "Test 87: (T-stress. surface cracked plate. face loading)" )
  m.append( "Quit" )
#
  print("> Select a problem to run:\n")
  if menu_displayed == False :
    for number, entry in enumerate( m ):
       print(" ",number+1,")",entry)
    menu_displayed = True
#
  choice = int( input("Enter your choice (0 for menu): " ) )
  if choice == 0 : 
    for number, entry in enumerate( m ):
       print(" ",number+1,")",entry)
    menu_displayed = True
    choice = int( input("Enter your choice (0 for menu): " ) )
  bad_choice = False
  if choice > len( m ) : 
      bad_choice = True
  last = False
  if choice == len( m ) : last = True
#
  return choice, bad_choice, last

#
#              function find_line
#              ------------------
#
def find_line( string, outfile, debug=False ):
#
   found = False
   nowline = -1
   cnt = -1
#
   for cnt, nowline in enumerate(outfile):
      if debug : print("Line {}: {}".format(cnt, nowline))
      if string in nowline :
        found = True
        break
#
   return nowline, cnt, found
#
#
#              function cleanup
#              ----------------
#
def cleanup( output_file, problem_dir ):
#
 global windows
#
 flist1 = 'wn* we* wm* *text *neutral *neut  states_header* *_db'
 flist2 =' *db *packets *tst states  energy step* model* ' 
 flist  = flist1 + flist2 + output_file
 if windows :
   delstr = 'del ' + flist1 + ' 2>nul'
   run_chk = subprocess.run( delstr, shell = True, cwd = problem_dir )
   delstr = 'del ' + flist2 + ' 2>nul'
   run_chk = subprocess.run( delstr, shell = True, cwd = problem_dir )
   delstr = 'del ' + output_file + ' 2>nul'
   run_chk = subprocess.run( delstr, shell = True, cwd = problem_dir )
 if windows == False : 
   delstr = '/bin/rm ' + flist + ' >/dev/null 2>&1'
   run_chk = subprocess.run( delstr, shell = True, cwd = problem_dir )
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
# 
 global debug, run_warp, count_diffs, failed_search, message_diff, \
        windows, linux, osx, menu_displayed
#	
 initialize()
#
 while True :  # over problem selection
#
   while True :  # get a single problem number to run/check
     choice, bad_choice, last = display_menu()
     if bad_choice :
        print("\n.... Bad menu entry. Please try again ...\n" )
     else:
        break
#
   if last : 
     print("\n >>> Number of tests w/ detected solution differences: ",
               count_diffs)
     print("\n >>> All done with tests ...\n\n" )
     exit(0)
#
   if choice == 1: run_all()
   if choice == 2: test_generic( 'test14' )
   if choice == 3: test_generic( 'test18' )
   if choice == 4: test_generic( 'test24' )
   if choice == 5: test_generic( 'test39' )
   if choice == 6: test_generic( "test41" ) 
   if choice == 7: test44()
   if choice == 8: test_generic( "test47" )
   if choice == 9: test_generic( "test86" )
   if choice == 10: test_generic( "test50" )
   if choice == 11: test_generic( "test51" )
   if choice == 12: test_generic( "test54" )
   if choice == 13: test_generic( "test57" )
   if choice == 14: test_generic( "test60" )
   if choice == 15: test_generic( "test61" )
   if choice == 16: test_generic( "test63" )
   if choice == 17: test_generic( "test67" )
   if choice == 18: test_generic( "test69" )
   if choice == 19: test_generic( "test70" )
   if choice == 20: test_generic( "test71" )
   if choice == 21: test_generic( "test72" )
   if choice == 22: test_generic( "test73" )
   if choice == 23: test_generic( "test74" )
   if choice == 24: test_generic( "test75" )
   if choice == 25: test_generic( "test76" )
   if choice == 26: test_generic( "test77" )
   if choice == 27: test_generic( "test78" )
   if choice == 28: test_generic( "test80" )
   if choice == 29: test81()
   if choice == 30: 
          test_generic( "test82/voche_model" )
          test_generic( "test82/mts_model" )
          test_generic( "test82/mts_model_multi" )
          test_generic( "test82/ornl_model" )
          test_generic( "test82/mrr_model" )
          test_generic( "test82/mrr_model_diff1B" )
          test_generic( "test82/mrr_model_diff2A" )
          test_generic( "test82/mrr_model_diff2B" )
          test_generic( "test82/mrr_model_diff3B" )
          test_generic( "test82/mrr_model_diff4B" )
          test_generic( "test82/djgm_model" )
          test_generic( "test82/djgm_overlap_taylor" )
          test_generic( "test82/djgm_taylor" )
          test_generic( "test82/djgm_hard_work")
   if choice == 31: test_generic( "test83" )
   if choice == 32: test_generic( "test84" )
   if choice == 33: test_generic( "test85" )
   if choice == 34: test_generic( "test86" )
   if choice == 35: test_generic( "test87" )
#
# ----------------------------------------------------------------------------
#
#    run main
#
# ----------------------------------------------------------------------------
#
if __name__ == "__main__": main()
exit(0)
