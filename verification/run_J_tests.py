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
 test_generic( 'testA' )
 test_generic( 'testB' )
 test_generic( 'testC' )
 test_generic( 'testD' )
# test_generic( 'testE' )   MPI test. run with other script
 test_generic( "testF" )  
 test_generic( "testG" )
 test_generic( "testH" )
 test_generic( "testI" )
 test_generic( "testJ" )
 test_generic( "testK" )
 test_generic( "testL" )
 test_generic( "testM" )
 test_generic( "testN" )
 test_generic( "testO" )
 test_generic( "testP" )
 test_generic( "testQ" )
 test_generic( "testR" )
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
 for i in range( num_problems ):
    cleanup( cpoutfiles[i], problem_dir )
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
 cygwin  = False
 if 'CYGWIN' in os_name: cygwin = True
 if "WARP3D_HOME" in os.environ:
 	 pass
 else:
 	 print("\n\n>>> Enviornment variable not set: WARP3D_HOME")
 	 print("    Please set this to define install directory of WARP3D and re-try...")
 	 print("    Program exiting...\n\n")
 	 exit(0)
#
 if windows :
   continuation = '&'
   warp_name = '%WARP3D_HOME%"/run_windows/warp3d.exe" '
   print(">>> Platform: Windows\n" )
#
 if linux :
   continuation = ';'
   print(">>> Platform: Linux\n" )
   choice = input(">>> Intel Fortran version (=0), GFortran version(=1):") 
   if choice == "0" :
    warp_name = '$WARP3D_HOME"/run_linux/warp3d_Intel.omp" '
   if choice == "1" :
    warp_name = '$WARP3D_HOME"/run_linux/warp3d_gfortran.omp" '
#
 if osx :
   continuation = ';'
   print(">>> Platform: Mac OSX\n" )
   choice = input(">>> Intel Fortran version (=0), GFortran version(=1):")
   if choice == "0" :
      warp_name = '$WARP3D_HOME"/run_mac_os_x/warp3d.omp" '
   if choice == "1" :
      warp_name = '$WARP3D_HOME"/run_mac_os_x/warp3d_gfortran.omp" '
#
 if cygwin :
   continuation = ';'
   print(">>> Platform: Cygwin (Windows). Only Intel Fortran version supported\n" )
   warp_name = '$WARP3D_HOME"/run_windows/warp3d.exe" '
#
 str_threads = str(input(">>> Number of threads to use: " )  )
 threads = 'set OMP_NUM_THREADS='+str_threads + continuation + \
           'set MKL_NUM_THREADS='+str_threads
 run_warp = threads + continuation + warp_name
# 
 print("\n>>> Note: comparison and output values may be" \
               " different in the last 1 or 2")
 print("          signficant digits with various number of threads\n")
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
  m = ["All Problems ( ~ 500 secs wallclock on 2018 Intel Xeon, 10 threads )"]
  m.append( "Test A: SE(T), LEFM, Thermal, 20-node, FGMs, Face Loading, P-strain")
  m.append( "Test B: Globally rotated model")
  m.append( "Test C: SE(T), Small Eps, Thermal, 20-node, Rotated Model, FGMs, P-strain")
  m.append( "Test D: SE(T), GEONL, Thermal, 20-node, Rotated Model, P-strain, FGMs, Focused Mesh")
  m.append( "Test E: <available>")
  m.append( "Test F: SE(T), LEFM, Thermal, 8-node, P-strain, FGMs, blunt front")
  m.append( "Test G: SE(T), Small Eps, Thermal, 8-node, P-strain, FGMs, blunt front")
  m.append( "Test H: SE(B), Small Eps, Bending, 8-node, deformation/mises, P. Strains, FGMs, blunt front")
  m.append( "Test I: SE(B), Small Eps, Bending, 20-node, deformation/mises, P-strain, FGMs, blunt front")
  m.append( "Test J: SE(B), Small/NLGEOM, Bending, 8-node, P-strain, FGMs, blunt front, refined mesh")
  m.append( "Test K: SE(T), Small/NLGEOM, Thermal, 8-node, P-strain, FGMs, blunt front, refined mesh")
  m.append( "Test L: SE(B), Small/NLGEOM, Bending, 20-node, P-strain, FGMs, blunt front, refined mesh")
  m.append( "Test M: SE(T), Small/NLGEOM, Thermal, 20-node, P-strain, FGMs, blunt front, refined mesh")
  m.append( "Test N: SE(B), Small/NLGEOM, Bending, 8-node, 3D:3-layers, FGMs, blunt front, refined mesh")
  m.append( "Test O: SE(B), Residual Stresses, Initial-State, 8-node, P. Strain")
  m.append( "Test P: Pipe, Weld, Eigenstrain, 8-node, P. Strain")
  m.append( "Test Q: Simulated weld, release nodes, load as M(T), 8-node P. Strain")
  m.append( "Test R: Crystal plasticity, SSY- P. strain")
  m.append( "Quit" )
#
  print("> Select a problem to run:\n")
  if menu_displayed == False :
    for number, entry in enumerate( m ):
       print(" ",number+1,")",entry)
    menu_displayed = True
#
  schoice = input("Enter your choice (<return> for menu): " ) 
  if len(schoice) == 0  : 
    while True:
      print(" ")
      for number, entry in enumerate( m ):
         print(" ",number+1,")",entry)
      menu_displayed = True
      schoice = input("Enter your choice (<return> for menu): " ) 
      if len(schoice) > 0:  break
  bad_choice = False
  choice = int( schoice )
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
   if choice == 2: test_generic( 'testA' )
   if choice == 3: test_generic( 'testB' )
   if choice == 4: test_generic( 'testC' )
   if choice == 5: test_generic( 'testD' )
   if choice == 6: test_generic( 'testE' )
   if choice == 7: test_generic( 'testF' )
   if choice == 8: test_generic( 'testG' )
   if choice == 9: test_generic( 'testH' )
   if choice == 10: test_generic( 'testI' )
   if choice == 11: test_generic( 'testJ' )
   if choice == 12: test_generic( 'testK' )
   if choice == 13: test_generic( 'testL' )
   if choice == 14: test_generic( 'testM' )
   if choice == 15: test_generic( 'testN' )
   if choice == 16: test_generic( 'testO' )
   if choice == 17: test_generic( 'testP' )
   if choice == 18: test_generic( 'testQ' )
   if choice == 19: test_generic( 'testR' )
#
# ----------------------------------------------------------------------------
#
#    run main
#
# ----------------------------------------------------------------------------
#
if __name__ == "__main__": main()
exit(0)
