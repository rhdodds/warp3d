#
filename = "./commands.inp"

print ("\n\n... writing compute-output commands to file: ", filename)
f = open( filename, 'w' )

start_step = 5
max_step = 52
delta_step = 1

step = start_step
while ( step <= max_step ):
  f.write( "\n" )
  f.write( "compute displacements for loading press_temp_trans step  " + str(step) + "\n" )
  f.write( "*input \"domain_define_compute.inp\"" + "\n" )
#  f.write( "output noheaders wide strains 21 stresses 21" + "\n" )
  step = step + delta_step
f.close()
print( "\n... file written ...\n\n")

exit(0)
