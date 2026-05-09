#!/bin/bash  
/bin/rm -f  -f *_db   >& /dev/null
/bin/rm -f  -f *~   >& /dev/null
/bin/rm -f  w*   >& /dev/null
/bin/rm -f  core*  >& /dev/null
/bin/rm -f  energy*  >& /dev/null
/bin/rm -f  *messages  >& /dev/null
/bin/rm -f  *.warp   >& /dev/null
/bin/rm -f  packet*  >& /dev/null
/bin/rm -f  step_*  >& /dev/null
/bin/rm -f  solver*  >& /dev/null
/bin/rm -f  pardiso*  >& /dev/null
/bin/rm -f *.exo >& /dev/null
/bin/rm -f *.db >& /dev/null
/bin/rm -f *packets* >& /dev/null
/bin/rm -f *.text 
#  
#  
echo "> done..."
