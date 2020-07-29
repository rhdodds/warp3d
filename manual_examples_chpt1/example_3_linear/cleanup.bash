#!/bin/bash
#
#           cleanup directory after running analysis
#
/bin/rm -f  model.text wem* wes* wnd* wns* states_head*   >& /dev/null
/bin/rm -f  energy *~ *-e *neutral*   >& /dev/null
#
echo "> done..."
exit
