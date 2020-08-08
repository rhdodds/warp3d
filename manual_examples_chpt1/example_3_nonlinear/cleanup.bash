#!/bin/bash
#
#           cleanup directory after running analysis
#
/bin/rm -f  model.text wef* wn* wem* wes* wnd* wne* wns* states_head*   >& /dev/null
/bin/rm -f  energy *~ *-e *neutral*  *exo z* *bak *pyc* *.pdf >& /dev/null
#
echo "> done..."
exit
