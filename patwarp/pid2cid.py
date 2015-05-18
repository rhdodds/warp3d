#!/usr/bin/python

'''
Takes a patran neutral file and copies the PID in the element (02) card to the config ID
'''

import sys

# Check for number of CLA (2 = script name and file name)
if (len(sys.argv) != 2):
      raise ValueError("Invalid number of command line arguments.  Requires exactly one - the neutral file name.")

# Get the file name
filename = sys.argv[1]

# Open up our file, loop, writing to stdout, until we find a 02 card, do the switch on those
with open(filename) as o_file:
      packet_line = 0
      twopacket = False
      for line in o_file:
#           If this is zero we have an actual packet
            if packet_line == 0:
                  toks = line.split()
                  if int(toks[0]) == 2:
                        twopacket = True
                  else:
                        twopacket = False
                  tlines = int(toks[3])+1
                  packet_line = tlines
            
#           The line we want to edit, copy token 2 to 1 (zero index)
            if ( twopacket and (tlines-packet_line) == 1):
                  toks = line.split()
                  toks[1] = toks[2]
                  line = (' '*(8-len(str(toks[0]))))+toks[0]+(' '*(8-len(str(toks[1]))))+toks[1]+(' '*(8-len(str(toks[2]))))+toks[2]+(' '*(8-len(str(toks[3]))))+toks[3]+(' '*(16-len(str(toks[4]))))+toks[4]+(' '*(16-len(str(toks[5]))))+toks[5]+(' '*(16-len(str(toks[6]))))+toks[6]+'\n'

            print(line,end='')
            packet_line-=1

