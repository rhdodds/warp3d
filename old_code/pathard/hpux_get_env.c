/*
      C-code to support POLO on HP-UX workstations.
      --------------------------------------------

      trgevn   -- get the value of an environment variable
*/

/*************************************************************
 *                                                           *
 *    Function to return contents of an environment          *
 *    variable as a character variable. specified string     *
 *    with variable name must be null terminated. if env.    *
 *    variable does not exist, the length (trgevn) is        *
 *    returned as zero, otherwise trgevn is number of        *
 *    characters in value of variable                        *
 *    Routine is in C but can be referenced as a single-     *
 *    precision function by F77                              *
 *                                                           *
 **************************************************************/


#include <stdio.h>
#include <stdlib.h>

int trgevn( env_name, env_value, len_name, len_value )

char *env_name, *env_value;
int len_name, len_value;
{
char *env_ptr;
char *getenv();
int  num_chars;
env_ptr = getenv( env_name);
num_chars = 0;
if ( *env_ptr != 0 )
 {
  num_chars = strlen( env_ptr );
  if ( num_chars <= len_value )
    {
      strcpy( env_value, env_ptr );
    }
   else
    {
      num_chars = -1;
    }
 }
 return num_chars;
}
