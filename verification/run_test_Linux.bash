#!/bin/bash
#
#
#    for Linux, set environment varaible to inlude libraries
#    in the WARP3D distribution.
#
#    not needed for Windows or Macs
#
is_linux=`uname | grep Linux | wc -l`
#echo $is_linux
if [[ $is_linux -eq 0 ]]
then
  echo " "
  echo ".... This script only for LINUX ...."
  echo " "
  exit
fi
export LD_LIBRARY_PATH=$WARP3D_HOME/linux_packages/lib:$LD_LIBRARY_PATH
#
python3 run_tests.py
#
exit
