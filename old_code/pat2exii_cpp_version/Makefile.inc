# Makefile definitions for patran to EXODUSII converter

# Programs
CXX = g++
LD = g++

# Flags
CXXFLAGS = -c -g -O2

# Includes
INCLUDE = -I../include

# Libraries
LIBS = ../lib/libnetcdf.a ../lib/libexoIIv2c.a ../lib/libnetcdf.a ../lib/libhdf5_hl.a ../lib/libhdf5.a ../lib/libsz.a -lcurl -lz

# Other
EXE = pat2exII.linux
INSTDIR = ../run_linux_em64t
