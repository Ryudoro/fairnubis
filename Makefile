# Set the shell.
SHELL=/usr/bin/env bash

# Install directory prefixes.
PREFIX_BIN=pythia8310/bin
PREFIX_INCLUDE=pythia8310/include
PREFIX_LIB=pythia8310/lib
PREFIX_SHARE=pythia8310/share/Pythia8

# Compilation flags 
CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -I$(PREFIX_INCLUDE) $(GZIP_LIB)
CXX_COMMON+= -L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -lpythia8 -ldl
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so

PYTHIA=$(PREFIX_LIB)/libpythia8$(LIB_SUFFIX)


#-----------------------#
#-- Making the runner --#
#-----------------------#

all: 
	make run

run: $(PYTHIA) run.cc
	${CXX} $@.cc -o $@ $(CXX_COMMON) 

clean: 
	rm run
