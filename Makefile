# Set the shell.
SHELL=/usr/bin/env bash

# Build directory.
BUILD_DIR=build

# Install directory prefixes.
PREFIX_BIN=pythia8/bin
PREFIX_INCLUDE=pythia8/include
PREFIX_LIB=pythia8/lib
PREFIX_SHARE=pythia8/share/Pythia8

# HepMC3 install directory prefix.
HEPMC3_INSTALL_DIR=hepmc3-install

# Compilation flags 
CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -I$(PREFIX_INCLUDE) -I$(HEPMC3_INSTALL_DIR)/include $(GZIP_LIB)
CXX_COMMON+= -L$(PREFIX_LIB) -L$(HEPMC3_INSTALL_DIR)/lib -Wl,-rpath,$(PREFIX_LIB):$(HEPMC3_INSTALL_DIR)/lib -lpythia8 -ldl -lHepMC3
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so

PYTHIA=$(PREFIX_LIB)/libpythia8$(LIB_SUFFIX)

#-----------------------#
#-- Making the runner --#
#-----------------------#

all: $(BUILD_DIR)/run

$(BUILD_DIR)/run: $(PYTHIA) run_hep.cpp
	@mkdir -p $(BUILD_DIR)
	${CXX} run_hep.cpp -o $@ $(CXX_COMMON) 

clean: 
	rm -rf $(BUILD_DIR)

