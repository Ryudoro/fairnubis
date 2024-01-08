#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

import os
import sys
from argparse import ArgumentParser

import shipunit as u

DownScaleDiMuon = False

# Default HNL parameters
theHNLMass   = 1.0*u.GeV
theProductionCouplings = theDecayCouplings = None

# Default dark photon parameters
theDPmass    = 0.2*u.GeV

# Alpaca
motherMode = True

mcEngine     = "TGeant4"
simEngine    = "Pythia8"  # "Genie" # Ntuple

inclusive    = "c"    # True = all processes if "c" only ccbar -> HNL, if "b" only bbar -> HNL, if "bc" only Bc+/Bc- -> HNL, and for darkphotons: if meson = production through meson decays, pbrem = proton bremstrahlung, qcd = ffbar -> DP.

MCTracksWithHitsOnly   = False  # copy particles which produced a hit and their history
MCTracksWithEnergyCutOnly = True # copy particles above a certain kin energy cut
MCTracksWithHitsOrEnergyCut = False # or of above, factor 2 file size increase compared to MCTracksWithEnergyCutOnly

charmonly    = False  # option to be set with -A to enable only charm decays, charm x-sec measurement  
HNL          = True

inputFile    = "/eos/experiment/ship/data/Charm/Cascade-parp16-MSTP82-1-MSEL4-978Bpot.root"
defaultInputFile = True

default = '2023'

inactivateMuonProcesses = False   # provisionally for making studies of various muon background sources

parser = ArgumentParser()
group = parser.add_mutually_exclusive_group()

parser.add_argument("--Pythia8", dest="pythia8", help="Use Pythia8", required=False, action="store_true")
parser.add_argument("-A",        dest="A",       help="b: signal from b, c: from c (default), bc: from Bc, or inclusive", required=False, default='c')
parser.add_argument("--Ntuple",  dest="ntuple",  help="Use ntuple as input", required=False, action="store_true")
parser.add_argument("--MuonBack",dest="muonback",  help="Generate events from muon background file, --Cosmics=0 for cosmic generator data", required=False, action="store_true")
parser.add_argument("--Cosmics",   dest="cosmics",  help="Use cosmic generator, argument switch for cosmic generator 0 or 1", required=False,  default=None)
parser.add_argument("--MuDIS",     dest="mudis",  help="Use muon deep inelastic scattering generator", required=False, action="store_true")
parser.add_argument("--RpvSusy", dest="RPVSUSY",  help="Generate events based on RPV neutralino", required=False, action="store_true")
parser.add_argument("--DarkPhoton", dest="DarkPhoton",  help="Generate dark photons", required=False, action="store_true")
parser.add_argument("--SusyBench", dest="RPVSUSYbench",  help="Generate HP Susy", required=False, default=2)
parser.add_argument("-m", "--mass", dest="theMass",  help="Mass of hidden particle, default "+str(theHNLMass)+"GeV for HNL, "+str(theDPmass)+"GeV for DP", required=False, default=None, type=float)
parser.add_argument("-c", "--couplings", "--coupling", dest="thecouplings",  help="couplings \'U2e,U2mu,U2tau\' or -c \'U2e,U2mu,U2tau\' to set list of HNL couplings.\
 TP default for HNL, ctau=53.3km", required=False,default="0.447e-9,7.15e-9,1.88e-9")
parser.add_argument("-cp", "--production-couplings", dest="theprodcouplings",  help="production couplings \'U2e,U2mu,U2tau\' to set the couplings for HNL production only"\
                                            ,required=False,default=None)
parser.add_argument("-cd", "--decay-couplings", dest="thedeccouplings",  help="decay couplings  \'U2e,U2mu,U2tau\' to set the couplings for HNL decay only", required=False,default=None)
parser.add_argument("-e", "--epsilon", dest="theDPepsilon",  help="to set mixing parameter epsilon", required=False,default=0.00000008, type=float)
parser.add_argument("-n", "--nEvents",dest="nEvents",  help="Number of events to generate", required=False,  default=100, type=int)
parser.add_argument("-i", "--firstEvent",dest="firstEvent",  help="First event of input file to use", required=False,  default=0, type=int)
parser.add_argument("-s", "--seed",dest="theSeed",  help="Seed for random number. Only for experts, see TRrandom::SetSeed documentation", required=False,  default=0, type=int)
parser.add_argument("-S", "--sameSeed",dest="sameSeed",  help="can be set to an integer for the muonBackground simulation with specific seed for each muon, only for experts!"\
                                            ,required=False,  default=False, type=int)
group.add_argument("-f",        dest="inputFile",       help="Input file if not default file", required=False, default=False)
parser.add_argument("-o", "--output",dest="outputDir",  help="Output directory", required=False,  default=".")
parser.add_argument("-t", "--test", dest="testFlag",  help="quick test", required=False,action="store_true")
parser.add_argument("--dry-run", dest="dryrun",  help="stop after initialize", required=False,action="store_true")
parser.add_argument("--scName", help="The name of the SC shield in the database", default="sc_v6")
parser.add_argument("--MesonMother",   dest="MM",  help="Choose DP production meson source", required=False,  default=True)
parser.add_argument("--debug",  help="1: print weights and field 2: make overlap check", required=False, default=0, type=int, choices=range(0,3))

options = parser.parse_args()

if options.A != 'c':
     inclusive = options.A
     if options.A =='b': inputFile = "/eos/experiment/ship/data/Beauty/Cascade-run0-19-parp16-MSTP82-1-MSEL5-5338Bpot.root"
     if options.A.lower() == 'charmonly':
           charmonly = True
           HNL = False 
     if options.A not in ['b','c','bc','meson','pbrem','qcd']: inclusive = True
     
if options.inputFile:
  if options.inputFile == "none": options.inputFile = None
  inputFile = options.inputFile
  defaultInputFile = False
if options.RPVSUSY: HNL = False
if options.DarkPhoton: HNL = False
if not options.theMass:
  if options.DarkPhoton: options.theMass  = theDPmass
  else:                  options.theMass  = theHNLMass
if options.thecouplings:
  theCouplings = [float(c) for c in options.thecouplings.split(",")]
if options.theprodcouplings:
  theProductionCouplings = [float(c) for c in options.theprodcouplings.split(",")]
if options.thedeccouplings:
  theDecayCouplings = [float(c) for c in options.thedeccouplings.split(",")]
  
#sanity check
if (HNL and options.RPVSUSY) or (HNL and options.DarkPhoton) or (options.DarkPhoton and options.RPVSUSY): 
 print("cannot have HNL and SUSY or DP at the same time, abort")
 sys.exit(2)
 
if not os.path.exists(options.outputDir):
  os.makedirs(options.outputDir)
  
