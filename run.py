import sys, os
import matplotlib.pyplot as plt
import math

print(os.getcwd())
cfg = open("pythia8310/Makefile.inc")
lib = "../lib"
for line in cfg:
    if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
sys.path.insert(0, lib)
print(lib)
import argparse
import pythia8

parser=argparse.ArgumentParser()
parser.add_argument('-N','--events', type=int, help="Specify the number of events to generate", default=10000)
parser.add_argument('--mode', choices=["hnl"], type=str, help="Specify what type of LLP to generate", default="hnl")
parser.add_argument('-i','--inFile', type=str, help="Specify an input cmnd file", default="pythia_config.cmnd", required=False)
parser.add_argument('--suffix', type=str, help="Specify a suffix for the output Les Houches file", default="")
args=parser.parse_args()


pythia = pythia8.Pythia()
print(dir(pythia))
suffix = "output"


outFile=args.inFile.replace(".cmnd",f"_{suffix}.lhe")

pythia.readFile(args.inFile)

pythia.readString(f"Main:numberOfEvents = {args.events}")

nEvents=pythia.mode("Main:numberOfEvents")
if args.events !=nEvents:
	raise Exception(f"Number of events in pythia do not match requested input: {nEvents} =/= {args.events}")

pythia.init()

pythia.stat()

from collections import Counter
from collections import Counter

with open('pythia_config.cmnd', 'r') as f:
    if '9900015' in f.readlines():
        hnl = 9900015
    else : 
        hnl = 9900015

hnl_mother_ids = []
hnl_etas = []
hnl_phis = []
hnl_distances = []
nAbort = 0
maxErr = pythia.mode("Main:timesAllowErrors") + 10
for event in range(nEvents):
    if not pythia.next():
        nAbort += 1
        if nAbort >= maxErr:
            print(f"Event generation ended prematurely at Event {event}, due to too many errors ({maxErr})")
            break
        continue

    for i in range(pythia.event.size()):
        particle = pythia.event[i]
        if particle.id() == hnl:  

            if particle.mother1() != 0:
                mother_id = pythia.event[particle.mother1()].id()
                hnl_mother_ids.append(mother_id)
            if particle.mother2() != 0 and particle.mother2() != particle.mother1():
                mother_id = pythia.event[particle.mother2()].id()
                hnl_mother_ids.append(mother_id)


mother_counts = Counter(hnl_mother_ids)
mothers = list(mother_counts.keys())
counts = list(mother_counts.values())
print(mother_counts,  mothers)

positions = range(len(mothers))

hnl_energies = []

nAbort = 0
maxErr = pythia.mode("Main:timesAllowErrors") + 10
for event in range(nEvents):
    if not pythia.next():
        nAbort += 1
        if nAbort >= maxErr:
            print(f"Event generation ended prematurely at Event {event}, due to too many errors ({maxErr})")
            break
        continue

    for i in range(pythia.event.size()):
        particle = pythia.event[i]
        if particle.id() == hnl: 
            print(particle.e())
            hnl_energies.append(particle.e())
            hnl_etas.append(particle.eta())
            hnl_phis.append(particle.phi())

            px, py = particle.px(), particle.py()
            distance_R = math.sqrt(px**2 + py**2) 
            hnl_distances.append(distance_R)

# Plotting HNL energy histogram
plt.figure()
plt.hist(hnl_energies, bins=50, alpha=0.7, label='HNL Energy')
plt.xlabel('Energy [GeV]')
plt.ylabel('Number of Particles')
plt.title('Histogram of HNL Particle Energies')
plt.legend()
plt.show()

# Plotting HNL eta histogram
plt.figure()
plt.hist(hnl_etas, bins=50, alpha=0.7, label='HNL Eta')
plt.xlabel('Eta')
plt.ylabel('Number of Particles')
plt.title('Histogram of HNL Particle Eta')
plt.legend()
plt.show()

# Plotting HNL phi histogram
plt.figure()
plt.hist(hnl_phis, bins=50, alpha=0.7, label='HNL Phi')
plt.xlabel('Phi')
plt.ylabel('Number of Particles')
plt.title('Histogram of HNL Particle Phi')
plt.legend()
plt.show()

# Plotting HNL distance R from interaction point histogram
plt.figure()
plt.hist(hnl_distances, bins=50, alpha=0.7, label='Distance R from IP')
plt.xlabel('Distance R [m]')
plt.ylabel('Number of Particles')
plt.title('Histogram of Distance R from Interaction Point for HNLs')
plt.legend()
plt.show()



