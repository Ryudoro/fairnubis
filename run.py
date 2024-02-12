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
#Setup the output file name
outFile=args.inFile.replace(".cmnd",f"_{suffix}.lhe")
# outputLH = pythia8.LHAupFromPYTHIA8(pythia.process, pythia.info)
# outputLH.openLHEF(outFile)

# Read in the user-specified cmnd file with the pythia configuration
pythia.readFile(args.inFile)

# Set the number of events in this file directly - allowing more events to be created
pythia.readString(f"Main:numberOfEvents = {args.events}")

# Error Checking number of events
nEvents=pythia.mode("Main:numberOfEvents")
if args.events !=nEvents:
	raise Exception(f"Number of events in pythia do not match requested input: {nEvents} =/= {args.events}")

 # Initialise the relevant Settings and ParticleData for pythia - store and write to LHE file
pythia.init()
# outputLH.setInit()
# outputLH.initLHEF()

# # Generate pythia events one at a time
# nAbort=0
# maxErr=pythia.mode("Main:timesAllowErrors")+10
# for event in range(nEvents):
# 	# Generate events with pythia.next() - will output bool with 1 for success and 0 for a failure
# 	if not pythia.next():
# 		nAbort+=1 
# 		# Break event loop if maximum number of errors reached	
# 		if nAbort >= maxErr:
# 			print(f"Event generation ended prematurely at Event {event}, due to too many errors ({maxErr})")
# 			break
# 		continue 

	# Store event info in LHAup object and write to file
	# outputLH.setEvent()
	# outputLH.eventLHEF(verbose=True) # Can set verbose false for smaller file
	#----------------------------------------------------------------------------------------------
	# Event-by-event processing can occur here - e.g. filling kinematic histograms, cutting events. 

	#----------------------------------------------------------------------------------------------

# Get the final statistical information about the pythia events generated
pythia.stat()

# # Update X-section info base on MC integration during run
# outputLH.updateSigma()

# # Write and close the LH file - set to true to overwrite initialisation w/ new X-sections
# outputLH.closeLHEF(True)
from collections import Counter
from collections import Counter

with open('pythia_config.cmnd', 'r') as f:
    if '9900015' in f.readlines():
        hnl = 23
    else : 
        hnl = 9900015
# Préparation pour le bar plot des particules mères des HNLs
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

    # Parcourir toutes les particules de l'événement
    for i in range(pythia.event.size()):
        particle = pythia.event[i]
        if particle.id() == hnl:  # Filtrer pour trouver vos HNLs
            # Collecter les ID des particules mères
            if particle.mother1() != 0:
                mother_id = pythia.event[particle.mother1()].id()
                hnl_mother_ids.append(mother_id)
            if particle.mother2() != 0 and particle.mother2() != particle.mother1():
                mother_id = pythia.event[particle.mother2()].id()
                hnl_mother_ids.append(mother_id)

# Compter les occurrences de chaque ID de particule mère
mother_counts = Counter(hnl_mother_ids)
mothers = list(mother_counts.keys())
counts = list(mother_counts.values())
print(mother_counts,  mothers)
# Utiliser les indices des particules mères comme positions des barres sur l'axe des x
positions = range(len(mothers))

# # Créer un bar plot
# plt.figure(figsize=(12, 6))
# plt.bar(positions, counts, align='center', alpha=0.7)
# plt.xlabel('ID des particules mères')
# plt.ylabel('Nombre de HNLs')
# plt.title('Bar plot des particules mères des HNLs')
# plt.xticks(positions, mothers, rotation=45)
# plt.show()


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

    # Parcourir toutes les particules de l'événement
    for i in range(pythia.event.size()):
        particle = pythia.event[i]
        if particle.id() == hnl:  # Filtrer pour trouver vos HNLs
            print(particle.e())
            hnl_energies.append(particle.e())  # Ajouter l'énergie de la particule HNL à la liste
            hnl_etas.append(particle.eta())
            hnl_phis.append(particle.phi())

            # Calculer la distance depuis le point d'interaction en direction R
            px, py = particle.px(), particle.py()
            distance_R = math.sqrt(px**2 + py**2)  # Distance dans le plan transverse
            hnl_distances.append(distance_R)

# Affichage de l'histogramme des énergies HNL
plt.figure()
plt.hist(hnl_energies, bins=50, alpha=0.7, label='Énergie des HNL')
plt.xlabel('Énergie [GeV]')
plt.ylabel('Nombre de particules')
plt.title('Histogramme de l\'énergie des particules HNL')
plt.legend()
plt.show()

plt.figure()
plt.hist(hnl_etas, bins=50, alpha=0.7, label='Eta des HNL')
plt.xlabel('Eta')
plt.ylabel('Nombre de particules')
plt.title('Histogramme de l\'eta des particules HNL')
plt.legend()
plt.show()


plt.figure()
plt.hist(hnl_phis, bins=50, alpha=0.7, label='Phi des HNL')
plt.xlabel('Phi')
plt.ylabel('Nombre de particules')
plt.title('Histogramme du phi des particules HNL')
plt.legend()
plt.show()

plt.figure()
plt.hist(hnl_distances, bins=50, alpha=0.7, label='Distance R depuis IP')
plt.xlabel('Distance R [m]')
plt.ylabel('Nombre de particules')
plt.title('Histogramme de la distance R depuis le point d\'interaction pour les HNL')
plt.legend()
plt.show()



