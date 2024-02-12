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
import datetime

# Configuration de Pythia
pythia = pythia8.Pythia()
pythia.readFile("pythia_config.cmnd")  # Assurez-vous d'avoir le bon fichier de config
pythia.init()

# Préparer le fichier de sortie LHE
output_filename = "output.lhe"
with open(output_filename, "w") as file:
    now = datetime.datetime.now()
    file.write(f"<!--\n  File written by custom Pythia8 Python script on {now.strftime('%d %b %Y at %H:%M:%S')}\n-->\n")
    file.write("<LesHouchesEvents version=\"1.0\">\n")
    # Insérez ici les données d'initialisation
    file.write("<init>\n")
    file.write(" 2212  2212  6.500000e+03  6.500000e+03  0  0  0  0  3  1\n")
    file.write(" 8.218007e+04  3.901604e+02  1.000000e+00   9999\n")
    file.write("</init>\n")
    
# Nombre total d'événements à générer
    total_events = 10000  # Modifiez selon le nombre d'événements désiré

    for iEvent in range(total_events):
        if not pythia.next():
            continue

        # Commencer le bloc <event>
        file.write("<event>\n")
        # Le nombre de particules dans l'événement, l'ID du processus, etc.
        nParticles = pythia.event.size() - 2  # Exclure les particules du beam
        file.write(f"    {nParticles} 9999 1.000000e+00 WEIGHT SCALE PDF1 PDF2 ...\n")

        for i in range(pythia.event.size()):
            part = pythia.event[i]
            if part.isFinal() or part.idAbs() in [21, 22]:  # Filtrer pour les particules finales, gluons, et photons
                # Écrire les données de la particule
                file.write(f"    {part.id()} {part.status()} 0 0 0 0 {part.px()} {part.py()} {part.pz()} {part.e()} {part.m()} 0. 9.\n")

        file.write("</event>\n")
        
        
        # Finaliser le fichier LHE
    file.write("</LesHouchesEvents>\n")

# Statistiques Pythia
pythia.stat()
