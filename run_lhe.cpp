#include "Pythia8/Pythia.h"

using namespace Pythia8;

#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char* argv[]) {
    // Vérifier les arguments
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <commandFile> <outputFile> <nEvents>" << std::endl;
        return 1;
    }
    std::string commandFile = argv[1];
    std::string outputFile = argv[2];
    int nEvents = std::stoi(argv[3]);

    // Créer une instance de Pythia
    Pythia pythia;

    // Lire le fichier de commande
    pythia.readFile(commandFile);

    // Définir le nombre d'événements
    pythia.readString("Main:numberOfEvents = " + std::to_string(nEvents));

    // Initialiser Pythia
    pythia.init();

    // Créer une instance de LHAup pour écrire le fichier LHE
    LHAupFromPYTHIA8 up(pythia.event, pythia.info);
    up.openLHEF(outputFile);

    // Générer les événements
    int nAbort = 0;
    int maxErr = pythia.mode("Main:timesAllowErrors") + 10;
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        if (!pythia.next()) {
            nAbort++;
            if (nAbort >= maxErr) {
                std::cerr << "Event generation ended prematurely at Event " << iEvent
                          << ", due to too many errors (" << maxErr << ")" << std::endl;
                break;
            }
            continue;
        }

        // Écrire l'événement dans le fichier LHE
        up.setEvent();
        up.eventLHEF();

        // Ici, vous pouvez ajouter du code pour traiter chaque événement si nécessaire
    }

    // Imprimer les statistiques des événements générés
    pythia.stat();

    // Fermer le fichier LHE
    up.closeLHEF(true);

    return 0;
}
