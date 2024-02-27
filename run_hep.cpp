/* 
   Pythia runner script to produce both Les Houches Event Files in LHEF 3.0 format
   and HepMC event files.
   Uses LHEF3FromPYTHIA8 for LHEF output and HepMC3 interface for HepMC output.
*/

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"
#include <iostream>
#include <fstream>
#include <string>
#include <regex>
using namespace Pythia8;

int main(int argc, char *argv[]) {
    // Setup argument parsing: including defaults
    std::string inFile = "pythia_config.cmnd";
    std::string outFileNameLHE = "";
    std::string outFileNameHepMC = "";
    std::string suffix = "";
    std::string mode = "hnl";
    std::string totalEvents = "10000";

    // Extended argument parsing to include HepMC output file
    for (int i = 0; i < argc; i++) {
        if (std::string(argv[i]) == "--infile" || std::string(argv[i]) == "-i") {
            inFile = std::string(argv[i + 1]);
            i++;
        } else if (std::string(argv[i]) == "--outfileLHE" || std::string(argv[i]) == "-ol") {
            outFileNameLHE = std::string(argv[i + 1]);
            i++;
        } else if (std::string(argv[i]) == "--outfileHepMC" || std::string(argv[i]) == "-oh") {
            outFileNameHepMC = std::string(argv[i + 1]);
            i++;
        } else if (std::string(argv[i]) == "--suffix" || std::string(argv[i]) == "-s") {
            suffix = std::string(argv[i + 1]);
            i++;
        } else if (std::string(argv[i]) == "--mode" || std::string(argv[i]) == "-m") {
            mode = std::string(argv[i + 1]);
            i++;
        } else if (std::string(argv[i]) == "--nevents" || std::string(argv[i]) == "-n") {
            totalEvents = std::string(argv[i + 1]);
            if (!std::regex_match(totalEvents, std::regex("((\\+|-)?[[:digit:]]+)(\\.(([[:digit:]]+)?)?"))) {
                std::cout << "Input --nevents/-n (" << totalEvents << ") is not a number, please enter an integer value." << std::endl;
                return 1;
            }
            totalEvents = std::to_string(std::stoi(totalEvents));
            i++;
        } else if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h") {
            std::cout << "Run the pythia event generation. Usage:" << std::endl;
            std::cout << std::endl;
            std::cout << "Run [options]" << std::endl;
            std::cout << std::endl;
            std::cout << "--nevents/-n: The total number of pythia events to simulate." << std::endl;
            std::cout << "--infile/-i:  The input cmnd configuration file to use." << std::endl;
            std::cout << "--outfileLHE/-ol: The name of the Les Houches file output." << std::endl;
            std::cout << "--outfileHepMC/-oh: The name of the HepMC file output." << std::endl;
            std::cout << "--suffix/-s:  A suffix to add to the output file name." << std::endl;
            std::cout << "--mode/-m:    Specify which type of model to simulate." << std::endl;
            std::cout << "--help/-h:    Print this Message!." << std::endl;
            return 0;
        }
    }

    // Generator setup
    Pythia pythia;
    // Read in commands from external file
    pythia.readFile(inFile);

    // Prepare the output file names
    if (outFileNameLHE.empty())
        outFileNameLHE = inFile.replace(inFile.find(".cmnd"), std::string(".cmnd").length(), "_lhe.lhe");
    if (outFileNameHepMC.empty())
        outFileNameHepMC = inFile.replace(inFile.find(".cmnd"), std::string(".cmnd").length(), "_hepmc.hepmc");

    // Append suffix to file names if provided
    if (!suffix.empty()) {
        outFileNameLHE.insert(outFileNameLHE.find(".lhe"), "_" + suffix);
        outFileNameHepMC.insert(outFileNameHepMC.find(".hepmc"), "_" + suffix);
    }

    // Initialize Pythia and set initial conditions
    pythia.readString("Main:numberOfEvents = " + totalEvents);
    pythia.init();

    // LHEF output
    LHEF3FromPythia8 lhef3(pythia.event, pythia.info);
    lhef3.openLHEF(outFileNameLHE);
    lhef3.setInit();

    // HepMC output
    HepMC3::Pythia8ToHepMC toHepMC;
    HepMC3::WriterAscii hepMCWriter(outFileNameHepMC);

    // Event loop
    for (int iEvent = 0; iEvent < std::stoi(totalEvents); ++iEvent) {
        if (!pythia.next()) continue;

        // Write event to LHEF
        lhef3.setEvent();

        // Write event to HepMC
        HepMC3::GenEvent hepmcEvent;
        toHepMC.fill_next_event(pythia, hepmcEvent);
        hepMCWriter.write_event(hepmcEvent);
    }

    // Finalize and close files
    pythia.stat();
    lhef3.closeLHEF(true); // Update init block with final cross-section
    // HepMC file is closed automatically by its destructor

    return 0;
}
