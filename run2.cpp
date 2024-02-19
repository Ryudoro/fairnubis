/* 
   Pythia runner script to produce Les Houches Event Files in LHEF 3.0 format,
   including detailed information like momentum, angles, masses, etc.
   Uses LHEF3FromPYTHIA8 for comprehensive event output.
*/

#include "Pythia8/Pythia.h"

//#include "Pythia8Plugins/LHEF3.h"

// #include "Pythia8Plugins/HepMC3.h"
#include <iostream>
#include <regex>
#include <string>
using namespace Pythia8;

int main(int argc, char *argv[]) {
    // Setup argument parsing: including defaults
    std::string inFile = "pythia_config.cmnd";
    std::string outFileName = "";
    std::string suffix = "";
    std::string mode = "hnl";
    std::string totalEvents = "10000";

    for (int i = 0; i < argc; i++) {
        if (std::string(argv[i]) == "--infile" || std::string(argv[i]) == "-i") {
            inFile = std::string(argv[i + 1]);
            i++;
            continue;
        } else if (std::string(argv[i]) == "--outfile" || std::string(argv[i]) == "-o") {
            outFileName = std::string(argv[i + 1]);
            i++;
            continue;
        } else if (std::string(argv[i]) == "--suffix" || std::string(argv[i]) == "-s") {
            suffix = std::string(argv[i + 1]);
            i++;
            continue;
        } else if (std::string(argv[i]) == "--mode" || std::string(argv[i]) == "-m") {
            mode = std::string(argv[i + 1]);
            i++;
            continue;
        } else if (std::string(argv[i]) == "--nevents" || std::string(argv[i]) == "-n") {
            totalEvents = std::string(argv[i + 1]);
            if (!std::regex_match(totalEvents, std::regex("((\\+|-)?[[:digit:]]+)(\\.(([[:digit:]]+)?)?"))) {
                std::cout << "Input --nevents/-n (" << totalEvents << ") is not a number, please enter an integer value." << std::endl;
                return 1;
            }
            totalEvents = std::to_string(std::stoi(totalEvents));
            i++;
            continue;
        } else if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h") {
            std::cout << "Run the pythia event generation. Usage:" << std::endl;
            std::cout << std::endl;
            std::cout << "Run [options]" << std::endl;
            std::cout << std::endl;
            std::cout << "--nevents/-n: The total number of pythia events to simulate." << std::endl;
            std::cout << "--infile/-i:  The input cmnd configuration file to use." << std::endl;
            std::cout << "--outfile/-o: The name of the Les Houches file output." << std::endl;
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
    // pythia.readString("Beams:frameType = 4");
    // pythia.readString("Beams:LHEF = wbj_lhef3.lhe");
    // Create an LHEF3 object for output
    LHEF3FromPythia8 myLHEF3(&pythia.event,&pythia.info);


    // Prepare the output file name
    if (outFileName.empty())
        outFileName = inFile.replace(inFile.find(".cmnd"), std::string(".cmnd").length(), ".lhe");
    if (outFileName.find(".") == std::string::npos)
        outFileName += ".lhe"; // Ensure the file has .lhe extension
    std::string fileExtension = "." + outFileName.substr(outFileName.find_last_of(".") + 1);
    if (fileExtension != ".lhe")
        outFileName.replace(outFileName.find(fileExtension), fileExtension.length(), ".lhe");
    suffix = suffix + "_Events" + totalEvents;
    outFileName.replace(outFileName.find(".lhe"), std::string(".lhe").length(), suffix + ".lhe");

    // Open LHEF file for output
    myLHEF3.openLHEF(outFileName);

    // Initialize Pythia and set initial conditions in LHEF
    pythia.readString("Main:numberOfEvents = " + totalEvents);
    pythia.init();
    myLHEF3.setInit();

    // Begin event loop
    for (int iEvent = 0; iEvent < std::stoi(totalEvents); ++iEvent) {
        if (!pythia.next()) continue;

        // Write event to LHEF
        myLHEF3.setEvent();
    }

    // Finalize and close LHEF file
    pythia.stat();
    myLHEF3.closeLHEF(true); // Update init block with final cross-section

    return 0;
}
