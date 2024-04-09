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
    std::string outFileNameHepMC = "truc.hepmc";
    std::string suffix = "";
    std::string mode = "hnl";
    std::string totalEvents = "100";

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
            continue;
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

    std::cout << "Infile: " << inFile <<std::endl;
    std::cout << "Outfile: " << outFileNameLHE <<std::endl;
    std::cout << "Suffix: " << suffix << std::endl;
    std::cout << "Mode: " << mode << std::endl;
    std::cout << "Total Events: " << totalEvents << std::endl;

    // Generator setup
    Pythia pythia;
    Event& event = pythia.event;
    // Read in commands from external file
    pythia.readFile(inFile);

    LHEF3FromPythia8 lhef3(&pythia.event, &pythia.info);

    // Open a file on which LHEF events should be stored, and write header
    if (outFileNameLHE=="")
        outFileNameLHE = inFile.replace(inFile.find(".cmnd"),std::string(".cmnd").length(),".lhe");
    
    // Checking the file extension
    if (outFileNameLHE.find(".") == std::string::npos)
        outFileNameLHE = outFileNameLHE + ".lhe";  // Add .lhe if no extension exists.  

    std::string fileExtension = "." + outFileNameLHE.substr(outFileNameLHE.find_last_of(".")+1); 
    if (fileExtension != ".lhe"){
        std::cout << "Current file extension is: " << fileExtension <<std::endl;
        outFileNameLHE = outFileNameLHE.replace(outFileNameLHE.find(fileExtension),fileExtension.length(),".lhe");
    }
    
    // Add suffix
    suffix = suffix + "_Events" + totalEvents;
    outFileNameLHE = outFileNameLHE.replace(outFileNameLHE.find(".lhe"),std::string(".lhe").length(),suffix+".lhe");



    lhef3.openLHEF(outFileNameLHE);

    // Initialize Pythia and set initial conditions
    pythia.readString("Main:numberOfEvents = " + totalEvents);
    int nEvent   = pythia.mode("Main:numberOfEvents");
    int nAbort   = pythia.mode("Main:timesAllowErrors");
    pythia.init();
    
    // LHEF output
    
    // Store initialization info in the LHAup object.
    lhef3.setInit();

    // Write out this initialization info on the file.
    lhef3.initLHEF();

    // HepMC output
    HepMC3::Pythia8ToHepMC3 toHepMC;
    HepMC3::WriterAscii hepMCWriter(outFileNameHepMC);

    

    int iAbort = 0;
    // Event loop
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
        if (!pythia.next()) {
        event.list();
        if (++iAbort < nAbort) continue;
        cout << " Event generation aborted prematurely, owing to error!\n";
        break;
        }

        // Write event to LHEF
        lhef3.setEvent();
        // lhef3.eventLHEF();
        // // Write event to HepMC
        HepMC3::GenEvent hepmcEvent;
        toHepMC.fill_next_event(pythia, hepmcEvent);
        hepMCWriter.write_event(hepmcEvent);
    }

    // Finalize and close files
    pythia.stat();
    // lhef3.updateSigma();
    lhef3.closeLHEF(true); // Update init block with final cross-section
    // HepMC file is closed automatically by its destructor

    return 0;
}
