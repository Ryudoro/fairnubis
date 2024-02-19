/* 
   Pythia runner script to produce Les Houches Files from inputs of cmnd files produced by the python scripts.
   Makes use of code adapted from Pythia examples: - main01.cc, main20.cc, main76.cc 
 */

#include "Pythia8/Pythia.h"
#include <iostream>
#include <regex>
#include <string>
using namespace Pythia8;
int main(int argc, char *argv[]) {
  // Setup argument parsing: including defaults
  std::string inFile="pythia_config.cmnd";
  std::string outFileName="";
  std::string suffix="";
  std::string mode="hnl";
  std::string totalEvents = "10000";

  for(int i=0; i<argc; i++){
      if(std::string(argv[i]) == "--infile" || std::string(argv[i]) == "-i"){ 
	      inFile = std::string(argv[i+1]);
	      i++;
	      continue;
      }
      else if(std::string(argv[i]) == "--outfile" || std::string(argv[i]) == "-o"){ 
	      outFileName = std::string(argv[i+1]);
	      i++;
	      continue;
      }
      else if(std::string(argv[i]) == "--suffix" || std::string(argv[i]) == "-s"){ 
	      suffix = std::string(argv[i+1]);
	      i++;
	      continue;
      }
      else if(std::string(argv[i]) == "--mode" || std::string(argv[i]) == "-m"){ 
	      mode = std::string(argv[i+1]);
	      i++;
	      continue;
      }
      else if(std::string(argv[i]) == "--nevents" || std::string(argv[i]) == "-n"){ 
	      totalEvents = std::string(argv[i+1]);
	      
	      // Check totalEvents is an integer
	      if (!std::regex_match( totalEvents, std::regex( ( "((\\+|-)?[[:digit:]]+)(\\.(([[:digit:]]+)?))?" ) ) )){
	      	std::cout << "Input --nevents/-n (" << totalEvents << ") is not a number, please enter an integer value." <<std::endl;
		return 1;
	      }

	      totalEvents = std::to_string( std::stoi(totalEvents) );

	      i++;
	      continue;
      }
      else if(std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h"){ 
	      std::cout << "Run the pythia event generation. Usage:" << std::endl;
	      std::cout<< std::endl;
	      std::cout<< "Run [options]" << std::endl;
	      std::cout<< std::endl;
	      std::cout<< "--nevents/-n: The total number of pythia events to simulate." <<std::endl;
	      std::cout<< "--infile/-i:  The input cmnd configuration file to use." <<std::endl;
	      std::cout<< "--outfile/-o: The name of the Les Houches file output." <<std::endl;
	      std::cout<< "--suffix/-s:  A suffix to add to the output file name (if not specified) and any produced plots." <<std::endl;
	      std::cout<< "--mode/-m:    Specify which type of LLP model to simulate (Unused currently)." <<std::endl;
	      std::cout<< "--help/-h:    Print this Message!." <<std::endl;
	      return 0;
      }
  }

  std::cout << "Infile: " << inFile <<std::endl;
  std::cout << "Outfile: " << outFileName <<std::endl;
  std::cout << "Suffix: " << suffix << std::endl;
  std::cout << "Mode: " << mode << std::endl;
  std::cout << "Total Events: " << totalEvents << std::endl;

  // Generator.
  Pythia pythia;
  Event& event = pythia.event;

  // Read in commands from external file.
  //pythia.readFile("pythia_config.cmnd");
  pythia.readFile(inFile);

  // Create an LHAup object that can access relevant information in pythia.
  LHAupFromPYTHIA8 myLHA(&pythia.process, &pythia.info);

  // Open a file on which LHEF events should be stored, and write header
  if (outFileName=="")
	outFileName = inFile.replace(inFile.find(".cmnd"),std::string(".cmnd").length(),".lhe");
  
  // Checking the file extension
  if (outFileName.find(".") == std::string::npos)
	outFileName = outFileName + ".lhe";  // Add .lhe if no extension exists.  

  std::string fileExtension = "." + outFileName.substr(outFileName.find_last_of(".")+1); 
  if (fileExtension != ".lhe"){
	std::cout << "Current file extension is: " << fileExtension <<std::endl;
	outFileName = outFileName.replace(outFileName.find(fileExtension),fileExtension.length(),".lhe");
  }
 
  // Add suffix
  suffix = suffix + "_Events" + totalEvents;
  outFileName = outFileName.replace(outFileName.find(".lhe"),std::string(".lhe").length(),suffix+".lhe");

  //myLHA.openLHEF("test_output.lhe");
  myLHA.openLHEF(outFileName);

  pythia.readString("Main:numberOfEvents = " + totalEvents); 

  // Extract settings to be used in the main program.
  int nEvent   = pythia.mode("Main:numberOfEvents");
  int nAbort   = pythia.mode("Main:timesAllowErrors");

  // Initialize. 
  pythia.init();

  // Store initialization info in the LHAup object.
  myLHA.setInit();

  // Write out this initialization info on the file.
  myLHA.initLHEF();


  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      event.list();
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }


    // Perform any analysis within the loop itself


    //-------- Les Houches File Outputting

    // Store event info in the LHAup object.
    myLHA.setEvent();

    // Write out this event info on the file.
    // With optional argument (verbose =) false the file is smaller.
    myLHA.eventLHEF();

  // End of event loop.
  }

  // Final statistics and histogram output.
  pythia.stat();

  // Update the cross section info based on Monte Carlo integration during run.
  myLHA.updateSigma();

  // Write endtag. Overwrite initialization info with new cross sections.
  myLHA.closeLHEF(true);

  return 0;
}

