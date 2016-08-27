/** \file PINProcessor.cpp
 * PIN class takes the signals from the mcp detector
 *   and calculates a 2D position based on the readout
 * SNL 2-2-08, Modified DTM 9-09
 */
#include "DammPlotIds.hpp"
#include "PINProcessor.hpp"
#include "RawEvent.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>

using std::string;
using std::vector;

using namespace dammIds::pin;

namespace dammIds {// offset = 950; range = 10;
  namespace pin {      
    const int DD_POSENE = 0;
    const int D_ENE     = 1;// ID=951; YX
  }
}
void PINProcessor::PINData::Clear(void)
{
}

PINProcessor::PINProcessor(void) : EventProcessor(OFFSET, RANGE, "pin")
{
  associatedTypes.insert("pin");
}

void PINProcessor::DeclarePlots(void)
{
  
  const int posBins      = S4; //    16
  const int energyBins   = SE; // 16384
  DeclareHistogram2D(DD_POSENE, posBins, energyBins, "PIN Position and Energy");
  DeclareHistogram1D(D_ENE,energyBins,"Pin1"); // by YX
  DeclareHistogram1D(D_ENE+1,energyBins,"Pin2");
  
}

//std::ofstream outfile;

bool PINProcessor::Process(RawEvent &event)
{
  if (!EventProcessor::Process(event))
    return false;

  static const vector<ChanEvent*> &pinEvents = sumMap["pin"]->GetList();

  data.Clear();
  for (vector<ChanEvent*>::const_iterator it = pinEvents.begin();
       it != pinEvents.end(); it++) {
      ChanEvent *chan = *it;

      string subtype   = chan->GetChanID().GetSubtype();
      int number 	= chan->GetChanID().GetLocation();
      double calEnergy = chan->GetCalEnergy(); // should be the same with StripEvent.E; by YX
      double pinTime   = chan->GetTime(); // count at the same time with DSSD? by YX
      using std::cout;
      using std:: endl; 
      plot(DD_POSENE,number,calEnergy);

      //      cout << number << endl;
      if(number==0){
	plot(D_ENE,calEnergy); 
      }
      if(number==2){
	//------------- by Yongchi; 03/12/2015
	//std::ofstream outfile;
	/*
	outfile.open("resultPIN2.txt", std::ofstream::out | std::ofstream::app);
        outfile << std::setw(10) << calEnergy 
		<< std::setw(14) << pinTime
		<< endl;
	outfile.close();
	*/
	//-----------------------------------
      plot(D_ENE+1,calEnergy);
      }
  }
      
      
      EndProcess();
      return true;
}


