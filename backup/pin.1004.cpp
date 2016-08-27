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
  DeclareHistogram1D(D_ENE,energyBins,"Pin1"); // by YX; 951
  DeclareHistogram1D(D_ENE+1,energyBins,"Pin2"); // 952
  DeclareHistogram1D(3, 1024*8, "Time interval: Pin - NaI(10ns/ch)"); // 953
  DeclareHistogram2D(4, energyBins, 1024*4, "PIN - E vs. NaI - E"); // 954
  
}

//std::ofstream outfile;


struct NaiEv
{
  int channel; 
  double energy;
  double time;
  NaiEv()
  {
    channel = -1; energy = -1; time = -1;
  }
  void Clear()
  {
    channel = -1; energy = -1; time = -1;
  }
  ~NaiEv(){};
};

bool PINProcessor::Process(RawEvent &event)
{
  if (!EventProcessor::Process(event))
    return false;

  static const vector<ChanEvent*> &pinEvents = sumMap["pin"]->GetList();

  //---------------------- HANDLING NAI SIGNALS IN THE SAME EVENT ---------------------------------------------//
  vector<ChanEvent*> NaIEvents = event.GetSummary("nai:nai", true)->GetList(); // by Yongchi Xiao; 09/25/2015
  vector<NaiEv> vecNaiEv;
  for(vector<ChanEvent*>::const_iterator it = NaIEvents.begin(); it != NaIEvents.end(); it++)
    {
      ChanEvent* chan = *it;
      string subtype   = chan->GetChanID().GetSubtype();
      int number 	= chan->GetChanID().GetLocation();
      double calEnergy = chan->GetCalEnergy();
      double naiTime   = chan->GetTime();

      NaiEv naiev;
      naiev.channel = number;
      naiev.energy = calEnergy; 
      naiev.time = naiTime; 

      vecNaiEv.push_back(naiev);
    } // by Yongchi Xiao; 09/25/2015
  //------------------------------------------------------------------------------------------------------------//


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
      if(number==2){ // GOOD FOR USE??
	plot(D_ENE+1,calEnergy);
	//----------- do NaI-PIN correlation ------------//
	for(int i = 0; i < vecNaiEv.size(); i++)
	  {
	    if(vecNaiEv.at(i).channel < 4){
	      plot(3, vecNaiEv.at(i).time - pinTime + 2000); // 953, offset = 2000
	      plot(4, calEnergy, vecNaiEv.at(i).energy); // 954
	    }
	  }
	//----------------------------------------------//
	
      }
  }
      
      
      EndProcess();
      return true;
}


