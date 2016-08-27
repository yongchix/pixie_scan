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

#include "YXEvent.hpp"

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
  //  DeclareHistogram2D(DD_POSENE, posBins, energyBins, "PIN Position and Energy");
  DeclareHistogram1D(D_ENE,energyBins,"Pin1"); // by YX
  DeclareHistogram1D(D_ENE+1,energyBins,"Pin2");

  /*
  DeclareHistogram2D(3, 1000, 1000, "PIN-EF vs. PIN-EB"); // by Yongchi Xiao; 05/05/2016
  DeclareHistogram2D(4, 1000, 500, "PIN-EF vs. Dt"); // 05/05/2016
  DeclareHistogram2D(5, 1000, 1000, "PIN-EF vs. PIN-EB"); // 05/10/2016
  */

  // 05/10/2016, PIN vs. NaI
  /*
  DeclareHistogram2D(10, 400, 1000, "PIN-EB vs. NaI-plugin"); // 960  
  DeclareHistogram2D(11, 400, 1000, "PIN-EF vs. NaI-plugin"); // 961
  DeclareHistogram2D(20, 400, 1000, "PIN-EB vs. NaI-barrel"); // 970
  DeclareHistogram2D(21, 400, 1000, "PIN-EF vs. NaI-barrel"); // 971
  */
}


bool PINProcessor::Process(RawEvent &event)
{
  if (!EventProcessor::Process(event))
    return false;
  /*
  static const vector<ChanEvent*> &pinEvents = sumMap["pin"]->GetList();
  vector<ChanEvent*> naiEvents = event.GetSummary("nai:nai", true)->GetList();


  vector<SimpleEvent> vecPINf, vecPINb; // 05/05/2016
  SimpleEvent se;

  data.Clear();
  for (vector<ChanEvent*>::const_iterator it = pinEvents.begin();
       it != pinEvents.end(); it++) {
      ChanEvent *chan = *it;

      string subtype   = chan->GetChanID().GetSubtype();
      int number 	= chan->GetChanID().GetLocation();
      double calEnergy = chan->GetCalEnergy(); // should be the same with StripEvent.E; by YX
      double pinTime   = chan->GetTime(); // count at the same time with DSSD? by YX

      if(number == 0){ // back pin
	se.AssignValue(pinTime, calEnergy, number, "pin-b");
	vecPINb.push_back(se);
      }
      else if(number > 0){ // front pin
	se.AssignValue(pinTime, calEnergy, number, "pin-f");
	vecPINf.push_back(se);
      }

      using std::cout;
      using std:: endl; 
      plot(DD_POSENE,number,calEnergy);

      //      cout << number << endl;
     
      if(number == 0){
	plot(D_ENE,calEnergy); 
      }
      else if(number > 0){
	plot(D_ENE+1,calEnergy);
      }
      }
  */

  /*
  bool found511 = false;
  vector<pair<int, double> > ref;
  for(vector<ChanEvent*>::const_iterator itn = naiEvents.begin();
      itn != naiEvents.end();
      itn++){
    ChanEvent *chan = *itn;
    string subtype   = chan->GetChanID().GetSubtype();
    int number        = chan->GetChanID().GetLocation();
    double calEnergy = chan->GetCalEnergy();
    double naiTime   = chan->GetTime(); 
    switch(number){
    case 0: if(calEnergy < 295 && calEnergy > 149){ found511 = true; } break;
    case 1: if(calEnergy < 295 && calEnergy > 149){ found511 = true; } break;
    case 2: if(calEnergy < 295 && calEnergy > 149){ found511 = true; } break;
    case 3: if(calEnergy < 295 && calEnergy > 149){ found511 = true; } break;
    case 4: if(calEnergy < 146 && calEnergy > 65){ found511 = true; } break;
    case 5: if(calEnergy < 146 && calEnergy > 65){ found511 = true; } break;
    case 6: if(calEnergy < 146 && calEnergy > 65){ found511 = true; } break;
    case 7: if(calEnergy < 146 && calEnergy > 65){ found511 = true; } break;
    }
    //    if(found511) break;
    ref.push_back(make_pair(number, calEnergy)); 
  }
  */

  // look for PIN-f/b coincidence; 05/05/2016

  /*
  if( vecPINf.size() > 0 && vecPINb.size() > 0){
    for(vector<SimpleEvent>::iterator ipf = vecPINf.begin();
	ipf != vecPINf.end();
	ipf++){
      //      SimpleEvent se;
      for(vector<SimpleEvent>::iterator ipb = vecPINb.begin();
	  ipb != vecPINb.end();
	  ipb++){
	if(abs(ipb->time - ipf->time) < 6){
	//	if(1){
	  plot(3, ipf->energy, ipb->energy); // 953
	  plot(4, ipf->energy, ipb->time - ipf->time + 100); // 954
	  // plot pin vs. nai
	  // -- pin-f vs. nai
	  for(vector<pair<int, double> >::iterator irf = ref.begin(); 
	      irf != ref.end();
	      irf++){
	    if(irf->first < 4){ // plug-in 
	      plot(10, ipb->energy, irf->second); // 960
	      plot(11, ipf->energy, irf->second); // 961
	    }else{
	    plot(20, ipb->energy, irf->second); // 970
	    plot(21, ipf->energy, irf->second); // 971
	    }
	  }
	  // if gated on 511 keV gamma-ray
	  if(found511)
	    plot(5, ipf->energy, ipb->energy); // 955
	  
	}
      }
    }
  }
  */

      
  EndProcess();
  return true;
}


