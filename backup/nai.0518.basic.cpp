/** \file NaIProcessor.cpp
 * NaI class takes the signals from the mcp detector
 *   and calculates a 2D position based on the readout
 * SNL 2-2-08, Modified DTM 9-09
 */
#include "DammPlotIds.hpp"
#include "NaIProcessor.hpp"
#include "RawEvent.hpp"
#include "YXEvent.hpp" // by Yongchi Xiao; 10/07/2015
#include<iostream>
#include<fstream>
#include<iomanip>


using std::string;
using std::vector;
using std::cout;
using std::endl;

CorrEv naiPair; // by Yongchi Xiao; 10/06/2015

using namespace dammIds::nai;
using namespace std;

namespace dammIds { // offset = 930; range = 20
  namespace nai {	
    const int DD_POSENE = 0;
    const int D_SUM1    = 1;
    const int D_SUM2    = 2;
    const int D_SUMALL  = 3;
    const int D_ENE     = 4;
    
  }
}
void NaIProcessor::NaIData::Clear(void)
{
}

NaIProcessor::NaIProcessor(void) : EventProcessor(OFFSET, RANGE, "nai")
{
  associatedTypes.insert("nai");
}

void NaIProcessor::DeclarePlots(void)
{  
  const int posBins      = S4; //    16
  const int energyBins   = SE; // 16384
  const int energyBins2  = SD; // 8192
  const int energyBins3  = SB; // 2048

  //  DeclareHistogram2D(DD_POSENE, posBins, energyBins, "NaI Position and Energy"); // 930
  
  DeclareHistogram1D(D_SUM1,energyBins2,"Nai sum 0-3");//931
  DeclareHistogram1D(D_SUM2,energyBins2,"Nai sum 4-7");// 932
  //  DeclareHistogram1D(D_SUMALL,energyBins2,"Nai sum 0-8");// 933

  /*  
  DeclareHistogram1D(D_ENE,energyBins2,"NaI ch 0"); //934
  DeclareHistogram1D(D_ENE+1,energyBins2,"NaI ch 1");
  DeclareHistogram1D(D_ENE+2,energyBins2,"NaI ch 2");
  DeclareHistogram1D(D_ENE+3,energyBins2,"NaI ch 3");
  DeclareHistogram1D(D_ENE+4,energyBins2,"NaI ch 4");
  DeclareHistogram1D(D_ENE+5,energyBins2,"NaI ch 5");
  DeclareHistogram1D(D_ENE+6,energyBins2,"NaI ch 6");
  DeclareHistogram1D(D_ENE+7,energyBins2,"NaI ch 7");// 941
  */

  //  DeclareHistogram1D(14, 10, "Types of NaI-correlation"); // 944
  //  DeclareHistogram2D(15, 1024, 1024, "NaI-E(in) vs. NaI-E(out)"); // 945

  /*
  DeclareHistogram1D(14, 1000, "tdiff. 1 vs. 6"); // 944
  DeclareHistogram2D(15, 1000, 1000, "NaI-1E vs. tdiff. "); // 944  
  DeclareHistogram2D(16, 1000, 1000, "NaI-6E vs. tdiff. "); // 944  
  DeclareHistogram2D(17, 1000, 1000, "NaI-1E vs. NaI-6E"); // 947
  */
}

bool NaIProcessor::PreProcess(RawEvent &event)
{
  if (!EventProcessor::PreProcess(event))
    return false;
  data.Clear();

  vector<ChanEvent*> naiEvents = event.GetSummary("nai:nai", true)->GetList(); // test
  //  vector<ChanEvent*> pinEvents = event.GetSummary("pin:pin", true)->GetList(); // by Yongchi Xiao; 10/06/2015
  vector<SimpleEvent> vecNaiEventsAll;

  if(naiPair.GetCorrType().compare("nainai") == 0) // by Yongchi Xiao; used to be "nainaipin"
    naiPair.AddAuxNai(naiEvents.size()); // accumulating NaI hits

  // handle NaI events
  for (vector<ChanEvent*>::const_iterator it = naiEvents.begin(); it != naiEvents.end(); it++) 
    {
      ChanEvent *chan = *it;     
      string subtype   = chan->GetChanID().GetSubtype();
      int number 	= chan->GetChanID().GetLocation();
      double calEnergy = chan->GetCalEnergy();
      double naiTime   = chan->GetTime();

      SimpleEvent naiev; 
      naiev.AssignValue(naiTime, calEnergy, number, subtype);
      vecNaiEventsAll.push_back(naiev);
    }
      
  // catch qualified BARREL signals, 1st

  // catch the second signal(BARREL/PLUG-IN), 2nd


  return true;
}




bool NaIProcessor::Process(RawEvent &event)
{
  if (!EventProcessor::Process(event))
    return false;

   //    cout << "\n >> Now processing a new event \n";

  //static const vector<ChanEvent*> &naiEvents = sumMap["nai"]->GetList();
  vector<ChanEvent*> naiEvents = event.GetSummary("nai:nai", true)->GetList(); // test
    
  double sum1=0;
  double sum2=0;
  double sumall=0;

  data.Clear();

  
  // handle NaI events
  vector<SimpleEvent> vecNai1, vecNai6;
  SimpleEvent se;
  for (vector<ChanEvent*>::const_iterator it = naiEvents.begin(); it != naiEvents.end(); it++) 
    {
      ChanEvent *chan = *it;
      
      string subtype   = chan->GetChanID().GetSubtype();
      int number 	= chan->GetChanID().GetLocation();
      double calEnergy = chan->GetCalEnergy();
      double naiTime   = chan->GetTime();
      /*
      if(number<=8)
	{
	  //	  plot(DD_POSENE,number,calEnergy); // 930
	  plot(D_ENE+number,calEnergy);
	  sumall+=calEnergy;
	}
      */
      if(number<=3)
	sum1+=calEnergy; // sum over ch 0-3
      if(number>=4 && number <=7)
	sum2+=calEnergy; // sum over ch 4-7
      if(sum1>0)
	plot(D_SUM1,sum1); // 931
      if(sum2>0)
	plot(D_SUM2, sum2); // 932
    }

  // 05/16/2016
  /*
  for( vector<SimpleEvent>::iterator it1 = vecNai1.begin(); 
       it1 != vecNai1.end(); it1++){
    for( vector<SimpleEvent>::iterator it6 = vecNai6.begin();
	 it6 != vecNai6.end(); it6++){
      plot(14, 200 + it1->time - it6->time); // 944
      plot(15, it1->energy, 200 + it1->time - it6->time); // 945
      plot(16, it6->energy, 200 + it1->time - it6->time); // 946
      if( (it1->time - it6->time) < 0
	  && (it1->time - it6->time > -50) 
	  )
	plot(17, it1->energy, it6->energy); // 947
    }
  }
  */
      
  //    cout << "Event processed! << \n";
  
  
  EndProcess();
  return true;
}
