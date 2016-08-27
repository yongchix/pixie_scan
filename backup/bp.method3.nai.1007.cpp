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

using namespace dammIds::nai;
using namespace std;

namespace dammIds { // offset = 930; range = 20
  namespace nai {	
    const int DD_POSENE = 0;
    const int D_SUM1    = 1;
    const int D_SUM2    = 2;
    const int D_SUMALL  = 3;
    const int D_ENE     = 4;
    const int DD_ENE_V_T = 15;
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

  DeclareHistogram2D(DD_POSENE, posBins, energyBins, "NaI Position and Energy"); // 930
  
  DeclareHistogram1D(D_SUM1,energyBins2,"Nai sum 0-3");//931
  DeclareHistogram1D(D_SUM2,energyBins2,"Nai sum 4-7");// 932
  DeclareHistogram1D(D_SUMALL,energyBins2,"Nai sum 0-8");// 933
  
  DeclareHistogram1D(D_ENE,energyBins2,"NaI ch 0"); //934
  DeclareHistogram1D(D_ENE+1,energyBins2,"NaI ch 1");
  DeclareHistogram1D(D_ENE+2,energyBins2,"NaI ch 2");
  DeclareHistogram1D(D_ENE+3,energyBins2,"NaI ch 3");
  DeclareHistogram1D(D_ENE+4,energyBins2,"NaI ch 4");
  DeclareHistogram1D(D_ENE+5,energyBins2,"NaI ch 5");
  DeclareHistogram1D(D_ENE+6,energyBins2,"NaI ch 6");
  DeclareHistogram1D(D_ENE+7,energyBins2,"NaI ch 7");// 941

  DeclareHistogram2D(DD_ENE_V_T,energyBins2,energyBins3,"NaI Energy vs. time");
 
}

YXEvent::CorrEv naiPair; // added by Yongchi Xiao; 10/06/2015

static double zongshu = 0; 
static double teshu = 0;

bool NaIProcessor::Process(RawEvent &event)
{
  if (!EventProcessor::Process(event))
    return false;

  zongshu++;
  //    cout << "\n >> Now processing a new event \n";

  //static const vector<ChanEvent*> &naiEvents = sumMap["nai"]->GetList();
  vector<ChanEvent*> naiEvents = event.GetSummary("nai:nai", true)->GetList(); // test
  vector<ChanEvent*> pinEvents = event.GetSummary("pin:pin", true)->GetList(); // by Yongchi Xiao; 10/06/2015
  
  double sum1=0;
  double sum2=0;
  double sumall=0;

  data.Clear();

  vector<YXEvent::SimpleEv> vecNaiEventsAll, vecPinEventsAll;
  bool hasPin = false;

  // handle PIN events
  for (vector<ChanEvent*>::const_iterator it = pinEvents.begin(); it != pinEvents.end(); it++) 
    {
      ChanEvent *chan = *it;
      
      string subtype   = chan->GetChanID().GetSubtype();
      int number        = chan->GetChanID().GetLocation();
      double calEnergy = chan->GetCalEnergy();
      double pinTime   = chan->GetTime();
      if(number==0 || number ==2) // channel switch
	{
	  YXEvent::SimpleEv pinev;
	  pinev.AssignValue(pinTime, calEnergy, number, subtype);
	  vecPinEventsAll.push_back(pinev);
	  //	  	  cout << setprecision(15) << pinTime << endl;
	}
    }if(vecPinEventsAll.size() > 0) hasPin = true;


  // handle NaI events
  for (vector<ChanEvent*>::const_iterator it = naiEvents.begin(); it != naiEvents.end(); it++) 
    {
      ChanEvent *chan = *it;
      
      string subtype   = chan->GetChanID().GetSubtype();
      int number 	= chan->GetChanID().GetLocation();
      double calEnergy = chan->GetCalEnergy();
      double naiTime   = chan->GetTime();

      YXEvent::SimpleEv naiev; 
      naiev.AssignValue(naiTime, calEnergy, number, subtype);
      vecNaiEventsAll.push_back(naiev);

      if(number<=8)
	{
	  plot(DD_POSENE,number,calEnergy);
	  plot(D_ENE+number,calEnergy);
	  sumall+=calEnergy;
	}
      if(number<=3)
	sum1+=calEnergy; // sum over ch 0-3
      if(number>=4 && number <=7)
	sum2+=calEnergy; // sum over ch 4-7
      if(sum1>0)
	plot(D_SUM1,sum1);
      if(sum2>0)
	plot(D_SUM2,sum2);
      if(sumall>0)
	plot(D_SUMALL,sumall);
  }

  /*
  if(vecNaiEventsAll.size() > 0 && vecPinEventsAll.size() > 0)
    {
      teshu++;
      cout << "\n >> Processing..." << endl;
      for(int i = 0; i < vecPinEventsAll.size(); i++)
	{
	  YXEvent::SimpleEv pine = vecPinEventsAll.at(i);
	  cout << pine.channel << ", "
	       << setprecision(15) << pine.time << endl;
	}
      cout << "------------------" << endl;
      for(int i = 0; i < vecNaiEventsAll.size(); i++)
	{
	  YXEvent::SimpleEv naie = vecNaiEventsAll.at(i);
	  cout << naie.channel << ", "
	       << setprecision(15) << naie.time << endl;
	}
      cout << teshu/zongshu * 100 << "%" << endl;
    }
  */

  // up to this point, vecNaiEventsAll and vecPinEventsAll have been done
  // now looking for 511 pairs

  double time_p = -1;
  double energy_o = -1, energy_i = -1;
  bool foundPairP[6] = {false, false, false, false, false, false};
  
  YXEvent::SimpleEv groupNaiEventsP[8];

  for(int i = 0; i < vecNaiEventsAll.size(); i++)
    {
      int channel_ = vecNaiEventsAll.at(i).channel;
      double eCheck = vecNaiEventsAll.at(i).energy;
      
      bool naiOk_p = false;
      
      if( channel_ < 4 // check gamma energy qualification
	  && eCheck < 275.53//323//344 
	  && eCheck > 220.57//216//166
	  )
	naiOk_p = true;
      else if( channel_ > 3
	       && eCheck < 122.12//180//123.7
	       && eCheck > 29.59//60//28.6
	       )
	naiOk_p = true;
      
      if(naiOk_p)
	groupNaiEventsP[channel_] = vecNaiEventsAll.at(i); // fill a possible 511keV photon
      
      // --------------------- handling any 511-pairs -------------------------------------------//
      if( (groupNaiEventsP[0].time > 0 && groupNaiEventsP[5].time > 0) 
	  && abs(groupNaiEventsP[0].time - groupNaiEventsP[5].time) <= 5
	  ) foundPairP[0] = true; 
      
      if( (groupNaiEventsP[1].time > 0 && groupNaiEventsP[6].time > 0)
	  && abs(groupNaiEventsP[1].time - groupNaiEventsP[6].time) <= 5 
	  ) foundPairP[1] = true;
      
      if( (groupNaiEventsP[2].time > 0 && groupNaiEventsP[7].time > 0)
	  && abs(groupNaiEventsP[2].time - groupNaiEventsP[7].time) <= 5 
	  ) foundPairP[2] = true;
      
      if( (groupNaiEventsP[3].time > 0 && groupNaiEventsP[4].time > 0)
	  && abs(groupNaiEventsP[3].time - groupNaiEventsP[4].time) <= 5 
	  ) foundPairP[3] = true;
      
      if( (groupNaiEventsP[4].time > 0 && groupNaiEventsP[6].time > 0)
	  && abs(groupNaiEventsP[4].time - groupNaiEventsP[6].time) <= 5 // almost at the same time;
	  && abs( (groupNaiEventsP[4].energy - groupNaiEventsP[6].energy)/groupNaiEventsP[4].energy*100. < 3.5 ) // energy almost the same
	  ) foundPairP[4] = true;
      
      if( (groupNaiEventsP[5].time > 0 && groupNaiEventsP[7].time > 0)
	  && abs(groupNaiEventsP[5].time - groupNaiEventsP[7].time) <= 5 // almost at the same time;
	  && abs( (groupNaiEventsP[5].energy - groupNaiEventsP[7].energy)/groupNaiEventsP[5].energy*100. < 3.5 ) // energy almost the same
	  ) foundPairP[5] = true;
      //------------------------------------------------------------------------------------------//
      
      if( foundPairP[0] == true)
	{
	  if(hasPin){
	    naiPair.AssignValue(groupNaiEventsP[0], groupNaiEventsP[5]);
	    groupNaiEventsP[0].Clear(); groupNaiEventsP[5].Clear(); foundPairP[0] = false;
	  }else{
	    groupNaiEventsP[0].Clear(); groupNaiEventsP[5].Clear(); foundPairP[0] = false;
	  }
	} // 0, 5
      
      if( foundPairP[1] == true)
	{
	  if(hasPin){
	    naiPair.AssignValue(groupNaiEventsP[1], groupNaiEventsP[6]);
	    groupNaiEventsP[1].Clear(); groupNaiEventsP[6].Clear(); foundPairP[1] = false;
	  }else{
	    groupNaiEventsP[1].Clear(); groupNaiEventsP[6].Clear(); foundPairP[1] = false;
	  }
	} // 1, 6
      
      if( foundPairP[2] == true)
	{
	  if(hasPin){
	    naiPair.AssignValue(groupNaiEventsP[2], groupNaiEventsP[7]);
	    groupNaiEventsP[2].Clear(); groupNaiEventsP[7].Clear(); foundPairP[2] = false;
	  }else{
	    groupNaiEventsP[2].Clear(); groupNaiEventsP[7].Clear(); foundPairP[2] = false;
	  }
	} // 2, 7
      
      if( foundPairP[3] == true)
	{
	  if(hasPin){
	    naiPair.AssignValue(groupNaiEventsP[3], groupNaiEventsP[4]);
	    groupNaiEventsP[3].Clear(); groupNaiEventsP[4].Clear(); foundPairP[3] = false;
	  }else{
	    groupNaiEventsP[3].Clear(); groupNaiEventsP[4].Clear(); foundPairP[3] = false;
	  }
	} // 3, 4
      
      if( foundPairP[4] == false )
	{
	  if(hasPin){
	    naiPair.AssignValue(groupNaiEventsP[4], groupNaiEventsP[6]);
	    groupNaiEventsP[4].Clear(); groupNaiEventsP[6].Clear(); foundPairP[4] = false;
	  }else{
	    groupNaiEventsP[4].Clear(); groupNaiEventsP[6].Clear(); foundPairP[4] = false;
	  }
	} // 4, 6
      
      if( foundPairP[5] == false)
	{
	  if(hasPin){
	    naiPair.AssignValue(groupNaiEventsP[5], groupNaiEventsP[7]);
	    groupNaiEventsP[5].Clear(); groupNaiEventsP[7].Clear(); foundPairP[5] = false;
	  }else{
	    groupNaiEventsP[5].Clear(); groupNaiEventsP[7].Clear(); foundPairP[5] = false;
	  }
	} // 5, 7
      
      //----------------------------------- end of handling preceding 511-pairs ---------------------------------------//
    }
  
  /*
  if(naiPair.corrType.compare("null") != 0) // if 511 pair is found
    {
      teshu++;
      cout << "\n >> Got a 511-pair " 
	   << teshu/zongshu*100 << " % \n";
      naiPair.Clear(); // remember to initialize it after it is used 
    }
    
  if( ((int)zongshu) % 1000000 == 0)
    cout << "\n" << teshu << "/" << zongshu << "\n";
  */        
      
  //    cout << "Event processed! << \n";
  
      EndProcess();
      return true;
}
