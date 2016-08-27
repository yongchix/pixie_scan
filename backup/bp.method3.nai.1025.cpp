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
  //  DeclareHistogram1D(D_SUM2,energyBins2,"Nai sum 4-7");// 932
  //  DeclareHistogram1D(D_SUMALL,energyBins2,"Nai sum 0-8");// 933
  
  DeclareHistogram1D(D_ENE,energyBins2,"NaI ch 0"); //934
  DeclareHistogram1D(D_ENE+1,energyBins2,"NaI ch 1");
  DeclareHistogram1D(D_ENE+2,energyBins2,"NaI ch 2");
  DeclareHistogram1D(D_ENE+3,energyBins2,"NaI ch 3");
  DeclareHistogram1D(D_ENE+4,energyBins2,"NaI ch 4");
  DeclareHistogram1D(D_ENE+5,energyBins2,"NaI ch 5");
  DeclareHistogram1D(D_ENE+6,energyBins2,"NaI ch 6");
  DeclareHistogram1D(D_ENE+7,energyBins2,"NaI ch 7");// 941

  //  DeclareHistogram1D(14, 10, "Types of NaI-correlation"); // 944
  //  DeclareHistogram2D(15, 1024, 1024, "NaI-E(in) vs. NaI-E(out)"); // 945
 
}

static double zongshu = 0; 
static double teshu = 0;

bool NaIProcessor::PreProcess(RawEvent &event)
{
  if (!EventProcessor::PreProcess(event))
    return false;
  data.Clear();

  zongshu++;

  vector<ChanEvent*> naiEvents = event.GetSummary("nai:nai", true)->GetList(); // test
  vector<ChanEvent*> pinEvents = event.GetSummary("pin:pin", true)->GetList(); // by Yongchi Xiao; 10/06/2015
  vector<SimpleEv> vecNaiEventsAll;
  bool foundPairP[6] = {false, false, false, false, false, false};
  SimpleEv groupNaiEventsP[8];

  bool hasPin = false;
  if(pinEvents.size() > 0) hasPin = true;// 511-pairs gated on PIN1/2

  
  if(naiPair.GetCorrType().compare("nainai") == 0)
    naiPair.AddAuxNai(naiEvents.size()); // accumulating NaI hits
  


  // handle NaI events
  for (vector<ChanEvent*>::const_iterator it = naiEvents.begin(); it != naiEvents.end(); it++) 
    {
      ChanEvent *chan = *it;     
      string subtype   = chan->GetChanID().GetSubtype();
      int number 	= chan->GetChanID().GetLocation();
      double calEnergy = chan->GetCalEnergy();
      double naiTime   = chan->GetTime();

      SimpleEv naiev; 
      naiev.AssignValue(naiTime, calEnergy, number, subtype);
      vecNaiEventsAll.push_back(naiev);
    }
      
  for(int i = 0; i < vecNaiEventsAll.size(); i++)
    {
      int channel_ = vecNaiEventsAll.at(i).GetChannel();
      double eCheck = vecNaiEventsAll.at(i).GetEnergy();
      bool naiOk_p = false;
      
      if( channel_ < 4 // check gamma energy qualification
	  && eCheck < 299//275.53//323//344 
	  && eCheck > 143//220.57//216//166
	  )
	naiOk_p = true;
      else if( channel_ > 3
	       && eCheck < 633//122.12//180//123.7
	       && eCheck > 113//29.59//60//28.6
	       )
	naiOk_p = true;
      
      if(naiOk_p){
	groupNaiEventsP[channel_] = vecNaiEventsAll.at(i); // fill a possible 511keV photon
	//	cout << setprecision(15) << groupNaiEventsP[channel_].GetTime() << endl;
	naiOk_p = false;
      } // good up to this point;
      
      
      // --------------------- handling any 511-pairs -------------------------------------------//
      if( (groupNaiEventsP[0].GetTime() > 0 && groupNaiEventsP[5].GetTime() > 0) 
	  && abs(groupNaiEventsP[0].GetTime() - groupNaiEventsP[5].GetTime()) <= 5
	  && i == vecNaiEventsAll.size() - 1
	  ){ foundPairP[0] = true; break; }
      
      if( (groupNaiEventsP[1].GetTime() > 0 && groupNaiEventsP[6].GetTime() > 0)
	  && abs(groupNaiEventsP[1].GetTime() - groupNaiEventsP[6].GetTime()) <= 5 
	  && i == vecNaiEventsAll.size() - 1
	  ) {foundPairP[1] = true; break; }
      
      if( (groupNaiEventsP[2].GetTime() > 0 && groupNaiEventsP[7].GetTime() > 0)
	  && abs(groupNaiEventsP[2].GetTime() - groupNaiEventsP[7].GetTime()) <= 5 
	  && i == vecNaiEventsAll.size() - 1
	  ) {foundPairP[2] = true; break; }
      
      if( (groupNaiEventsP[3].GetTime() > 0 && groupNaiEventsP[4].GetTime() > 0)
	  && abs(groupNaiEventsP[3].GetTime() - groupNaiEventsP[4].GetTime()) <= 5 
	  && i == vecNaiEventsAll.size() - 1
	  ) {foundPairP[3] = true; break; }
      
      if( (groupNaiEventsP[4].GetTime() > 0 && groupNaiEventsP[6].GetTime() > 0)
	  && abs(groupNaiEventsP[4].GetTime() - groupNaiEventsP[6].GetTime()) <= 5 // almost at the same time;
	  && abs( (groupNaiEventsP[4].GetEnergy() - groupNaiEventsP[6].GetEnergy())/groupNaiEventsP[4].GetEnergy()*100. < 3.5 ) // energy almost the same
	  && i == vecNaiEventsAll.size() - 1
	  ) {foundPairP[4] = true; break; }
      
      if( (groupNaiEventsP[5].GetTime() > 0 && groupNaiEventsP[7].GetTime() > 0)
	  && abs(groupNaiEventsP[5].GetTime() - groupNaiEventsP[7].GetTime()) <= 5 // almost at the same time;
	  && abs( (groupNaiEventsP[5].GetEnergy() - groupNaiEventsP[7].GetEnergy())/groupNaiEventsP[5].GetEnergy()*100. < 3.5 ) // energy almost the same
	  && i == vecNaiEventsAll.size() - 1
	  ) {foundPairP[5] = true; break;}
      //------------------------------------------------------------------------------------------//
    }// still good up to this point

  if(foundPairP[0] && hasPin){
    naiPair.Clear();
    naiPair.AddEvent(groupNaiEventsP[0]); 
    naiPair.AddEvent(groupNaiEventsP[5]);
    //    plot(14, 1); // 944
  }
  if(foundPairP[1] && hasPin){
    naiPair.Clear();
    naiPair.AddEvent(groupNaiEventsP[1]); 
    naiPair.AddEvent(groupNaiEventsP[6]);
    //    plot(14, 2); // 944
  }
  if(foundPairP[2] && hasPin){
    naiPair.Clear();
    naiPair.AddEvent(groupNaiEventsP[2]); 
    naiPair.AddEvent(groupNaiEventsP[7]);
    //    plot(14, 3); // 944
  }
  if(foundPairP[3] && hasPin){
    naiPair.Clear();
    naiPair.AddEvent(groupNaiEventsP[3]); 
    naiPair.AddEvent(groupNaiEventsP[4]);
    //    plot(14, 4); // 944
  }
  if(foundPairP[4] && hasPin){
    naiPair.Clear();
    naiPair.AddEvent(groupNaiEventsP[4]); 
    naiPair.AddEvent(groupNaiEventsP[6]);
    //    plot(14, 5); // 944
  }
  if(foundPairP[5] && hasPin){
    naiPair.Clear();
    naiPair.AddEvent(groupNaiEventsP[5]); 
    naiPair.AddEvent(groupNaiEventsP[7]);
    //    plot(14, 6); // 944
  }

    



  // test
  /*
  for(int j = 0; j < 6; j++)
    if(foundPairP[j] == true && hasPin){
      for(int k = 0; k < 8; k++)
	cout << groupNaiEventsP[k].GetType() << ", ";
      cout << endl;
    }
  */

  // test
  /*
    if(naiPair.GetCorrType().compare("nainai") == 0)
    {
      teshu++;
      cout << "\n Got A 511 keV Pair! " 
	   << teshu/zongshu*100 << "% \n \n";
      naiPair.Clear();
    }
  */
  
  
  
  //  cout << "\n preprocessing finished << \n";  


  return true;
}




bool NaIProcessor::Process(RawEvent &event)
{
  if (!EventProcessor::Process(event))
    return false;

   //    cout << "\n >> Now processing a new event \n";

  //static const vector<ChanEvent*> &naiEvents = sumMap["nai"]->GetList();
  vector<ChanEvent*> naiEvents = event.GetSummary("nai:nai", true)->GetList(); // test
  bool foundPairP[6] = {false, false, false, false, false, false};
  SimpleEv groupNaiEvents[8];
    
  double sum1=0;
  double sum2=0;
  double sumall=0;

  data.Clear();

  
  // handle NaI events
  for (vector<ChanEvent*>::const_iterator it = naiEvents.begin(); it != naiEvents.end(); it++) 
    {
      ChanEvent *chan = *it;
      
      string subtype   = chan->GetChanID().GetSubtype();
      int number 	= chan->GetChanID().GetLocation();
      double calEnergy = chan->GetCalEnergy();
      double naiTime   = chan->GetTime();

      if(number<=8)
	{
	  //	  plot(DD_POSENE,number,calEnergy); // 930
	  plot(D_ENE+number,calEnergy);
	  sumall+=calEnergy;
	}
      if(number<=3)
	sum1+=calEnergy; // sum over ch 0-3
      if(number>=4 && number <=7)
	sum2+=calEnergy; // sum over ch 4-7
      if(sum1>0)
	plot(D_SUM1,sum1);
      //      if(sum2>0)
	//	plot(D_SUM2,sum2); //933
      //      if(sumall>0)
	//	plot(D_SUMALL,sumall); // 932

      /*
      int ch1 = 4, ch2 = 6;
      groupNaiEvents[number].AssignValue(naiTime, calEnergy, number, subtype);
      if(abs(groupNaiEvents[ch1].GetTime() - groupNaiEvents[ch2].GetTime()) < 5
	 && groupNaiEvents[ch1].GetTime() > 0
	 )
	{
	  plot(15, groupNaiEvents[ch1].GetEnergy(), groupNaiEvents[ch2].GetEnergy()); // 945
	  groupNaiEvents[ch1].Clear();
	  groupNaiEvents[ch2].Clear();
	}
      */
             
    }
  
      
  //    cout << "Event processed! << \n";
  
  
  EndProcess();
  return true;
}
