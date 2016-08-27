/** \file NaIProcessor.cpp
 * NaI class takes the signals from the mcp detector
 *   and calculates a 2D position based on the readout
 * SNL 2-2-08, Modified DTM 9-09
 */
#include "DammPlotIds.hpp"
#include "NaIProcessor.hpp"
#include "RawEvent.hpp"
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

  DeclareHistogram2D(12, 1042*8, 8, "Time diff. between NaI channels(10ns/ch)"); // 942

  DeclareHistogram2D(DD_ENE_V_T,energyBins2,energyBins3,"NaI Energy vs. time");
 
}


bool NaIProcessor::Process(RawEvent &event)
{
  double timeStamp[4] = {0, 0, 0, 0};

  if (!EventProcessor::Process(event))
    return false;

  static const vector<ChanEvent*> &naiEvents = sumMap["nai"]->GetList();
  
  double sum1=0;
  double sum2=0;
  double sumall=0;
  
  data.Clear();
  for (vector<ChanEvent*>::const_iterator it = naiEvents.begin();
       it != naiEvents.end(); it++) {
      ChanEvent *chan = *it;

      string subtype   = chan->GetChanID().GetSubtype();
      int number 	= chan->GetChanID().GetLocation();
      double calEnergy = chan->GetCalEnergy();
      double naiTime   = chan->GetTime();
      using std::cout;
      using std:: endl; 
      if(number<=8){
	plot(DD_POSENE,number,calEnergy);
	plot(D_ENE+number,calEnergy);
	sumall+=calEnergy;
      }
      if(number<=3){
	sum1+=calEnergy; // sum over ch 0-3
      }
      if(number>=4 && number <=7){
	sum2+=calEnergy; // sum over ch 4-7
      }
   
      if(sum1>0){
	plot(D_SUM1,sum1);
      }
      if(sum2>0){
	plot(D_SUM2,sum2);
      }
      if(sumall>0){
	plot(D_SUMALL,sumall);
	
      }

      if(number < 4)
	{
	  if(number == 3){
	    if(timeStamp[0] > 0)
	      plot(12, naiTime - timeStamp[0] + 2000, 3);} // ch3 - ch0
	  else 
	    if(timeStamp[number + 1] > 0)
	      plot(12, naiTime - timeStamp[number + 1] + 2000, number); // ch(i) - ch(i+1)
	  timeStamp[number] = naiTime;// by Yongchi Xiao; 10/04/2015
	}

    //  cout << (naiTime - 1.0414e11)/1e7 << endl; 
     // plot(DD_ENE_V_T,(naiTime - 1.0414e11)/1e8,calEnergy/5);

      


  }
      
        
      EndProcess();
      return true;

      

}


