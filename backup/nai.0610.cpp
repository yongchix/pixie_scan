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

CorrFlag naiPair; // 05/20/2016


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
  
	DeclareHistogram1D(D_SUM1, energyBins2, "Nai sum 0-3, averaged");//931
	DeclareHistogram1D(D_SUM2, energyBins2, "Nai sum 0-3");// 932
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
  
  
	//  DeclareHistogram1D(14, 1000, "tdiff. 1 vs. 6"); // 944
	//  DeclareHistogram2D(15, 1000, 1000, "NaI-1E vs. tdiff. "); // 945
	//  DeclareHistogram2D(16, 1000, 1000, "NaI-6E vs. tdiff. "); // 946  
	//  DeclareHistogram2D(17, 1000, 1000, "NaI-1E vs. NaI-6E"); // 947
  
}

bool NaIProcessor::PreProcess(RawEvent &event)
{
	if (!EventProcessor::PreProcess(event))
		return false;
	data.Clear();
	
	// clear correlation flag
	if( event.GetSummary("nai", true)->GetMult() > 0) {
		naiPair.Clear(); 
	}

	vector<ChanEvent*> naiEvents = event.GetSummary("nai:nai", true)->GetList();
	// my containers
	vector<SimpleEvent> vecNaiEventsAll;
	vector<SimpleEvent> vecNaiCh4, vecNaiCh6, vecNaiCh5, vecNaiCh7;
	// references
	double plugTime = -1;
	double plugEnergySum = -1;
	int firedCh[4] = {0, 0, 0, 0};
	int numFiredCh = 4;
	// flags
	bool foundPair = false;

	// handle NaI events
	for (vector<ChanEvent*>::const_iterator it = naiEvents.begin(); it != naiEvents.end(); it++) 
		{
			ChanEvent *chan = *it;     
			string subtype   = chan->GetChanID().GetSubtype();
			int number 	= chan->GetChanID().GetLocation();
			double calEnergy = chan->GetCalEnergy();
			double naiTime   = chan->GetTime();
			// ---
			SimpleEvent naiev;
			naiev.AssignValue(naiTime, calEnergy, number, subtype);
			if(number < 4){ // plug signals
				plugEnergySum += calEnergy;
				plugTime = naiTime;
				firedCh[number]++;
			}else{ 
				if(calEnergy < 332 && calEnergy > 200){
					switch(number){
					case 4: vecNaiCh4.push_back(naiev); break;
					case 6: vecNaiCh6.push_back(naiev); break;
					case 5: vecNaiCh5.push_back(naiev); break;
					case 7: vecNaiCh7.push_back(naiev); break;
					}
				}
			}
			vecNaiEventsAll.push_back(naiev);
		}
	for(int i = 0; i < 4; i++) {
		if(firedCh[i] == 0) numFiredCh--;
	}
	plugEnergySum /= numFiredCh;
	plot(D_SUM1, plugEnergySum); // 931
	// look for paired 511 keV gamma-rays
	// 1st step, plug vs. barrel
	if( abs(plugEnergySum - 236) < 34 ) {
		for(vector<SimpleEvent>::iterator it = vecNaiEventsAll.begin();
			it != vecNaiEventsAll.end();
			it++){
			if( it->channel > 3 
				&& it->energy > 200 && it->energy < 332
				&& it->time - plugTime > -24 && it->time - plugTime < 7
				){
				foundPair = true;
				naiPair.Mark(it->time, foundPair); 
				break;
			}
		}
	}
	// 2nd step, ch 4 vs. ch 6
	if(foundPair == false){
		for(vector<SimpleEvent>::iterator it4 = vecNaiCh4.begin();
			it4 != vecNaiCh4.end();
			it4++){
			for(vector<SimpleEvent>::iterator it6 = vecNaiCh6.begin();
				it6 != vecNaiCh6.end();
				it6++){
				if(abs(it4->time - it6->time) < 13) {
					foundPair = true; 
					naiPair.Mark(it4->time, foundPair);
					break;
				}
			}if(foundPair == true) break;
		}
	}
	// 3rd step, ch 5 vs. ch 7
	if(foundPair == false){
		for(vector<SimpleEvent>::iterator it5 = vecNaiCh5.begin();
			it5 != vecNaiCh5.end();
			it5++){
			for(vector<SimpleEvent>::iterator it7 = vecNaiCh7.begin();
				it7 != vecNaiCh7.end();
				it7++){
				if(abs(it5->time - it7->time) < 13) {
					foundPair = true; 
					naiPair.Mark(it5->time, foundPair);
					break;
				}
			}if(foundPair == true) break;
		}
	}
  
	return true;
}




bool NaIProcessor::Process(RawEvent &event)
{
	if (!EventProcessor::Process(event))
		return false;

	//static const vector<ChanEvent*> &naiEvents = sumMap["nai"]->GetList();
	vector<ChanEvent*> naiEvents = event.GetSummary("nai:nai", true)->GetList(); // test    
	double sum1=0;

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
			
			if(number<=3)
				sum1+=calEnergy; // sum over ch 0-3
			//			if(number>=4 && number <=7)
			//				sum2+=calEnergy; // sum over ch 4-7
			//			if(sum1>0)
			//				plot(D_SUM1,sum1); // 931
			if(sum1>0)
				plot(D_SUM2, sum1); // 932
		}
      
  
	EndProcess();
	return true;
}
