/** \file NaIProcessor.cpp
 * NaI class takes the signals from the mcp detector
 *   and calculates a 2D position based on the readout
 * SNL 2-2-08, Modified DTM 9-09
 */
#include<iostream>
#include<fstream>
#include<iomanip>

#include "DammPlotIds.hpp"
#include "DetectorDriver.hpp"
#include "NaIProcessor.hpp"
#include "PINProcessor.hpp"
#include "RawEvent.hpp"
#include "YXEvent.hpp" // by Yongchi Xiao; 10/07/2015


using std::string;
using std::vector;
using std::cout;
using std::endl;

CorrFlag corrNaiPin; // 06/17/2016

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
	//	DeclareHistogram1D(D_SUM2, energyBins2, "Nai sum 0-3");// 932
	//	DeclareHistogram2D(3, 2000, 2000,  "gamma-gamma coincidence"); // 933
	//  DeclareHistogram1D(D_SUMALL, energyBins2, "Nai sum 0-8");// 933
	//	DeclareHistogram1D(4, 2000, "plug-gamma-E"); // 934
	DeclareHistogram1D(7, 500, "Dt. Pin-F"); // 937
	DeclareHistogram1D(8, 500, "Energy Pin-F"); // 938
	DeclareHistogram1D(9, 500, "Dt. Pin-B"); // 939
	DeclareHistogram1D(10, 500, "Energy Pin-B"); // 940
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
  
	//	DeclareHistogram2D(15, 2000, 500, "plug-E vs. tdiff. "); // 945
	//	DeclareHistogram2D(16, 2000, 500, "barrel-E vs. tdiff. "); // 946  
}

bool NaIProcessor::PreProcess(RawEvent &event)
{
	if (!EventProcessor::PreProcess(event))
		return false;
	data.Clear();
	// clear correlation flag
	/* It does not seem to be necessary to reset the start point 
	 * everytime when a photon is observed. 
	 */
	//	naiPair.Clear(); 
	//	corrNaiPin.Clear();
	
	vector<ChanEvent*> naiEvents = event.GetSummary("nai:nai", true)->GetList();
	vector<ChanEvent*> pinBackEvents, pinFrontEvents;
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
	bool hasPinFront = false;
	bool hasPinBack = false;

	if( event.GetSummary("pin:pin_back", true)->GetMult() > 0) {
		hasPinBack = true;
		pinBackEvents = event.GetSummary("pin:pin_back", true)->GetList();
		//		pinBackEvents = ( (PINProcessor*)DetectorDriver::get()->GetProcessor("PINProcessor") )->GetBackPinEvents();
	} 
	if( event.GetSummary("pin:pin_front", true)->GetMult() > 0) {
		hasPinFront = true;
		pinFrontEvents = event.GetSummary("pin:pin_front", true)->GetList();
		//		pinFrontEvents = ( (PINProcessor*)DetectorDriver::get()->GetProcessor("PINProcessor") )->GetFrontPinEvents();
	}

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
				//			if(calEnergy < 560 && calEnergy > 450){
				if(true) {
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
	// handle plug energy
	for(int i = 0; i < 4; i++) {
		if(firedCh[i] == 0) numFiredCh--;
	}
	plugEnergySum /= numFiredCh;
	plugEnergySum = 1.641*plugEnergySum + 114.920; // my own calibration
	plot(D_SUM1, plugEnergySum); // 931


	/* Correlation between NaI and PIN-front
	 */
	if( abs(plugEnergySum - 505) < 55 ) {
		/*
		//		corrNaiPin.Mark(plugTime, true, plugEnergySum, 100); 
		for (vector<ChanEvent*>::iterator it = pinFrontEvents.begin();
			 it != pinFrontEvents.end();
			 it++) {
			ChanEvent *chan = *it;
			plot(5, chan->GetCalEnergy()); // 935
		}
		for (vector<ChanEvent*>::iterator it = pinBackEvents.begin();
			 it != pinBackEvents.end();
			 it++) {
			ChanEvent *chan = *it;
			plot(6, chan->GetCalEnergy()); // 936
		}
		*/
		if(hasPinFront) {
			for(vector<ChanEvent*>::iterator it = pinFrontEvents.begin();
				it != pinFrontEvents.end();
				it++) {
				ChanEvent *chan = *it;
				if( abs(chan->GetCalEnergy() - 10) < 5) {
					corrNaiPin.Clear();
					corrNaiPin.Mark(plugTime, true, plugEnergySum, chan->GetCalEnergy());
					plot(7, 200 + chan->GetTime() - plugTime); // 937
					plot(8, chan->GetCalEnergy()); // 938
				}
			}					 
		}
		
		if(hasPinBack) {
			for(vector<ChanEvent*>::iterator it = pinBackEvents.begin();
				it != pinBackEvents.end();
				it++) {
				ChanEvent *chan = *it;
				if( abs(chan->GetCalEnergy() - 10) < 5) {
					corrNaiPin.Clear();
					corrNaiPin.Mark(plugTime, true, plugEnergySum, chan->GetCalEnergy());
					plot(9, 200 + chan->GetTime() - plugTime); // 939
					plot(10, chan->GetCalEnergy()); // 940
				}
			}
		}
		
	}
	

	/*  plug vs. barrel, barrel vs. barrel
	 */

	/*	// 1st step, plug vs. barrel
	//	if( abs(plugEnergySum - 505) < 55 ) {
	if(true) {
		for(vector<SimpleEvent>::iterator it = vecNaiEventsAll.begin();
			it != vecNaiEventsAll.end();
			it++) {
			if( it->channel > 3 // barrel
				//				&& it->energy > 450 && it->energy < 560 // qualified E. in barrel
				&& it->time - plugTime > -16 && it->time - plugTime < 34 // in coincidence
				) {
				foundPair = true;
				naiPair.Mark(it->time, foundPair, plugEnergySum, it->energy); 
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
				if( abs(it4->time - it6->time) < 200) {
					foundPair = true; 
					naiPair.Mark(it4->time, foundPair);
					plot(15, it4->energy, it4->time - it6->time + 200); // 945
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
				if(abs(it5->time - it7->time) < 200) {
					foundPair = true; 
					naiPair.Mark(it5->time, foundPair);
					plot(16, it5->energy, it5->time - it7->time + 200);// 946
					break;
				}
			}if(foundPair == true) break;
		}
	}
	*/
	
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
		}
      
  
	EndProcess();
	return true;
}
