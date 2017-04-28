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
	DeclareHistogram1D(D_ENE,energyBins,"Pin1"); // by YX
	DeclareHistogram1D(D_ENE+1,energyBins,"Pin2");
}


bool PINProcessor::Process(RawEvent &event)
{
	if (!EventProcessor::Process(event))
		return false;

	pinEvents = event.GetSummary("pin", true)->GetList();
	pinBackEvents = event.GetSummary("pin:pin_back", true)->GetList();
	pinFrontEvents = event.GetSummary("pin:pin_front", true)->GetList();

	//	vector<SimpleEvent> vecPINf, vecPINb; // 05/05/2016
	//	SimpleEvent se;

	data.Clear();


	for(vector<ChanEvent*>::const_iterator it = pinEvents.begin();
		it != pinEvents.end();
		it++) {
		ChanEvent *chan = *it;
		string subtype   = chan->GetChanID().GetSubtype();
		int number 	= chan->GetChanID().GetLocation();
		double calEnergy = chan->GetCalEnergy();
		double pinTime   = chan->GetTime(); 

		//		se.AssignValue(pinTime, calEnergy, number, subtype);
		//		vecPINb.push_back(se);

		if(subtype.compare("pin_front") == 0)
			plot(D_ENE, calEnergy); // 951;
		else if(subtype.compare("pin_back") == 0)
			plot(D_ENE+1, calEnergy); // 952
	}

	      
	EndProcess();
	return true;
}


