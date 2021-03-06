/*! \file Correlator.cpp
 *
 *  The correlator class keeps track of where and when implantation and decay
 *  events have occurred, and then correlates each decay with its corresponding
 *  implant. A decay will only be validate if it occurred close enough in time
 *  to the implant, and the implant was well separated in time with regard to
 *  all other implants at the same location
 *  
 *  This file is derived from previous "correlator.cpp"
 *  
 *  David Miller, April 2010
 */

#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <cmath>
#include <ctime>

#include "DammPlotIds.hpp"
#include "Globals.hpp"
#include "DetectorDriver.hpp"
#include "LogicProcessor.hpp"
#include "RawEvent.hpp"
#include "Correlator.hpp"
#include "param.h"
using namespace std;
using namespace dammIds::correlator;


ListData::ListData(double t, double e, LogicProcessor *lp) : time(t), energy(e)
{
  
  logicBits[3] = '\0';

    if (lp) {       
	if (!lp->LogicStatus(3)) {
	    logicBits[0]='0'; 
	    offTime = lp->TimeOff(3, time);
	} else logicBits[0]='1'; 
	if (!lp->LogicStatus(4)) {
	    logicBits[1]='0';
	    if (logicBits[0] == '1')
	 	offTime = lp->TimeOff(4, time) + 300e-6 / pixie::clockInSeconds;
	} else logicBits[1]='1';
	if (!lp->LogicStatus(5)) {
	    logicBits[2]='0';
	    if (logicBits[0] == '1' && logicBits[1] == '1') 
		offTime = lp->TimeOff(5, time) + 600e-6 / pixie::clockInSeconds;
	} else logicBits[2]='1';

	clockCount = lp->StartCount(2);
    } else {
	logicBits[0] = logicBits[1] = logicBits[2] = 'X'; // NO LOGIC
	offTime = 0.;
	clockCount = 0;
    }
 
}

namespace dammIds {
    namespace correlator {
        const int D_CONDITION            = 0;
        const int D_TIME_BW_IMPLANTS     = 1;
        const int D_TIME_BW_ALL_IMPLANTS = 2;
    }
} // correlator namespace

// all in seconds
const double Correlator::minImpTime = 5e-3;
const double Correlator::corrTime   = 60; // used to be 3300
const double Correlator::fastTime   = 40e-6;

Correlator::Correlator() : histo(OFFSET, RANGE, "correlator"), 
			   lastImplant(NULL), lastDecay(NULL), condition(UNKNOWN_CONDITION)
{
 
}

EventInfo::EventInfo()
{
    flagged = false;
    pileUp  = false;
    hasVeto = false;
    hasTof  = false;
    beamOn  = false;

    foilTime = offTime = energy = time = position = NAN;
    boxMult = mcpMult = impMult = 0;

    type = UNKNOWN_EVENT;
    logicBits[0] = 'X'; logicBits[1] = '\0'; 
    logicBits[dammIds::logic::MAX_LOGIC] = '\0';
    generation = 0;
}

CorrelationList::CorrelationList() : std::vector<EventInfo>()
{
    flagged = false;
}

double CorrelationList::GetDecayTime() const
{
    if (empty() || back().type == EventInfo::IMPLANT_EVENT) {
	return NAN;
    } else {
	return back().dtime;
    }
}

double CorrelationList::GetImplantTime() const
{
    if (empty() || front().type != EventInfo::IMPLANT_EVENT) {
	return NAN;
    } else {
	return front().time;
    }
}
/*
void CorrelationList::Flag() 
{
    if (!empty())
	back().flagged = true;
    flagged = true;
}
*/
bool CorrelationList::IsFlagged() const
{
    return true;
}

void CorrelationList::clear()
{
    flagged = false;
    vector<EventInfo>::clear();
}

void CorrelationList::PrintDecayList() const
{

}



void Correlator::Init(RawEvent& rawev)
{
    // do nothing
}

Correlator::~Correlator()
{
    // dump any flagged decay lists which have not been output
    for (unsigned int i=0; i < arraySize; i++) {
	for (unsigned int j=0; j < arraySize; j++) {
	    if (IsFlagged(i,j))
		PrintDecayList(i,j);
	}
    }    
}

void Correlator::DeclarePlots()
{
    using namespace dammIds::correlator;
    static bool done = false;

    if (done)  {
	return;
    }

    DeclareHistogram1D(D_CONDITION, S9, "Correlator condition");
    DeclareHistogram1D(D_TIME_BW_IMPLANTS, S9, "time between implants, 100 ms/bin"); 
    DeclareHistogram1D(D_TIME_BW_ALL_IMPLANTS, SA, "time between all implants, 1 us/bin"); 

    done = true;
}


void Correlator::Correlate(EventInfo &event, 
			   unsigned int fch, unsigned int bch)
{
    using namespace dammIds::correlator;
    /*
    if (fch < 0 || fch >= arraySize || bch < 0  || bch >= arraySize) {
	plot(D_CONDITION, INVALID_LOCATION);
	return;
    }

    CorrelationList &theList = decaylist[fch][bch];
    
    double lastTime = NAN;
    double clockInSeconds = Globals::get()->clockInSeconds();

    switch (event.type) {
	case EventInfo::IMPLANT_EVENT:
	    if (theList.IsFlagged()) {
		PrintDecayList(fch, bch);
	    }	
	    lastTime = GetImplantTime(fch, bch);
	    theList.clear();
	    
	    condition = VALID_IMPLANT;
	    if ( lastImplant != NULL ) {
		double dt = event.time - lastImplant->time;
		plot(D_TIME_BW_ALL_IMPLANTS, dt * clockInSeconds / 1e-6);
	    } 
	    if ( !isnan(lastTime) ) {
		condition = BACK_TO_BACK_IMPLANT;
		event.dtime = event.time - lastTime;
		plot(D_TIME_BW_IMPLANTS, event.dtime * clockInSeconds / 100e-3);
	    } else {
		event.dtime = INFINITY;
	    }
	    event.generation = 0;
	    theList.push_back(event);
	    lastImplant = &theList.back();

	    break;
	default:	    
 	    if ( theList.empty() ) {
		break;
	    }
	    if ( isnan(theList.GetImplantTime()) ) {
		cout << "No implant time for decay list" << endl;
		break;
	    }
	    if (event.type == EventInfo::UNKNOWN_EVENT) {
		condition = UNKNOWN_CONDITION;
	    } else {
		condition = VALID_DECAY;
	    }
	    condition = VALID_DECAY; // tmp -- DTM
	    lastTime = theList.back().time;

	    double dt = event.time - theList.GetImplantTime();
	    if (dt < 0) {
		if ( dt < -5e11 && event.time < 1e9 ) {
		    cout << "Decay following pixie clock reset, clearing decay lists!" << endl;
		    cout << "  Event time: " << event.time
			 << "\n  Implant time: " << theList.GetImplantTime()
			 << "\n  DT: " << dt << endl;
		    // PIXIE's clock has most likely been zeroed due to a file marker
		    //   no chance of doing correlations
		    for (unsigned int i=0; i < arraySize; i++) {
			for (unsigned int j=0; j < arraySize; j++) {
			    if (IsFlagged(i,j)) {
				PrintDecayList(i, j);
			    }
			    decaylist[i][j].clear();
			}
		    }
		} else if (event.type != EventInfo::GAMMA_EVENT) {
		    // since gammas are processed at a different time than everything else
		    cout << "negative correlation time, DECAY: " << event.time
			 << " IMPLANT: " << theList.GetImplantTime()
			 << " DT: " << dt << endl;
		}

		event.dtime = NAN;
		break;
	    } // negative correlation itme
	    if (theList.front().dtime * clockInSeconds >= minImpTime) {
		if (dt * clockInSeconds < corrTime) {
		    // event.dtime = event.time - lastTime; // (FOR CHAINS)
		    event.dtime = event.time - theList.front().time; // FOR LERIBSS
		    if (event.dtime * clockInSeconds < fastTime && event.dtime > 0) {
			// event.flagged = true; 
		    }
		} else {
		    // event.dtime = event.time - lastTime; // (FOR CHAINS)
		    event.dtime = event.time - theList.front().time; // FOR LERIBSS
		    condition = DECAY_TOO_LATE;
		}
	    } else {
		condition = IMPLANT_TOO_SOON;
	    }
	    if (condition == VALID_DECAY) {
		event.generation = theList.back().generation + 1;
	    }
	    theList.push_back(event);
	    if (event.energy == 0 && isnan(event.time))
		cout << " Adding zero decay event " << endl;

	    if (event.flagged)
		theList.Flag();

	    if (condition == VALID_DECAY) {
		lastDecay = &theList.back();
	    } else if (condition == DECAY_TOO_LATE) {
		theList.clear();
	    }

	    break;
    }
    
    plot(D_CONDITION, condition);
    */
}

/**
 *  This correlates an event with all positions in the setup. 
 *    This is useful for cases where the event is interesting but can not be assigned
 *    to a particular implant location as in the case of external gamma-ray detectors
 */

void Correlator::CorrelateAll(EventInfo &event)
{
    for (unsigned int fch=0; fch < arraySize; fch++) {
	for (unsigned int bch=0; bch < arraySize; bch++) {
	    if (decaylistold[fch][bch].size() == 0)
		continue;
	    if (event.time - decaylistold[fch][bch].back().time <
                10e-6 / Globals::get()->clockInSeconds()) {
		// only correlate fast events for now
		Correlate(event, fch, bch);
	    }
	}
    }
}



void Correlator::CorrelateOld(RawEvent &rawev, EEventType type, unsigned int frontCh, // EEventType type; by YX
			      unsigned int backCh, double time,double energy){
  
  using namespace dammIds::correlator;
 

  if (frontCh < 0 || frontCh >= MAX_STRIP ||
      backCh < 0  || backCh >= MAX_STRIP) {
    plot(D_CONDITION, INVALID_LOCATION);
    return;
  }
  
  ImplantData &imp = implant[frontCh][backCh];
  DecayData   &dec = decay[frontCh][backCh];
  
  if (type == IMPLANT_EVENT) {
      if (imp.flagged) {
	    PrintDecayList(frontCh, backCh);
      }
      
      decaylistold[frontCh][backCh].clear();
      // decaylistold[frontCh][backCh].push_back( ListData(time, energy, logicProc) );
	
	condition = VALID_IMPLANT;
	if (lastImplantold != NULL) {
	    double dt = time - lastImplantold->time;
	    plot(D_TIME_BW_ALL_IMPLANTS, dt * pixie::clockInSeconds / 1e-6);
	}
	if (imp.implanted) {
	    condition = BACK_TO_BACK_IMPLANT;
	    imp.dtime = time - imp.time;
	    imp.tacValue = 0;
	    imp.flagged = false;
	    plot(D_TIME_BW_IMPLANTS, imp.dtime * pixie::clockInSeconds / 100e-3 );
	} else {
	    imp.implanted = true;
	    imp.dtime = INFINITY;
	}
	lastImplantold = &imp;
	imp.time = time;
    } else if (type == DECAY_EVENT && imp.implanted) {
      condition = VALID_DECAY;
      //      decaylistold[frontCh][backCh].push_back( ListData(time, energy, logicProc) );
      
      if (time < imp.time ) {	 
	double dt = time - imp.time;
	
	if (dt > 5e11 && time < 1e9 ) {
	  // PIXIE's clock has most likely been zeroed due to file marker, no chance of doing correlations
		for (unsigned int i=0; i < MAX_STRIP; i++) {
		    for (unsigned int j=0; j < MAX_STRIP; j++) {			
			if (implant[i][j].time > time) {
			    if (implant[i][j].flagged)
				PrintDecayList(i,j);
			    implant[i][j].Clear();
			}
		    }
		}
	    }
	    cout << "negative correlation time, DECAY: " << time 
		 << " IMPLANT: " << imp.time 
		 << " DT:" << dt << endl;

      }// negative correlation time
	if ( imp.dtime >= minImpTime ) {
	    if (time - imp.time < corrTime) {
		double lastTime = NAN;
		if (decaylistold[frontCh][backCh].size() == 1) {
		    lastTime = imp.time;
		} else if (decaylistold[frontCh][backCh].size() > 1) {
		    lastTime = dec.time;
		}
		double dt = (dec.time - lastTime) * pixie::clockInSeconds;
		if (!isnan(lastTime) && dt< fastTime && dt > 0) {
		    cout << "flagging " << frontCh << " , " << backCh << " in correlator." << endl;
		    Flag(frontCh, backCh);
		}
		dec.time    = time;
		dec.dtime   = time - imp.time;
		
		lastDecayold = &dec;
	    } else {
		condition = DECAY_TOO_LATE;
	    }
	} else {
	    condition = IMPLANT_TOO_SOON;
	}
    } // is decay with implant 
    else {
	condition = UNKNOWN_CONDITION;
    }
   
  
}


	
/*
 *  This correlates an event with all positions in the setup. 
 *    This is useful for cases where the event is interesting but can not be assigned
 *    to a particular implant location as in the case of external gamma-ray detectors
 */

void Correlator::CorrelateAllX(EventInfo &event, unsigned int bch)
{
  /*
    for (unsigned int fch = 0; fch < arraySize; fch++) {
	Correlate(event, fch, bch);
    }
  */
}
  
void Correlator::CorrelateAllY(EventInfo &event, unsigned int fch)
{
  /*
    for (unsigned int bch = 0; bch < arraySize; bch++) {
	Correlate(event, fch, bch);
    }
  */
}

double Correlator::GetDecayTime(void) const
{
  if (lastDecay == NULL) {
    return NAN;
  } else {
    return lastDecay->dtime;
  }
}

 /*
double Correlator::GetDecayTime(int fch, int bch) const
{
    return decaylistold[fch][bch].GetDecayTime();
}
 */
  /*
double Correlator::GetImplantTime(void) const
{
    if (lastImplant == NULL) {
	return NAN;
    } else {
	return lastImplant->time;
    }
}

double Correlator::GetImplantTime(int fch, int bch) const
{
    return decaylistold[fch][bch].GetImplantTime();
}
*/
void Correlator::Flag(int fch, int bch) 
{
  // if (!decaylistold[fch][bch].empty())
  // decaylistold[fch][bch].Flag();
  implant[fch][bch].flagged=true;
}

bool Correlator::IsFlagged(int fch, int bch)
{
  //return decaylistold[fch][bch].IsFlagged();
  return true;
}

void Correlator::PrintDecayList(unsigned int fch, unsigned int bch) const
{
    // this channel has a lot of noise

  //    cout << "Current decay list for " << fch << " , " << bch << " : " << endl;
  //  decaylistold[fch][bch].PrintDecayList();
}
