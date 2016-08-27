/** \file LiquidScintProcessor.cpp
 *
 * implementation for scintillator processor
 */
#include <vector>
#include <sstream>

#include <cmath>

#include "DammPlotIds.hpp"
#include "RawEvent.hpp"
#include "LiquidScintProcessor.hpp"
#include "TimingInformation.hpp"
#include "Trace.hpp"

using namespace std;
using namespace dammIds::liquid_scint;

namespace dammIds {
    namespace liquid_scint {
        const int DD_TQDCLIQUID       = 0;
        const int DD_MAXLIQUID        = 1;
        const int DD_DISCRIM          = 2;
        const int DD_TOFLIQUID        = 3;
        const int DD_TRCLIQUID        = 4;
        const int DD_TQDCVSDISCRIM    = 5;
        const int DD_TOFVSDISCRIM     = 6;
        const int DD_NEVSDISCRIM      = 8;
        const int DD_TQDCVSLIQTOF     = 10;
        const int DD_TQDCVSENERGY     = 12;
    }
} 

LiquidScintProcessor::LiquidScintProcessor() : 
    EventProcessor(OFFSET, RANGE, "liquid_scint")
{
    associatedTypes.insert("liquid_scint");
}

void LiquidScintProcessor::DeclarePlots(void)
{
    /** WARNING
     * This part was commented in the old ScintProcessor and is 
     * copied as is.
     */

    //To handle Liquid Scintillators
    // DeclareHistogram2D(DD_TQDCLIQUID, SC, S3, "Liquid vs. Trace QDC");
    // DeclareHistogram2D(DD_MAXLIQUID, SC, S3, "Liquid vs. Maximum");
    // DeclareHistogram2D(DD_DISCRIM, SA, S3, "N-Gamma Discrimination");
    // DeclareHistogram2D(DD_TOFLIQUID, SE, S3,"Liquid vs. TOF");
    // DeclareHistogram2D(DD_TRCLIQUID, S7, S7, "LIQUID TRACES");

    // for(unsigned int i=0; i < 2; i++) { 
    // 	DeclareHistogram2D(DD_TQDCVSDISCRIM+i, SA, SE,"Trace QDC vs. NG Discrim");
    // 	DeclareHistogram2D(DD_TOFVSDISCRIM+i, SA, SA, "TOF vs. Discrim");
    // 	DeclareHistogram2D(DD_NEVSDISCRIM+i, SA, SE, "Energy vs. Discrim");
    // 	DeclareHistogram2D(DD_TQDCVSLIQTOF+i, SC, SE, "Trace QDC vs. Liquid TOF");
    // 	DeclareHistogram2D(DD_TQDCVSENERGY+i, SD, SE, "Trace QDC vs. Energy");
    // }
}

bool LiquidScintProcessor::PreProcess(RawEvent &event){
    if (!EventProcessor::PreProcess(event))
        return false;
    return true;
}

/**
 * WARNING!
 * This part was the LiquidAnalysis function in the old ScintProcessor. 
 * It looks like written for some older version of code.
 * Before using, examine it carefully!
 */
bool LiquidScintProcessor::Process(RawEvent &event)
{
    if (!EventProcessor::Process(event))
        return false;

    static const vector<ChanEvent*> &liquidEvents = 
	event.GetSummary("liquid_scint:liquid")->GetList();
    static const vector<ChanEvent*> &betaStartEvents = 
	event.GetSummary("liquid_scint:beta:start")->GetList();
    static const vector<ChanEvent*> &liquidStartEvents = 
	event.GetSummary("liquid_scint:liquid:start")->GetList();

    vector<ChanEvent*> startEvents;
    startEvents.insert(startEvents.end(), betaStartEvents.begin(),
		       betaStartEvents.end());
    startEvents.insert(startEvents.end(), liquidStartEvents.begin(),
		       liquidStartEvents.end());
    
    for(vector<ChanEvent*>::const_iterator itLiquid = liquidEvents.begin();
	itLiquid != liquidEvents.end(); itLiquid++) {
        unsigned int loc = (*itLiquid)->GetChanID().GetLocation();
        TimingData liquid((*itLiquid));

        //Graph traces for the Liquid Scintillators
        if(liquid.discrimination == 0) {
            for(Trace::const_iterator i = liquid.trace.begin(); 
            i != liquid.trace.end(); i++)
                plot(DD_TRCLIQUID, int(i-liquid.trace.begin()), 
                    counter, int(*i)-liquid.aveBaseline);
            counter++;
        }
        
        if(liquid.dataValid) {
            plot(DD_TQDCLIQUID, liquid.tqdc, loc);
            plot(DD_MAXLIQUID, liquid.maxval, loc);

            double discrimNorm = 
            liquid.discrimination/liquid.tqdc;	    
            
            double discRes = 1000;
            double discOffset = 100;
            
            TimingCal calibration =
            GetTimingCal(make_pair(loc, "liquid"));
            
            if(discrimNorm > 0)
                plot(DD_DISCRIM, discrimNorm*discRes+discOffset, loc);
            plot(DD_TQDCVSDISCRIM, discrimNorm*discRes+discOffset,
            liquid.tqdc);
            
            if((*itLiquid)->GetChanID().HasTag("start"))
                continue;
            
            for(vector<ChanEvent*>::iterator itStart = startEvents.begin(); 
            itStart != startEvents.end(); itStart++) { 
                unsigned int startLoc = (*itStart)->GetChanID().GetLocation();
                TimingData start((*itStart));
                int histLoc = loc + startLoc;
                const int resMult = 2;
                const int resOffset = 2000;
                
                if(start.dataValid) {
                    double tofOffset;
                    if(startLoc == 0)
                    tofOffset = calibration.tofOffset0;
                    else
                    tofOffset = calibration.tofOffset1;
                    
                    double TOF = liquid.highResTime - 
                    start.highResTime - tofOffset; //in ns
                    double nEnergy = CalcEnergy(TOF, calibration.r0);
                    
                    plot(DD_TOFLIQUID, TOF*resMult+resOffset, histLoc);
                    plot(DD_TOFVSDISCRIM+histLoc, 
                    discrimNorm*discRes+discOffset, TOF*resMult+resOffset);
                    plot(DD_NEVSDISCRIM+histLoc, discrimNorm*discRes+discOffset, nEnergy);
                    plot(DD_TQDCVSLIQTOF+histLoc, TOF*resMult+resOffset, 
                    liquid.tqdc);
                    plot(DD_TQDCVSENERGY+histLoc, nEnergy, liquid.tqdc);
                }
            } //Loop over starts
        } // Good Liquid Check
    }//end loop over liquid events
    EndProcess();
    return true;
}
