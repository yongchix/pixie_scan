/** \file Dssd4SHEProcessor.hpp
 *
 */

#ifndef __DSSD4SHE_PROCESSOR_HPP_
#define __DSSD4SHE_PROCESSOR_HPP_

#include <vector>
#include <utility>
#include "EventProcessor.hpp"
#include "RawEvent.hpp"
#include "SheCorrelator.hpp"
#include "Trace.hpp"

namespace dammIds { 
    namespace dssd4she {
        const int D_ENERGY_X = 0;
        const int D_ENERGY_Y = 1;

        const int D_DTIME = 2;
        const int D_MWPC_MULTI = 3;
        const int D_ENERGY_CORRELATED_SIDE = 4;
        const int D_DTIME_SIDE = 5;

        const int DD_ENERGY_DT__DSSD_MWPC = 6;
        const int DD_DE_E__DSSD_VETO = 7;
	
	const int D_BACKE_W_FRONT_SAT= 8;
	const int D_FRONTE_W_BACK_SAT=9;
	
	const int D_ENERGY_IMPLANT = 10;
        const int D_ENERGY_DECAY = 11;
        const int D_ENERGY_LIGHT = 12;
        const int D_ENERGY_UNKNOWN = 13;
        const int D_ENERGY_FISSION = 14;
        const int D_ENERGY_DECAY_BEAMSTOP = 15;

        const int D_ENERGY_WITH_VETO = 16;
        const int D_ENERGY_WITH_MWPC = 17;
        const int D_ENERGY_WITH_VETO_MWPC = 18;
        const int D_ENERGY_NO_VETO_MWPC = 19;

        const int DD_EVENT_POSITION = 21;
        const int DD_EVENT_POSITION_FROM_E = 22;
	const int DD_LIGHT_POSITION = 23;        
        const int DD_UNKNOWN_POSITION = 24;
        const int DD_FISSION_POSITION = 25;
	const int DD_IMPLANT_POSITION = 26;
        const int DD_DECAY_POSITION = 27;

        const int DD_ENERGY_POS_X_TRACE = 29; 
        const int DD_ENERGY_POS_Y_TRACE = 30;

        const int DD_EVENT_ENERGY__X_POSITION_IMP = 31;
        const int DD_EVENT_ENERGY__Y_POSITION_IMP = 32;
        const int DD_EVENT_ENERGY__X_POSITION_DEC = 33;
        const int DD_EVENT_ENERGY__Y_POSITION_DEC = 34;
        const int DD_EVENT_ENERGY_COMP_X_POSITION_IMP = 35;
        const int DD_EVENT_ENERGY_COMP_Y_POSITION_IMP = 36;
        const int DD_EVENT_ENERGY_COMP_X_POSITION_DEC = 37;
        const int DD_EVENT_ENERGY_COMP_Y_POSITION_DEC = 38;  
	//39
        const int DD_FRONTE__BACKE = 40;
        const int DD_MAXEVENT_ENERGY__X_POSITION_IMP = 41;
        const int DD_MAXEVENT_ENERGY__Y_POSITION_IMP = 42;
        const int DD_MAXEVENT_ENERGY__X_POSITION_DEC = 43;
        const int DD_MAXEVENT_ENERGY__Y_POSITION_DEC = 44;
    	const int DD_MAXEVENT_ENERGY_COMP_X_POSITION_IMP = 45;
        const int DD_MAXEVENT_ENERGY_COMP_Y_POSITION_IMP = 46;
	    const int DD_MAXEVENT_ENERGY_COMP_X_POSITION_DEC = 47;
        const int DD_MAXEVENT_ENERGY_COMP_Y_POSITION_DEC = 48;
    	//49
	    const int DD_ENERGY_DECAY_TIME_GRANX = 50;
        
	/** Diagnostic **/
        const int DD_ENERGY__POSX_T_MISSING = 60;
        const int DD_ENERGY__POSY_T_MISSING = 61; 
        const int DD_DENERGY__DPOS_X_CORRELATED = 62; 
        const int DD_DENERGY__DPOS_Y_CORRELATED = 63;

	/** Plotting Pixel Correlated events **/
	const int DD_CHAINS_ENERGY_V_TIME = 70;	
	const int DD_CHAIN_NUM_ALPHA = 71;
	const int DD_TOF_A_EVENT = 72;
	const int DD_MWPC_ENERGY_A_EVENT = 73;
	const int DD_CHAIN_ALPHA_V_ALPHA = 74;

	const int DD_CHAIN_NUM_FISSION = 75;
	const int DD_TOF_SF_EVENT = 76;
	const int DD_MWPC_ENERGY_SF_EVENT =77;

        const int D_TOF_HEAVY = 80;
	const int D_TOF_LIGHT = 81;
	const int D_TOF_UNK = 82;
        const int D_MWPC_ENERGY_HEAVY = 83;
        const int D_MWPC_ENERGY_LIGHT = 84;
        const int D_MWPC_ENERGY_UNK = 85;
	/** Plotting Correlated Traces **/

	const int DD_TRACES_DECAY_CORR = 90;
	
    }
}


class Dssd4SHEProcessor : public EventProcessor {
public:
    Dssd4SHEProcessor(double frontBackTimeWindow, 
                      double deltaEnergy,
		      double recoilEnergyCut, 
                      double highEnergyCut, 
                      double lowEnergyCut, 
                      double fisisonEnergyCut, 
                      int numFrontStrips, 
                      int numBackStrips);
    virtual void DeclarePlots();
    virtual bool PreProcess(RawEvent &event);
    virtual bool Process(RawEvent &event);

protected:
    bool pickEventType(SheEvent& event);

    struct StripEvent {
        StripEvent() {
            Trace emptyTrace;
            t = 0;
            E = 0;
            pos = -1;
            sat = false;
            pileup = false;
	    tr = emptyTrace;
        }

        StripEvent(double energy, double time, int position,
                   bool saturated, Trace trace) {
            E = energy;
            t = time;
            pos = position;
            sat = saturated;
            pileup = false;
	    tr = trace;
        }

        double t;
        double E;
        int pos;
        bool sat;
        bool pileup;
	Trace tr;
    };

    SheCorrelator correlator_;
    /** Events matched based on energy (MaxEvent) **/
    std::vector<std::pair<StripEvent, StripEvent> > xyEventsEMatch_; 

    /** Events matched based on timing correlation  **/
    std::vector<std::pair<StripEvent, StripEvent> > xyEventsTMatch_; 

    /**Limit in seconds for the time difference between front and
     * back to be correlated. Also to find Si Side detectors correlated
     * events (escapes)*/
    double timeWindow_;

    /**Limit in keV of difference between front and back events to
     * be considered a good event */
    double deltaEnergy_;

    /**Limit in keV of difference between front and back events to
     * be considered a good event */
    double recoilEnergyCut_;

    /** Energy cut to differentiate high-energy events (fission?)
     * from implantation and alpha decays (in keV) **/
    double highEnergyCut_;

    /** Low Energy cut for interesting alphas (in keV) **/
    double lowEnergyCut_;

    /** Fission Energy cut (in keV) **/
    double fissionEnergyCut_;

};

#endif 
