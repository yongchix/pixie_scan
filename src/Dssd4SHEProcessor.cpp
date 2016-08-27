/*! \file DssdProcessor.cpp
 *
 * The DSSD processor handles detectors of type dssd_front and dssd_back and
 *   determines whether the events are implants or decays and informs the
 *   correlator accordingly
 */

#include <algorithm>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <signal.h>

#include "Dssd4SHEProcessor.hpp"
#include "DammPlotIds.hpp"
#include "Globals.hpp"
#include "Messenger.hpp"
#include "Notebook.hpp"
#include "RawEvent.hpp"

#include "CfdAnalyzer.hpp"
#include "DoubleTraceAnalyzer.hpp"
#include "FittingAnalyzer.hpp"
#include "TauAnalyzer.hpp"
#include "TraceAnalyzer.hpp"
#include "TraceExtracter.hpp"
#include "WaveformAnalyzer.hpp"
#include "DetectorDriver.hpp"
#include "Notebook.hpp"
using namespace dammIds::dssd4she;
using namespace std;

Dssd4SHEProcessor::Dssd4SHEProcessor(double timeWindow,
                                     double deltaEnergy,
				     double recoilEnergyCut,
                                     double highEnergyCut,
                                     double lowEnergyCut,
                                     double fissionEnergyCut,
                                     int numBackStrips,
                                     int numFrontStrips) :
    EventProcessor(OFFSET, RANGE, "dssd4she"),

    correlator_(numBackStrips, numFrontStrips)
{
    timeWindow_ = timeWindow;
    deltaEnergy_ = deltaEnergy;
    recoilEnergyCut_ = recoilEnergyCut;
    highEnergyCut_ = highEnergyCut;
    lowEnergyCut_ = lowEnergyCut;
    fissionEnergyCut_ = fissionEnergyCut;
    name = "dssd";
    associatedTypes.insert("dssd_front");
    associatedTypes.insert("dssd_back");

    stringstream ss;
    ss << fixed 
       << "#T" 
       << " " << setw(12) << "E (keV)"
       << " " << setw(12) << "t (ms)"  
       << " M" << " "
       << " B" << " "
       << " V" << " "
       << " E" << " "
       << endl;
    Notebook::get()->report(ss.str());
}


void Dssd4SHEProcessor::DeclarePlots(void)
{
    using namespace dammIds::dssd;

    const int energyBins = SE;
    const int energyBins2 = SB;
    const int xBins = S7;
    const int yBins = S6;
    const int timeBins = S8;

    DeclareHistogram1D(D_ENERGY_X, energyBins, "Energy/10 dssd X strips");
    DeclareHistogram1D(D_ENERGY_Y, energyBins, "Energy/10 dssd Y strips");
    DeclareHistogram1D(D_DTIME, S8, "Pairs time diff in 10 ns (+ 1 bin)");
    DeclareHistogram1D(D_MWPC_MULTI, S5, "MWPC multiplicity");
    DeclareHistogram1D(D_ENERGY_CORRELATED_SIDE, energyBins, 
                       "Energy Side corr. with DSSD");
    DeclareHistogram1D(D_DTIME_SIDE, S8, 
                        "Side det. time diff in 10 ns (+ 1 bin)");

    DeclareHistogram2D(DD_ENERGY_DT__DSSD_MWPC, 
		       SB, S8, "DSSD energy/100 vs DT (10 ns) to MWPC");

    DeclareHistogram2D(DD_DE_E__DSSD_VETO, 
		       SB, SB, "DSSD energy/100 vs veto/100");

    /**Declare Saturation diagnostics **/
    DeclareHistogram1D(D_BACKE_W_FRONT_SAT, energyBins, "Back Energy with Front Saturation");
    DeclareHistogram1D(D_FRONTE_W_BACK_SAT, energyBins, "Front Energy with Back Saturation");

    DeclareHistogram1D(D_ENERGY_IMPLANT,
		       energyBins2, "DSSD energy/100 implant");
    DeclareHistogram1D(D_ENERGY_DECAY,
		       energyBins2, "DSSD energy/100 decay");
    DeclareHistogram1D(D_ENERGY_LIGHT,
		       energyBins2, "DSSD energy/100 light ion");
    DeclareHistogram1D(D_ENERGY_UNKNOWN,
		       energyBins2, "DSSD energy/100 unknown");
    DeclareHistogram1D(D_ENERGY_FISSION,
		       energyBins2, "DSSD energy/100 fission");
    DeclareHistogram1D(D_ENERGY_DECAY_BEAMSTOP,
		       energyBins, "DSSD energy*1 alpha beam stopped");
    DeclareHistogram1D(D_ENERGY_WITH_VETO, energyBins, 
                      "Energy dssd/10 coin. veto");
    DeclareHistogram1D(D_ENERGY_WITH_MWPC, energyBins, 
                      "Energy dssd/10 coin. mwpc");
    DeclareHistogram1D(D_ENERGY_WITH_VETO_MWPC, energyBins, 
                      "Energy dssd/10 coin. veto and mwpc");
    DeclareHistogram1D(D_ENERGY_NO_VETO_MWPC, energyBins, 
                      "Energy dssd/10 coin. no veto and mwpc");
    DeclareHistogram1D(D_MWPC_ENERGY_WITH_DSSD, energyBins, "Energy of MWPC coin. with Dssd (real Recoil)");

    DeclareHistogram2D(DD_EVENT_POSITION, 
		       xBins, yBins, "DSSD all events positions");
    DeclareHistogram2D(DD_EVENT_POSITION_FROM_E, 
		       xBins, yBins, "DSSD position all max event");
    DeclareHistogram2D(DD_LIGHT_POSITION, 
		       xBins, yBins, "DSSD position light ion");
    DeclareHistogram2D(DD_UNKNOWN_POSITION, 
		       xBins, yBins, "DSSD position unknown");
    DeclareHistogram2D(DD_FISSION_POSITION, 
		       xBins, yBins, "DSSD position fission");
    DeclareHistogram2D(DD_IMPLANT_POSITION, 
		       xBins, yBins, "DSSD position implant");
    DeclareHistogram2D(DD_DECAY_POSITION, 
		       xBins, yBins, "DSSD position decay");

    /** Trace Information **/
    DeclareHistogram2D(DD_ENERGY_POS_X_TRACE,
		       energyBins, xBins, "DSSD E vs X from Traces");
    DeclareHistogram2D(DD_ENERGY_POS_Y_TRACE,
		       energyBins, yBins, "DSSD E vs Y from Traces");
    DeclareHistogram2D(DD_FRONTE__BACKE, energyBins2, energyBins2,
            "Front vs Back energy (calib / 100)");

    /** Check Gain match alpha region events **/
    DeclareHistogram2D(DD_EVENT_ENERGY__X_POSITION_IMP,
		       energyBins, xBins, "DSSD X strips E vs. position");
    DeclareHistogram2D(DD_EVENT_ENERGY__Y_POSITION_IMP,
		       energyBins, yBins, "DSSD Y strips E vs. position");
    DeclareHistogram2D(DD_EVENT_ENERGY__X_POSITION_DEC,
		       energyBins, xBins, "DSSD X strips E vs. position");
    DeclareHistogram2D(DD_EVENT_ENERGY__Y_POSITION_DEC,
		       energyBins, yBins, "DSSD Y strips E vs. position");
    /** Check Gain match for HE events **/
    DeclareHistogram2D(DD_EVENT_ENERGY_COMP_X_POSITION_IMP,
		       energyBins2, xBins, "DSSD X strips E/100 vs. position");
    DeclareHistogram2D(DD_EVENT_ENERGY_COMP_Y_POSITION_IMP,
		       energyBins2, yBins, "DSSD Y strips E/100 vs. position");
    DeclareHistogram2D(DD_EVENT_ENERGY_COMP_X_POSITION_DEC,
		       energyBins2, xBins, "DSSD X strips E/100 vs. position");
    DeclareHistogram2D(DD_EVENT_ENERGY_COMP_Y_POSITION_DEC,
		       energyBins2, yBins, "DSSD Y strips E/100 vs. position");

    /** Check Gain match spectra via MaxEvents routine **/

    DeclareHistogram2D(DD_MAXEVENT_ENERGY__X_POSITION_IMP,
		       energyBins, xBins, "MAXDSSD X strips E vs. position");
    DeclareHistogram2D(DD_MAXEVENT_ENERGY__Y_POSITION_IMP,
		       energyBins, yBins, "MAXDSSD Y strips E vs. position");
    DeclareHistogram2D(DD_MAXEVENT_ENERGY__X_POSITION_DEC,
		       energyBins, xBins, "MAXDSSD X strips E vs. position");
    DeclareHistogram2D(DD_MAXEVENT_ENERGY__Y_POSITION_DEC,
		       energyBins, yBins, "MAXDSSD Y strips E vs. position");

    DeclareHistogram2D(DD_MAXEVENT_ENERGY_COMP_X_POSITION_IMP,
		       energyBins2, xBins, "MAXDSSD X strips E/100 vs. position");
    DeclareHistogram2D(DD_MAXEVENT_ENERGY_COMP_Y_POSITION_IMP,
		       energyBins2, yBins, "MAXDSSD Y strips E/100 vs. position");
    DeclareHistogram2D(DD_MAXEVENT_ENERGY_COMP_X_POSITION_DEC,
		       energyBins2, xBins, "MAXDSSD X strips E/100 vs. position");
    DeclareHistogram2D(DD_MAXEVENT_ENERGY_COMP_Y_POSITION_DEC,
		       energyBins2, yBins, "MAXDSSD Y strips E/100 vs. position");


    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 0, energyBins, timeBins,
		       "DSSD Ty,Ex (10ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 1, energyBins, timeBins,
		       "DSSD Ty,Ex (100ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 2, energyBins, timeBins,
		       "DSSD Ty,Ex (400ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 3, energyBins, timeBins,
		       "DSSD Ty,Ex (1us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 4, energyBins, timeBins,
		       "DSSD Ty,Ex (10us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 5, energyBins, timeBins,
		       "DSSD Ty,Ex (100us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 6, energyBins, timeBins,
		       "DSSD Ty,Ex (1ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 7, energyBins, timeBins,
		       "DSSD Ty,Ex (10ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 8, energyBins, timeBins,
		       "DSSD Ty,Ex (100ms/ch)(xkeV)");
    /** Diagnostics **/ 
    DeclareHistogram2D(DD_ENERGY__POSX_T_MISSING,
		       energyBins, xBins, "DSSD T missing X strips E vs. position");
    DeclareHistogram2D(DD_ENERGY__POSY_T_MISSING,
		       energyBins, yBins, "DSSD T missing Y strips E vs. position");

    /** Check how many strips and how far fired **/
    DeclareHistogram2D(DD_DENERGY__DPOS_X_CORRELATED,
		       energyBins, xBins, "DSSD dE dX correlated events");
    DeclareHistogram2D(DD_DENERGY__DPOS_Y_CORRELATED,
		       energyBins, yBins, "DSSD dE dY correlated events");
    /** Pixel Correlated Events **/
    DeclareHistogram2D(DD_CHAIN_NUM_FISSION, SB, SA, 
            "Event Number vs. Fission Chain Energy/10");
    DeclareHistogram2D(DD_CHAIN_NUM_ALPHA, SB, SA, 
            "Event Number vs. Alpha Chain Energy");
    DeclareHistogram2D(DD_CHAIN_ALPHA_V_ALPHA, SB, SB,
	    "First Alpha vs. Alpha Chain Energy");
}


bool Dssd4SHEProcessor::PreProcess(RawEvent &event) {
    if (!EventProcessor::PreProcess(event))
        return false;

    xyEventsTMatch_.clear();
    xyEventsEMatch_.clear();

    vector<ChanEvent*> xEvents = 
        event.GetSummary("dssd_back:dssd_back", true)->GetList();
    vector<ChanEvent*> yEvents = 
        event.GetSummary("dssd_front:dssd_front", true)->GetList();

    /**
     * Matching the front-back by the time correlations
     */
    vector< pair<StripEvent, bool> > xEventsTMatch;
    vector< pair<StripEvent, bool> > yEventsTMatch;
    StripEvent ev2x;
    StripEvent ev2y;
    for (vector<ChanEvent*>::iterator itx = xEvents.begin();
         itx != xEvents.end();
         ++itx) {
        StripEvent ev((*itx)->GetCalEnergy(), 
                      (*itx)->GetTime(),
                      (*itx)->GetChanID().GetLocation(),
                      (*itx)->IsSaturated());
        pair<StripEvent, bool> match(ev, false);
        xEventsTMatch.push_back(match);

        const Trace& trace = (*itx)->GetTrace();
	
        /** Handle additional pulses (no. 2, 3, ...) */
        int pulses = trace.GetValue("numPulses");
        for (int i = 1; i < pulses; ++i) {
            stringstream energyCalName;
            energyCalName << "filterEnergy" << i + 1 << "Cal";
            stringstream timeName;
            timeName << "filterTime" << i + 1;

            ev.pileup = true;

            
            ev2x.E = trace.GetValue(energyCalName.str());
            ev2x.t = (trace.GetValue(timeName.str()) - 
                     trace.GetValue("filterTime") + ev.t);
            ev2x.pos = ev.pos;
            ev2x.sat = false;
            ev2x.pileup = true;
	    
            pair<StripEvent, bool> match2(ev2x, false);
            xEventsTMatch.push_back(match2);
	    
           /* if (i > 1 && ev2x.E > 0) {
                stringstream ss;
                ss << "DSSD X, " << i + 1 << " pulse"
                << ", E = " << ev2x.E
                << ", dt = " << ev2x.t - ev.t;
    	  	
                Messenger m;
                m.run_message(ss.str());
            } */
        }

        for (vector<ChanEvent*>::iterator itx2 = itx;
            itx2 != xEvents.end();
            ++itx2) {
            int dx = abs( ev.pos -  (*itx2)->GetChanID().GetLocation());
            double dE = abs( ev.E - (*itx2)->GetCalEnergy());
	 // plot(DD_DENERGY__DPOS_X_CORRELATED, dE, dx);
	   
        }
    }

    for (vector<ChanEvent*>::iterator ity = yEvents.begin();
         ity != yEvents.end();
         ++ity) {
        StripEvent ev((*ity)->GetCalEnergy(), 
                      (*ity)->GetTime(),
                      (*ity)->GetChanID().GetLocation(),
                      (*ity)->IsSaturated());
        pair<StripEvent, bool> match(ev, false);
        yEventsTMatch.push_back(match);

        const Trace& trace = (*ity)->GetTrace();

        int pulses = trace.GetValue("numPulses");

        for (int i = 1; i < pulses; ++i) {
            stringstream energyCalName;
            energyCalName << "filterEnergy" << i + 1 << "Cal";
            stringstream timeName;
            timeName << "filterTime" << i + 1;

            ev.pileup = true;


            ev2y.E = trace.GetValue(energyCalName.str());
            ev2y.t = (trace.GetValue(timeName.str()) - 
                     trace.GetValue("filterTime") + ev.t);
            ev2y.pos = ev.pos;
            ev2y.sat = false;
            ev2y.pileup = true;
            pair<StripEvent, bool> match2(ev2y, false);
            yEventsTMatch.push_back(match2);
	    
            if (i > 1 && abs(1-ev2x.E/ev2y.E) < 0.3) {
               // stringstream ss;
               // ss << "DSSD Y, " << i + 1 << " pulse"
               // << ", E = " << ev2y.E <<  ", Ex = " << ev2x.E
               // << ", dt = " << ev2y.t - ev.t <<" , dtx = " << ev2x.t - ev.t;
                plot(DD_ENERGY_POS_Y_TRACE, ev2y.E, ev2y.pos);
		plot(DD_ENERGY_POS_X_TRACE, ev2x.E, ev2x.pos);
		//plot(DD_DENERGY__DPOS_X_CORRELATED, ev2x.E, ev2x.t-ev2y.t+20);
		//plot(DD_DENERGY__DPOS_Y_CORRELATED, abs(ev2y.E-ev2x.E), ev2x.t-ev2y.t+20);

               // Messenger m;
               // m.run_message(ss.str());
            }
        }

        for (vector<ChanEvent*>::iterator ity2 = ity;
            ity2 != yEvents.end();
            ++ity2) {
            int dy = abs( ev.pos  - (*ity2)->GetChanID().GetLocation());
            double dE = abs( ev.E - (*ity2)->GetCalEnergy());
  	   // plot(DD_DENERGY__DPOS_Y_CORRELATED, dE, dy);
            
        }
    }

    for (vector< pair<StripEvent, bool> >::iterator itx = xEventsTMatch.begin();
         itx != xEventsTMatch.end();
         ++itx) {
        double bestDtime = numeric_limits<double>::max();
        vector< pair<StripEvent, bool> >::iterator bestMatch =
            yEventsTMatch.end();
        for (vector< pair<StripEvent, bool> >::iterator ity = 
                                                     yEventsTMatch.begin();
            ity != yEventsTMatch.end();
            ++ity) 
        {
            // If already matched, skip
            if ((*ity).second)
                continue;

            double energyX = (*itx).first.E;
            double energyY = (*ity).first.E;

            /** If energies are in lower range and/or not satured
             *  check if delta energy condition is not met, 
             *  if not, skip this event
             *
             *  For high energy events and satured set 20 MeV
             *  energy for difference check. The calibration in this
             *  range is most likely imprecise, so one cannot correlate
             *  by energy difference.
             **/
            if ( (*itx).first.sat || energyX > highEnergyCut_ )
                energyX = 20000.0;
            if ( (*ity).first.sat || energyY > highEnergyCut_ )
                energyY = 20000.0;
            if ( abs(energyX - energyY) > deltaEnergy_)
                continue;

            double dTime = abs((*itx).first.t - (*ity).first.t) *
	      Globals::get()->clockInSeconds();
            if (dTime < bestDtime) {
                bestDtime = dTime;
                bestMatch = ity;
            }
        }

	
        if (bestDtime < timeWindow_) {
            xyEventsTMatch_.push_back(
                pair<StripEvent, StripEvent>((*itx).first, (*bestMatch).first));
            (*itx).second = true;
            (*bestMatch).second = true;
            plot(D_DTIME, int(bestDtime / 1.0e-8) + 1);
        } else {
            bestDtime = int(bestDtime / 1.0e-8);
            if (bestDtime > S8)
                bestDtime = S8 - 1;
            else if (bestDtime < 0)
                bestDtime = 0;
            plot(D_DTIME, bestDtime);
        }
    }
    double ev2xpos, ev2ypos, ev2yE, ev2xE, ev2xt, ev2yt;
    for (vector< pair<StripEvent, bool> >::iterator itx = xEventsTMatch.begin();
         itx != xEventsTMatch.end();
         ++itx) {
        if ((*itx).second)
            continue;
        int position = (*itx).first.pos;
            ev2xpos=position;
        double energy = (*itx).first.E;
            ev2xE = energy;
            ev2xt = (*itx).first.t;
        plot(DD_ENERGY__POSX_T_MISSING, energy, position);
    }

    for (vector< pair<StripEvent, bool> >::iterator ity = yEventsTMatch.begin();
         ity != yEventsTMatch.end();
         ++ity) {
        if ((*ity).second)
            continue;
        int position = (*ity).first.pos;
            ev2ypos=position;
        double energy = (*ity).first.E;
            ev2yE=energy;
            ev2yt = (*ity).first.t;
            plot(DD_ENERGY__POSY_T_MISSING, energy, position);
    }
     

    /**
     * Matching the front-back by the Energy of the event
     * Using the old style GetMaxEvent for comparison
     */
    if (xEvents.size() > 0 && yEvents.size() > 0) {
            ChanEvent* maxFront =
                event.GetSummary("dssd_back:dssd_back")->GetMaxEvent(true);
            ChanEvent* maxBack = 
                event.GetSummary("dssd_front:dssd_front")->GetMaxEvent(true);
            StripEvent evf(maxFront->GetCalEnergy(), 
                           maxFront->GetTime(),
                           maxFront->GetChanID().GetLocation(),
                           maxFront->IsSaturated());
            StripEvent evb(maxBack->GetCalEnergy(), 
                           maxBack->GetTime(),
                           maxBack->GetChanID().GetLocation(),
                           maxBack->IsSaturated());
        xyEventsEMatch_.push_back(pair<StripEvent, StripEvent>(evf, evb));
    }

    return true; 
}


bool Dssd4SHEProcessor::Process(RawEvent &event)
{
    using namespace dammIds::dssd4she;
    //           cout << "test" << endl;
    if (!EventProcessor::Process(event))
        return false;
    
    vector<ChanEvent*> vetoEvents = 
        event.GetSummary("si:veto", true)->GetList();
    vector<ChanEvent*> sideEvents = 
        event.GetSummary("si:si", true)->GetList();
    vector<ChanEvent*> mwpcEvents = 
        event.GetSummary("mcp:mcp", true)->GetList();
    int mwpc = event.GetSummary("mcp", true)->GetMult();
    bool hasBeam = TreeCorrelator::get()->place("Beam")->status();

   plot(D_MWPC_MULTI, mwpc);

   
    for (vector< pair<StripEvent, StripEvent> >::iterator it =
                                                 xyEventsTMatch_.begin();
         it != xyEventsTMatch_.end(); ++it)
    {
        double xEnergy = (*it).first.E;
        double yEnergy = (*it).second.E;
	int hasPileup = (*it).first.pileup;
        /** If saturated set to 200 MeV **/
        if ((*it).first.sat && !(*it).second.sat) {
	    plot(D_BACKE_W_FRONT_SAT, yEnergy);
            xEnergy = yEnergy;
	} else if (!(*it).first.sat && (*it).second.sat) {
	    plot(D_BACKE_W_FRONT_SAT, xEnergy);
            yEnergy = xEnergy;
	} else if ((*it).first.sat && (*it).second.sat) {
            xEnergy = 250000.0;
            yEnergy = 250000.0;
        }
	
	
	

        int xPosition = (*it).first.pos;
        int yPosition = (*it).second.pos;

        double time = min((*it).first.t, (*it).second.t);
	
	plot(D_ENERGY_X, xEnergy / 10.0); 
        plot(D_ENERGY_Y, yEnergy / 10.0);
        plot(DD_FRONTE__BACKE, xEnergy / 100.0, yEnergy / 100.0);
        plot(DD_EVENT_POSITION, xPosition, yPosition);
	
        double static mwpcTime;
        double MwpcE;

	for (vector<ChanEvent*>::iterator itm = mwpcEvents.begin();
	     itm != mwpcEvents.end();
	     ++itm) {
/*            double dt = (time - (*itm)->GetTime()) *
                     Globals::get()->clockInSeconds();
           // if (dt < mwpcTime) 
                mwpcTime = dt;*/
		mwpcTime = (*itm)->GetTime();
		
		MwpcE = (*itm) -> GetCalEnergy();
	}
        // Plot up to 3 us only
       /* if (mwpcTime < 3.0e-6) {
            int timeBin = int(mwpcTime / 1.0e-8);
            int energyBin = xEnergy / 100.0;
            plot(DD_ENERGY_DT__DSSD_MWPC, energyBin, timeBin);
        }*/
        if (vetoEvents.size() > 0) {
	  for (vector<ChanEvent*>::iterator itv = vetoEvents.begin();
	       itv != vetoEvents.end();
	       ++itv) {
	    double vetoEnergy = (*itv)->GetCalEnergy();
	    plot(DD_DE_E__DSSD_VETO, (vetoEnergy + xEnergy) / 100.0, xEnergy / 100.0);
	  }
        }
	
        double bestSiTime = numeric_limits<double>::max();
        ChanEvent* correlatedSide = 0;
        bool hasEscape = false;
        double escapeEnergy = 0.0;
        for (vector<ChanEvent*>::iterator its = sideEvents.begin();
            its != sideEvents.end();
            ++its) {
            double dt = abs(time - (*its)->GetTime()) *
                        Globals::get()->clockInSeconds();
            if (dt < bestSiTime) {
                bestSiTime = dt;
                correlatedSide = *its;
            }
        }

        if (correlatedSide != 0) {
            int siTime = int(bestSiTime / 1.0e-8) + 1;
            if (siTime > S8)
                siTime = S8 - 1;
            else if (siTime < 0)
                siTime = 0;
            plot(D_DTIME_SIDE, siTime);

            if (bestSiTime < timeWindow_) {
                plot(D_ENERGY_CORRELATED_SIDE, correlatedSide->GetCalEnergy());
                hasEscape = true;
                escapeEnergy = correlatedSide->GetCalEnergy();
            }
        }
	/** if pileup second signal has no mwpc **/
	/*if (hasPileup) {
	    mwpc=0; //Can be done better.
        }*/



        bool hasVeto = false;
        if (vetoEvents.size() > 0)
            hasVeto = true;

        if (hasVeto) {
            plot(D_ENERGY_WITH_VETO, xEnergy / 10.0);
	}
        if (mwpc > 0 && !hasVeto) {
            plot(D_ENERGY_WITH_MWPC, xEnergy / 10.0);
            plot(D_MWPC_ENERGY_WITH_DSSD, MwpcE);
	}
        if (hasVeto && mwpc > 0) {
	    plot(D_ENERGY_WITH_VETO_MWPC, xEnergy / 10.0);
	}
        if (!hasVeto && mwpc == 0) {
            plot(D_ENERGY_NO_VETO_MWPC, xEnergy / 10.0);
	}
        
	
	double tCalEnergy=xEnergy;
	double tCalRes=10e-9;
	double tCaltime=((time-mwpcTime)*Globals::get()->clockInSeconds());
	DetectorDriver* driver = DetectorDriver::get();
	time_t theTime;
	theTime = driver -> GetWallTime(time);

	SheEvent event = SheEvent(tCalEnergy + escapeEnergy, time, mwpc, mwpcTime,
                                  hasBeam, hasVeto, hasEscape, unknown);
        //cout << time << endl;
	pickEventType(event);
	
	if (!event.get_beam())  
	plot(D_ENERGY_DECAY_BEAMSTOP, event.get_energy());
	if (event.get_type() == heavyIon) {
	  plot(DD_IMPLANT_POSITION, xPosition, yPosition);
	  plot(D_ENERGY_IMPLANT, event.get_energy());
	  plot(DD_EVENT_ENERGY__X_POSITION_IMP, event.get_energy(), xPosition);
	  plot(DD_EVENT_ENERGY__Y_POSITION_IMP, event.get_energy(), yPosition);
	  plot(DD_EVENT_ENERGY_COMP_X_POSITION_IMP, event.get_energy()/100, xPosition);
	  plot(DD_EVENT_ENERGY_COMP_Y_POSITION_IMP, event.get_energy()/100, yPosition);
	}
        else if (event.get_type() == alpha) {
	  plot(D_ENERGY_DECAY, event.get_energy() / 10.0);
	  plot(DD_DECAY_POSITION, xPosition, yPosition);
	  plot(DD_EVENT_ENERGY__X_POSITION_DEC,event.get_energy(), xPosition);
	  plot(DD_EVENT_ENERGY__Y_POSITION_DEC, event.get_energy(), yPosition);
	  plot(DD_EVENT_ENERGY_COMP_X_POSITION_DEC, event.get_energy()/100, xPosition);
	  plot(DD_EVENT_ENERGY_COMP_Y_POSITION_DEC, event.get_energy()/100, yPosition);
	}
        else if (event.get_type() == lightIon) {
	  plot(DD_LIGHT_POSITION, xPosition, yPosition);
	  plot(D_ENERGY_LIGHT, event.get_energy() / 10.0);
	  plot(DD_EVENT_ENERGY__X_POSITION_IMP, event.get_energy(), xPosition);
	  plot(DD_EVENT_ENERGY__Y_POSITION_IMP, event.get_energy(), yPosition);
	  plot(DD_EVENT_ENERGY_COMP_X_POSITION_IMP, event.get_energy()/100, xPosition);
	  plot(DD_EVENT_ENERGY_COMP_Y_POSITION_IMP, event.get_energy()/100, yPosition);
	}
        else if (event.get_type() == unknown) {
	  plot(DD_UNKNOWN_POSITION, xPosition, yPosition);
	  plot(D_ENERGY_UNKNOWN, event.get_energy() / 10.0);
        }
        else if (event.get_type() == fission) {
	  plot(DD_FISSION_POSITION, xPosition, yPosition);
	  plot(D_ENERGY_FISSION, event.get_energy() / 10.0);
	  plot(DD_EVENT_ENERGY__X_POSITION_DEC, event.get_energy(), xPosition);
	  plot(DD_EVENT_ENERGY__Y_POSITION_DEC, event.get_energy(), yPosition);
	  plot(DD_EVENT_ENERGY_COMP_X_POSITION_DEC, event.get_energy()/100, xPosition);
	  plot(DD_EVENT_ENERGY_COMP_Y_POSITION_DEC, event.get_energy()/100, yPosition);
	}
	
	
	if (event.get_type()== alpha) {
	  const unsigned int NumGranularities = 9;
	  // time resolution in seconds per bin
	  const double timeResolution[NumGranularities] = 
	    {10e-9, 100e-9, 300e-9, 1e-6, 10e-6, 100e-6, 1e-3, 10e-3, 100e-3};
	  for (unsigned int i = 0; i < NumGranularities; i++) {
	    double timeBin = (event.get_time()
			      *Globals::get()->clockInSeconds());
	    double energyBin = event.get_energy();
	    //    plot(DD_ENERGY_DECAY_TIME_GRANX + i,energyBin,timeBin);
	    }
        }
        correlator_.add_event(event, xPosition, yPosition, histo);
	
    }

    /** Old style max event for comparison */
    for (vector< pair<StripEvent, StripEvent> >::iterator it =
                                                 xyEventsEMatch_.begin();
         it != xyEventsEMatch_.end();
         ++it) {
	double time = min((*it).first.t, (*it).second.t);
        double xEnergy = (*it).first.E;
        double yEnergy = (*it).second.E;

        int xPosition = (*it).first.pos;
        int yPosition = (*it).second.pos;
	int mwpc = event.GetSummary("mcp", true)->GetMult();
	vector<ChanEvent*> vetoEvents = event.GetSummary("si:veto", true)->GetList();
        bool hasVeto = false;

        if (vetoEvents.size() > 0){
            hasVeto = true;
	}

        plot(DD_EVENT_POSITION_FROM_E, xPosition, yPosition);

        if (mwpc > 0 && !hasVeto) {
            plot(DD_MAXEVENT_ENERGY__X_POSITION_IMP, xEnergy, xPosition);
            plot(DD_MAXEVENT_ENERGY__Y_POSITION_IMP, yEnergy, yPosition);
            plot(DD_MAXEVENT_ENERGY_COMP_X_POSITION_IMP, xEnergy/100, xPosition);
            plot(DD_MAXEVENT_ENERGY_COMP_Y_POSITION_IMP, yEnergy/100, yPosition);
	}
        if (!hasVeto && mwpc == 0) {
            plot(DD_MAXEVENT_ENERGY__X_POSITION_DEC, xEnergy, xPosition);
            plot(DD_MAXEVENT_ENERGY__Y_POSITION_DEC, yEnergy, yPosition);
            plot(DD_MAXEVENT_ENERGY_COMP_X_POSITION_DEC, xEnergy/100, xPosition);
            plot(DD_MAXEVENT_ENERGY_COMP_Y_POSITION_DEC, yEnergy/100, yPosition);
	}

    }        

    EndProcess();
    return true;
}

bool Dssd4SHEProcessor::pickEventType(SheEvent& event) {
    /**
     * Logic table (V - veto, M - mwpc, B - beam )
     * Logic state is converted into a numerical value N
     * like a binary number:
     *
     * V M B | N | decision
     * --------------------
     * 0 0 0 | 0 | unknown / alpha / fission (depending on energy)
     * 0 0 1 | 1 | -"-
     * 0 1 0 | 2 | unknown
     * 0 1 1 | 3 | heavyIon
     * 1 0 0 | 4 | unknown
     * 1 0 1 | 5 | lightIon
     * 1 1 0 | 6 | unknown
     * 1 1 1 | 7 | lightIon
     *
     **/
    int condition = 0;
    if (event.get_beam()) 
        condition += 1;
    if (event.get_mwpc() > 0) 
        condition += 2;
    if (event.get_veto()) 
        condition += 4;

    if (condition == 0) {
        double energy = event.get_energy();
        if (energy < highEnergyCut_)
            event.set_type(alpha);
        else if (energy < fissionEnergyCut_)
            event.set_type(unknown);
        else 
            event.set_type(fission);
    } 
    else if (condition == 1) {
        double energy = event.get_energy();
        if (energy < lowEnergyCut_)
            event.set_type(unknown);
        else if (energy < highEnergyCut_)
            event.set_type(alpha);
        else if (energy < fissionEnergyCut_)
            event.set_type(unknown);
        else 
            event.set_type(fission);
    } 
    else if (condition == 2 || 
             condition == 4 ||
             condition == 6) {
        event.set_type(unknown);
    } 
    else if (condition == 3 && event.get_energy() > recoilEnergyCut_) {
        event.set_type(heavyIon);
    }
    else if (condition == 5 || condition == 7) {
        event.set_type(lightIon);
    }
    else
        event.set_type(unknown);

    return true;
}
