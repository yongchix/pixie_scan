/*! \file Dssd4JAEAProcessor.cpp
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
#include <limits.h>
#include "Dssd4JAEAProcessor.hpp"
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
#include "TraceFilterer.hpp"
#include "WaveformAnalyzer.hpp"
#include "DetectorDriver.hpp"
#include "Notebook.hpp"

using namespace dammIds::dssd4jaea;
using namespace std;


Dssd4JAEAProcessor::Dssd4JAEAProcessor(double timeWindow,
                                       double deltaEnergy,
				       double recoilEnergyCut,
                                       double highEnergyCut,
                                       double lowEnergyCut,
                                       double fissionEnergyCut,
                                       int numBackStrips,
                                       int numFrontStrips) :
  EventProcessor(OFFSET, RANGE, "dssd4jaea"),

    correlator_(numBackStrips, numFrontStrips)
{
    timeWindow_ = timeWindow;
    deltaEnergy_ = deltaEnergy;
    recoilEnergyCut_ = recoilEnergyCut;
    highEnergyCut_ = highEnergyCut;
    lowEnergyCut_ = lowEnergyCut;
    fissionEnergyCut_ = fissionEnergyCut;
    name = "dssd";
    associatedTypes.insert("dssd_front_jaea");
    associatedTypes.insert("dssd_back_jaea");
    numDoubleTraces=0;
    
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


void Dssd4JAEAProcessor::DeclarePlots(void)
{
  using namespace dammIds::dssd;
  //using namespace dammIds::trace::doubletraceanalyzer;
 
  //   TraceFilterer::DeclarePlots();
  
  const int energyBins = SE;  // 16384   
  const int implantEnergyBins = SF; //SE
  const int xBins = S6;       //    64
  const int timeBins = S8;    //   256
  const int decayEnergyBins = SD; // 8192
  const int decayEnergyBins2= SA; // 1024
  const int positionBins = S6;
  const int posBins = SE;
  unsigned short numTraces = Globals::get()->numTraces();
  
  // Diagnostic plot in PreProcessor 
  // 700 : DSSD decay pattern
  DeclareHistogram2D(DD_DEVENT_POSITION, 
		     xBins, xBins, "DSSD decay pattern");
  // 701 : DSSD implant pattern
  DeclareHistogram2D(DD_IEVENT_POSITION, 
		     xBins, xBins, "DSSD implant pattern");
  // 702  DSSD ch vs decayEnergy at front
  DeclareHistogram2D(DD_DSSDDFRONT_POSENERGY, 
		     decayEnergyBins,xBins, "DSSD front ch vs decay energy");
  // 703  Dssd ch vs decayEnergy at back
  DeclareHistogram2D(DD_DSSDDBACK_POSENERGY, 
		     decayEnergyBins,xBins, "DSSD back ch vs decay energy");
  // 704 : Dssd ch vs implant Energy at front
  DeclareHistogram2D(DD_DSSDIFRONT_POSENERGY, 
		     energyBins, xBins, "DSSD front ch vs implant energy");
  // 705 : Dssd ch vs implant Energy at back
  DeclareHistogram2D(DD_DSSDIBACK_POSENERGY, 
		     energyBins, xBins, "DSSD back pos vs implant energy");
  // 725 checked
  DeclareHistogram2D(DD_IMPLANT_POSITION, 
		     positionBins, positionBins, "DSSD Hit Pattern - R");
  // 726 checked
  DeclareHistogram2D(DD_DECAY_POSITION, 
		     positionBins, positionBins, "DSSD Hit Pattern - D");
  // 741 checked
  DeclareHistogram2D(DD_IMPLANT_FRONT_ENERGY__POSITION,  
		     implantEnergyBins, positionBins, "DSSD Strip vs E - RF");
  // 742 checked 
    DeclareHistogram2D(DD_IMPLANT_BACK_ENERGY__POSITION,
		       implantEnergyBins, positionBins, "DSSD Strip vs E - RB");
    // 743 checked 
    DeclareHistogram2D(DD_DECAY_FRONT_ENERGY__POSITION,
		       decayEnergyBins, positionBins, "DSSD Strip vs E - DF");
    // 744 checked
    DeclareHistogram2D(DD_DECAY_BACK_ENERGY__POSITION, 
		       decayEnergyBins, positionBins, "DSSD Strip vs E - DB");
    
    // 750-759 
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 0, decayEnergyBins, timeBins,
		       "DSSD Ty,Ex (10ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 1, decayEnergyBins, timeBins, 
		       "DSSD Ty,Ex (100ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 2, decayEnergyBins, timeBins,
		       "DSSD Ty,Ex (400ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 3, decayEnergyBins, timeBins,
		       "DSSD Ty,Ex (1us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 4, decayEnergyBins, timeBins, 
		       "DSSD Ty,Ex (100us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 5, decayEnergyBins, timeBins,
		       "DSSD Ty,Ex (1ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 6, decayEnergyBins, timeBins,
		       "DSSD Ty,Ex (10ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 7, decayEnergyBins, timeBins,
		       "DSSD Ty,Ex (100ms/ch)(xkeV)");//ntb
    
    // 850-859 GRANX with VETO                                                        
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_VETO + 0, decayEnergyBins, timeBins,
                       "DSSD_VETO Ty,Ex (10ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_VETO + 1, decayEnergyBins, timeBins
		       ,                         "DSSD_VETO Ty,Ex (100ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_VETO + 2, decayEnergyBins, timeBins
		       ,
                       "DSSD_VETO Ty,Ex (400ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_VETO + 3, decayEnergyBins, timeBins
		       ,
                       "DSSD_VETO Ty,Ex (1us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_VETO + 4, decayEnergyBins, timeBins
		       ,
                       "DSSD_VETO Ty,Ex (100us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_VETO + 5, decayEnergyBins, timeBins
		       ,
                       "DSSD_VETO Ty,Ex (1ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_VETO + 6, decayEnergyBins, timeBins
		       ,
                       "DSSD_VETO Ty,Ex (10ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_VETO + 7, decayEnergyBins, timeBins
		       ,
                       "DSSD_VETO Ty,Ex (100ms/ch)(xkeV)");

    // 860-861 GRANX with VETO                                                        
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NOVETO + 0, decayEnergyBins, timeBins,
		       "DSSD_NOVETO Ty,Ex (10ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NOVETO + 1, decayEnergyBins, timeBins,
		       "DSSD_NOVETO Ty,Ex (100ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NOVETO + 2, decayEnergyBins, timeBins,
                       "DSSD_NOVETO Ty,Ex (400ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NOVETO + 3, decayEnergyBins, timeBins,
                       "DSSD_NOVETO Ty,Ex (1us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NOVETO + 4, decayEnergyBins, timeBins,
		       "DSSD_NOVETO Ty,Ex (100us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NOVETO + 5, decayEnergyBins, timeBins,
                       "DSSD_NOVETO Ty,Ex (1ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NOVETO + 6, decayEnergyBins, timeBins,
                       "DSSD_NOVETO Ty,Ex (10ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NOVETO + 7, decayEnergyBins, timeBins,
                       "DSSD_NOVETO Ty,Ex (100ms/ch)(xkeV)");

    
    // 860-861 GRANX with NAI
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NAI + 0, decayEnergyBins, timeBins,
		       "DSSD_NAI Ty,Ex (10ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NAI + 1, decayEnergyBins, timeBins,
                       "DSSD_NAI Ty,Ex (100ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NAI + 2, decayEnergyBins, timeBins,
                       "DSSD_NAI Ty,Ex (400ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NAI + 3, decayEnergyBins, timeBins,
		       "DSSD_NAI Ty,Ex (1us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NAI + 4, decayEnergyBins, timeBins,
		       "DSSD_NAI Ty,Ex (100us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NAI + 5, decayEnergyBins, timeBins,
		       "DSSD_NAI Ty,Ex (1ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NAI + 6, decayEnergyBins, timeBins,
                       "DSSD_NAI Ty,Ex (10ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NAI + 7, decayEnergyBins, timeBins,
                       "DSSD_NAI Ty,Ex (100ms/ch)(xkeV)");

    // 860-861 GRANX with no_NAI
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NONAI + 0, decayEnergyBins, timeBins,
		       "DSSD_NONAI Ty,Ex (10ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NONAI + 1, decayEnergyBins, timeBins,
                       "DSSD_NONAI Ty,Ex (100ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NONAI + 2, decayEnergyBins, timeBins,
                       "DSSD_NONAI Ty,Ex (400ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NONAI + 3, decayEnergyBins, timeBins,
		       "DSSD_NONAI Ty,Ex (1us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NONAI + 4, decayEnergyBins, timeBins,
		       "DSSD_NONAI Ty,Ex (100us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NONAI + 5, decayEnergyBins, timeBins,
		       "DSSD_NONAI Ty,Ex (1ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NONAI + 6, decayEnergyBins, timeBins,
                       "DSSD_NONAI Ty,Ex (10ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NONAI + 7, decayEnergyBins, timeBins,
                       "DSSD_NONAI Ty,Ex (100ms/ch)(xkeV)");


    // 860-861 GRANX with PIN
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_PIN + 0, decayEnergyBins, timeBins,
		       "DSSD_PIN Ty,Ex (10ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_PIN + 1, decayEnergyBins, timeBins,
                       "DSSD_PIN Ty,Ex (100ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_PIN + 2, decayEnergyBins, timeBins,
                       "DSSD_PIN Ty,Ex (400ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_PIN + 3, decayEnergyBins, timeBins,
		       "DSSD_PIN Ty,Ex (1us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_PIN + 4, decayEnergyBins, timeBins,
		       "DSSD_PIN Ty,Ex (100us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_PIN + 5, decayEnergyBins, timeBins,
		       "DSSD_PIN Ty,Ex (1ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_PIN + 6, decayEnergyBins, timeBins,
                       "DSSD_PIN Ty,Ex (10ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_PIN + 7, decayEnergyBins, timeBins,
                       "DSSD_PIN Ty,Ex (100ms/ch)(xkeV)");

    // 860-861 GRANX with noPIN
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NOPIN + 0, decayEnergyBins, timeBins,
		       "DSSD_NOPIN Ty,Ex (10ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NOPIN + 1, decayEnergyBins, timeBins,
                       "DSSD_NOPIN Ty,Ex (100ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NOPIN + 2, decayEnergyBins, timeBins,
                       "DSSD_NOPIN Ty,Ex (400ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NOPIN + 3, decayEnergyBins, timeBins,
		       "DSSD_NOPIN Ty,Ex (1us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NOPIN + 4, decayEnergyBins, timeBins,
		       "DSSD_NOPIN Ty,Ex (100us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NOPIN + 5, decayEnergyBins, timeBins,
		       "DSSD_NOPIN Ty,Ex (1ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NOPIN + 6, decayEnergyBins, timeBins,
                       "DSSD_NOPIN Ty,Ex (10ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_NOPIN + 7, decayEnergyBins, timeBins,
                       "DSSD_NOPIN Ty,Ex (100ms/ch)(xkeV)");


    

    
    // 759 checked
    DeclareHistogram1D(D_DECAY_FRONT_ENERGY,
                       decayEnergyBins, "DSSD Energy Front");  
    // 760 checked
    DeclareHistogram1D(D_DECAY_BACK_ENERGY,
                       decayEnergyBins, "DSSD Energy Back");
    // 761 
    DeclareHistogram2D(DD_IMPLANT_POSITION__MCP,
		       posBins, positionBins, "DSSD Front Strip vs MCP_Tac");
    // ** New Histograms for JAEA experiment (Dec 2014)
    // 763
    DeclareHistogram1D(D_DECAY_FRONT_ENERGY_NOVETO,
                       decayEnergyBins, "DSSD Front with noVETO");    
    DeclareHistogram1D(D_DECAY_BACK_ENERGY_NOVETO,
                       decayEnergyBins, "DSSD Back with noVETO");
    // 765
    DeclareHistogram1D(D_DECAY_FRONT_ENERGY_VETO,
                       decayEnergyBins, "DSSD Front with VETO");  
    DeclareHistogram1D(D_DECAY_BACK_ENERGY_VETO,
                       decayEnergyBins, "DSSD Back with VETO");
    // 767
    DeclareHistogram1D(D_DECAY_FRONT_ENERGY_NAI,
                       decayEnergyBins, "DSSD Front with NaI");   
    DeclareHistogram1D(D_DECAY_BACK_ENERGY_NAI,
                       decayEnergyBins, "DSSD Back with NaI");
    // 769
    DeclareHistogram1D(D_DECAY_FRONT_ENERGY_PIN,
                       decayEnergyBins, "DSSD Front with PIN");   
    DeclareHistogram1D(D_DECAY_BACK_ENERGY_PIN,
                       decayEnergyBins, "DSSD Back with PIN");
    // 771
    DeclareHistogram1D(D_DECAY_FRONT_ENERGY_NONAI,
                       decayEnergyBins, "DSSD Front with noNaI");   
    DeclareHistogram1D(D_DECAY_BACK_ENERGY_NONAI,
                       decayEnergyBins, "DSSD Back with noNaI");
    
    // 778
    DeclareHistogram1D(D_DECAY_FRONT_ENERGY_NOPIN,
                       decayEnergyBins, "DSSD Front with noPIN");   
    DeclareHistogram1D(D_DECAY_BACK_ENERGY_NOPIN,
                       decayEnergyBins, "DSSD Back with noPIN");
    
    

    // 773 
    DeclareHistogram2D(DD_ENERGY2F_DT, decayEnergyBins, decayEnergyBins2, 		       "DSSD Ty,Ex traces (10ns/ch)(xkeV)");
    // 774 
    DeclareHistogram2D(DD_ENERGY2F_ENERGY1F, decayEnergyBins, decayEnergyBins2, 
		       "E1 vs E2");
    // 775
    DeclareHistogram2D(DD_ENERGY2B_ENERGY1B, decayEnergyBins, decayEnergyBins2, 
		       "E1 vs E2");
    // 776 checked
    DeclareHistogram2D(DD_DOUBLETRACE_FRONT_WITHOUT_MWPC, 
		       SA, numTraces, "Double traces at front with mcp");
    // 777 checked  
    DeclareHistogram2D(DD_DOUBLETRACE_BACK_WITHOUT_MWPC, 
		       SA, numTraces, "Double traces at back with mcp");
    
    // 780 
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_SECOND + 0, decayEnergyBins, timeBins,
		       "2nd  DSSD Ty,Ex (10ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_SECOND + 1, decayEnergyBins, timeBins, 
		       "2nd DSSD Ty,Ex (100ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_SECOND + 2, decayEnergyBins, timeBins,
		       "2nd DSSD Ty,Ex (400ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_SECOND + 3, decayEnergyBins, timeBins,
		       "2nd DSSD Ty,Ex (1us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_SECOND + 4, decayEnergyBins, timeBins, 
		       "2nd DSSD Ty,Ex (100us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_SECOND + 5, decayEnergyBins, timeBins,
		       "2nd DSSD Ty,Ex (1ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_SECOND + 6, decayEnergyBins, timeBins,
		       "2nd DSSD Ty,Ex (10ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_SECOND + 7, decayEnergyBins, timeBins,
		       "2nd DSSD Ty,Ex (100ms/ch)(xkeV)");
    // DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX_SECOND + 8, decayEnergyBins, timeBins,
    //		       "2nd DSSD Ty,Ex (100ms/ch)(xkeV)");

 
    // 788
    DeclareHistogram2D(DD_ENERGY_DECAY12, decayEnergyBins, decayEnergyBins, 		       "Decay 1st & 2nd front ");
    // 789
    DeclareHistogram2D(DD_ENERGY_DECAY12_2, decayEnergyBins, decayEnergyBins, 		       "Decay 1st & 2nd back ");
    
    // 790
    DeclareHistogram2D(DD_MCPPOS_DSSDBACK,decayEnergyBins2,xBins,"position(MWPC vs DssdBack)");
    
    // 791
    DeclareHistogram2D(DD_MCPPOS_DSSDFRONT,decayEnergyBins2,xBins,"position(MWPC vs DssdFront)");
    
    // 792                                                                  
    DeclareHistogram2D(DD_MCPPOS_DSSDENERGY_BACK,decayEnergyBins2,energyBins," MWPC vs ImplantEnergy back");
    // 793                                                                  
    DeclareHistogram2D(DD_MCPPOS_DSSDENERGY_FRONT,decayEnergyBins2,energyBins," MWPC vs ImplantEnergy front");
    // 794                                                         
    DeclareHistogram2D(DD_MCP2D,SA,SA," MWPC gated by DSSD hit");
    
    
    // 795  Implantation pos vs Energy gated by Cs 
  DeclareHistogram2D(DD_IMPLANT_CSGATE, 
		     decayEnergyBins,xBins, "Implantation gated by Cs decay");
  // 796 Implantation pos vs Energy gated by I
  DeclareHistogram2D(DD_IMPLANT_IGATE, 
		     decayEnergyBins,xBins, "Implantation gated by I decay");



}


bool Dssd4JAEAProcessor::PreProcess(RawEvent &event) {
    if (!EventProcessor::PreProcess(event))
        return false;

    xyEventsTMatch_.clear();
    xyEventsEMatch_.clear();
    
    vector<ChanEvent*> xEvents = 
      event.GetSummary("dssd_front_jaea:dssd_front_jaea",true)->GetList();
    vector<ChanEvent*> yEvents = 
      event.GetSummary("dssd_back_jaea:dssd_back_jaea",true)->GetList();
    
    // only for the diagnostic
    double trace_energy1F=0,trace_time1F=0;
    double trace_energy2F=0,trace_time2F=0;
    double trace_energy1B=0,trace_time1B=0;
    double trace_energy2B=0,trace_time2B=0;

    unsigned int frontPos = INT_MAX, backPos = INT_MAX;
    double frontEnergy,backEnergy;
    //double frontTime = 0;
    Trace frontTrace,backTrace;
    
    //    bool hasFront = (event.GetSummary("dssd_front_jaea",true)->GetMult()>0);
    //bool hasBack  = (event.GetSummary("dssd_back_jaea",true)->GetMult()>0);
    
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
                      (*itx)->IsSaturated(),
		      (*itx)->GetTrace());
	pair<StripEvent, bool> match(ev, false);
        xEventsTMatch.push_back(match);

	
       
        const Trace& traceF = (*itx)->GetTrace();
	frontEnergy = ev.E;
	frontPos    = ev.pos;
	
	if(traceF.HasValue("filterEnergy")){
	  trace_energy1F = traceF.GetValue("filterEnergy");
	  trace_time1F   = traceF.GetValue("filterTime");
	}
	if(traceF.HasValue("filterEnergy2")){
	  trace_energy2F = traceF.GetValue("filterEnergy2");
	  trace_time2F   = traceF.GetValue("filterTime2");
	}
       

        /** Handle additional pulses (no. 2, 3, ...) */
        int pulses = traceF.GetValue("numPulses");
        for (int i = 1; i < pulses; ++i) {
            stringstream energyCalName;
            energyCalName << "filterEnergy" << i + 1 << "Cal";
            stringstream timeName;
            timeName << "filterTime" << i + 1;

            ev.pileup = true;

	    ev2x.E = traceF.GetValue(energyCalName.str());
            ev2x.t = (traceF.GetValue(timeName.str()) - 
                     traceF.GetValue("filterTime") + ev.t);
            ev2x.pos = ev.pos;
            ev2x.sat = false;
            ev2x.pileup = true;
	    
            pair<StripEvent, bool> match2(ev2x, false);
            xEventsTMatch.push_back(match2);
	    

        }
	
	for (vector<ChanEvent*>::iterator itx2 = itx;
	     itx2 != xEvents.end();
	     ++itx2) {
	  // int dx = abs( ev.pos -  (*itx2)->GetChanID().GetLocation());
	  // double dE = abs( ev.E - (*itx2)->GetCalEnergy());
	}
	
    }

    for (vector<ChanEvent*>::iterator ity = yEvents.begin();
         ity != yEvents.end();
         ++ity) {
        StripEvent ev((*ity)->GetCalEnergy(), 
                      (*ity)->GetTime(),
                      (*ity)->GetChanID().GetLocation(),
                      (*ity)->IsSaturated(),
		      (*ity)->GetTrace());
        pair<StripEvent, bool> match(ev, false);
        yEventsTMatch.push_back(match);
	
	//
	const Trace& traceB = (*ity)->GetTrace();
	backEnergy = ev.E;
	backPos    = ev.pos;
	if(traceB.HasValue("filterEnergy")){
	  trace_energy1B = traceB.GetValue("filterEnergy");
	  trace_time1B   = traceB.GetValue("filterTime");
	}
	if(traceB.HasValue("filterEnergy2")){
	  trace_energy2B = traceB.GetValue("filterEnergy2");
	  trace_time2B   = traceB.GetValue("filterTime2");
	}
	//

        int pulses = traceB.GetValue("numPulses");

        for (int i = 1; i < pulses; ++i) {
            stringstream energyCalName;
            energyCalName << "filterEnergy" << i + 1 << "Cal";
            stringstream timeName;
            timeName << "filterTime" << i + 1;

            ev.pileup = true;
	    
            ev2y.E = traceB.GetValue(energyCalName.str());
            ev2y.t = (traceB.GetValue(timeName.str()) - 
                     traceB.GetValue("filterTime") + ev.t);
            ev2y.pos = ev.pos;
            ev2y.sat = false;
            ev2y.pileup = true;
            pair<StripEvent, bool> match2(ev2y, false);
            yEventsTMatch.push_back(match2);
	    
	}
	
        for (vector<ChanEvent*>::iterator ity2 = ity;
	     ity2 != yEvents.end();
	     ++ity2) {
	  ///int dy = abs( ev.pos  - (*ity2)->GetChanID().GetLocation());
	  //double dE = abs( ev.E - (*ity2)->GetCalEnergy());
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
	  if ((*ity).second)
	    continue;
	  
	  //double energyX = (*itx).first.E;
	  //double energyY = (*ity).first.E;
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
      } else {
	bestDtime = int(bestDtime / 1.0e-8);
	if (bestDtime > S8)
	  bestDtime = S8 - 1;
	else if (bestDtime < 0)
	  bestDtime = 0;
      }
    }
    
    
      

    if (xEvents.size() > 0 && yEvents.size() > 0) {
      ChanEvent* maxFront =
	event.GetSummary("dssd_front_jaea:dssd_front_jaea")->GetMaxEvent(true);
      ChanEvent* maxBack = 
	event.GetSummary("dssd_back_jaea:dssd_back_jaea")->GetMaxEvent(true);
      
      StripEvent evf(maxFront->GetCalEnergy(), 
		     maxFront->GetTime(),
		     maxFront->GetChanID().GetLocation(),
		     maxFront->IsSaturated(),
		     maxFront->GetTrace());
      StripEvent evb(maxBack->GetCalEnergy(), 
		     maxBack->GetTime(),
		     maxBack->GetChanID().GetLocation(),
		     maxBack->IsSaturated(),
		     maxBack->GetTrace());
      xyEventsEMatch_.push_back(pair<StripEvent, StripEvent>(evf, evb));
      
      }
    


    
    return true; 
}


bool Dssd4JAEAProcessor::Process(RawEvent &event)
{
  
  using namespace dammIds::dssd4jaea;
  
  if (!EventProcessor::Process(event))
        return false;
  

  double cutoffEnergy=6500;

  vector<ChanEvent*> pinEvents = 
    event.GetSummary("pin:pin", true)->GetList();
  vector<ChanEvent*> NaIEvents = 
    event.GetSummary("nai:nai", true)->GetList();
  vector<ChanEvent*> mwpcEvents =
    event.GetSummary("mcp", true)->GetList();
  vector<ChanEvent*> xEvents = 
    event.GetSummary("dssd_front_jaea:dssd_front_jaea",true)->GetList();
  vector<ChanEvent*> yEvents = 
    event.GetSummary("dssd_back_jaea:dssd_back_jaea",true)->GetList();

  static Correlator &corr = event.GetCorrelator();

  bool hasPin = false;
  bool hasMcp  = false;
  bool hasNaI  = false;
  bool hasVETO = false;
  bool hasFront = false;
  bool hasBack  = false;
  
  int mult_pin  = event.GetSummary("pin", true)->GetMult();
  int mult_mwpc = event.GetSummary("mcp", true)->GetMult();
  int mult_nai  = event.GetSummary("nai", true)->GetMult();
  int mult_dssd_front = event.GetSummary("dssd_front_jaea",true)->GetMult();
  int mult_dssd_back  = event.GetSummary("dssd_back_jaea",true)->GetMult();
  int mwpc = event.GetSummary("mcp", true)->GetMult();
 
 
  if(pinEvents.size() > 0){ hasPin=true; }
  if(mult_mwpc > 0){ hasMcp=true;}
  if(mult_nai  > 0){ hasNaI=true;}
  if(hasNaI||hasPin){
    hasVETO=true;
  }
  if(mult_dssd_front>0){hasFront=true;}
  if(mult_dssd_back>0) {hasBack=true;}
  
  if(pinEvents.size()>0 || mult_nai>0){
    //   cout << " Pin " << pinEvents.size() << " NaI " << mult_nai <<endl;
  }

  //******Manual switch for tests **************
  //hasPin = true;
  //hasMcp  = true;
  //hasNaI  = true;
  //hasVETO = true;
  //hasFront = true;
  //hasBack  = true;
  //***********************************************
  

  double mwpcEnergy1=0;
  double mwpcEnergy2=0;
  double mwpcEnergy3=0;
  double mwpcEnergy4=0;
  string mwpcsubtype;
  int mwpcch=0;
  double mwpcPos=0;
  double mwpcCathode=0;
  double static mwpcTime;
  int inc=0;

  
  
  for (vector<ChanEvent*>::const_iterator it = pinEvents.begin();
       it != pinEvents.end(); it++) {
    ChanEvent *chan = *it;

    string subtype   = chan->GetChanID().GetSubtype();
    int number        = chan->GetChanID().GetLocation();
    double calEnergy = chan->GetCalEnergy();
    double pinTime   = chan->GetTime();
    using std::cout;
    using std:: endl;

    //      cout << number << endl;                                                                       
    if(number==0){
      //  cout << "PIN1 " << calEnergy << " " << pinTime << endl; 
    }
    if(number==2){
      //   cout << "PIN2 " << calEnergy <<" " << pinTime <<  endl;
    }
  }




  for (vector<ChanEvent*>::iterator itm = mwpcEvents.begin();
       itm != mwpcEvents.end();
       ++itm) {
    mwpcch=(*itm)->GetChanID().GetLocation();
    if(inc==0){
      mwpcEnergy1=(*itm)->GetCalEnergy();
      mwpcCathode=mwpcEnergy1;
    }
    if(inc==1){
      mwpcEnergy2=(*itm)->GetCalEnergy();
    }
    if(inc==2){
      mwpcEnergy3=(*itm)->GetCalEnergy();
    }
    if(inc==3){
      mwpcEnergy4=(*itm)->GetCalEnergy();
      mwpcPos=mwpcEnergy4/5;
    }
      
    mwpcTime=(*itm)->GetTime();	 
    
    //      mwpcPos=mwpcEnergy[2];
      //cout << inc << " " << mwpcEnergy[inc] << " " << mwpcTime<< endl;
    inc++;
    
  }
  
  if(hasFront && hasBack){
    plot(DD_MCP2D,mwpcPos,mwpcCathode);
  }
  
  
  
  //    cout << " test "<< endl;
  
  unsigned int frontPos=INT_MAX,backPos=INT_MAX;
  double frontEnergy,backEnergy,frontTime=0;
  
  double time;
  
  static int traceNum;
  static int traceNum_with_mcp;


  int loop=0;
  for (vector< pair<StripEvent, StripEvent> >::iterator it =
	 xyEventsTMatch_.begin();
       it != xyEventsTMatch_.end(); ++it)
    {
      int xPosition = 40-(*it).first.pos;
      int yPosition = 40-(*it).second.pos;
      
      //      cout << "loop " << loop << endl;
      loop++;
      double xEnergy = (*it).first.E;
      double yEnergy = (*it).second.E;
      double xTime   = (*it).first.t;
      double yTime   = (*it).second.t;
      Trace& xTrace  = (*it).first.tr;
      Trace& yTrace  = (*it).second.tr;
      int xpulses = xTrace.GetValue("numPulses"); 
      int ypulses = yTrace.GetValue("numPulses"); 
      
      //      xTime=xTime-100000000;
      
      //      cout << xTime <<" " << yTime << endl;
      
      double trace_energy1F=0,trace_time1F=0;
      double trace_energy2F=0,trace_time2F=0;
      double trace_energy1B=0,trace_time1B=0;
      double trace_energy2B=0,trace_time2B=0;
      int hasPileup = (*it).first.pileup;
      time = min((*it).first.t, (*it).second.t);
      


      // overflow cut of DSSD energy
      if(xEnergy<30000){
	hasFront=true;
      }else{
	hasFront=false;
      }
      if(yEnergy<30000){
	hasBack=true;
      }else{
	hasBack=false;
      }
      ///////////////////////////////

      
      if(hasMcp){
	if(hasFront && hasBack){

	  
	  plot(DD_IEVENT_POSITION,xPosition,yPosition);     	
	  plot(DD_DSSDIFRONT_POSENERGY,xEnergy/5,xPosition); 
	  plot(DD_DSSDIBACK_POSENERGY,yEnergy/5,yPosition); 
	  
	  plot(DD_MCPPOS_DSSDBACK,mwpcPos,yPosition);
	  plot(DD_MCPPOS_DSSDFRONT,mwpcPos,xPosition);
	  
	  plot(DD_MCPPOS_DSSDENERGY_BACK,mwpcPos,yEnergy/5);
	  plot(DD_MCPPOS_DSSDENERGY_FRONT,mwpcPos,xEnergy/5);
	  
	  //	  if(xPosition>=23 && xPosition<=25 && yPosition>=11 && yPosition<=13){
	 // if(xPosition==24 && yPosition==12){
	    //    cout << "Implant " <<xPosition << " " << yPosition << " " << xEnergy << " " << yEnergy  << " " <<  xTime << endl;	  
	 // }

	}
      }
      
     
      if(!hasMcp){
	
	//if(xPosition>=23 && xPosition<=25 && yPosition>=11 && yPosition<=13){
		if(xPosition==24 && yPosition==12){
		  //	  cout << "Decay " <<xPosition << " " << yPosition  << " " << xEnergy << " " << yEnergy <<  " " << xTime  <<endl;	  
	}
	

	plot(DD_DEVENT_POSITION, xPosition, yPosition);
	plot(DD_DSSDDFRONT_POSENERGY,xEnergy,xPosition);
	plot(DD_DSSDDBACK_POSENERGY,yEnergy,yPosition);
	
	
	if(pinEvents.size()>0){
	  for (vector<ChanEvent*>::iterator itv = pinEvents.begin();
	       itv != pinEvents.end();
	       ++itv) {
	    double pinEnergy = (*itv)->GetCalEnergy();
	  }
	}
	if(xpulses>1 && ypulses>1){
	  if(xTrace.HasValue("filterEnergy")){
	    trace_energy1F = xTrace.GetValue("filterEnergy");
	    trace_time1F   = xTrace.GetValue("filterTime");
	  };
	  if(xTrace.HasValue("filterEnergy2")){
	    trace_energy2F = xTrace.GetValue("filterEnergy2");
	    trace_time2F   = xTrace.GetValue("filterTime2");
	  };
	  if(yTrace.HasValue("filterEnergy")){
	    trace_energy1B = yTrace.GetValue("filterEnergy");
	    trace_time1B   = yTrace.GetValue("filterTime");
	  };
	  if(yTrace.HasValue("filterEnergy2")){
	    trace_energy2B = yTrace.GetValue("filterEnergy2");
	    trace_time2B   = yTrace.GetValue("filterTime2");
	  };
	  
	  double DT=trace_time2F-trace_time1F;
	  plot(DD_ENERGY2F_DT,trace_energy2F,DT);
	  plot(DD_ENERGY2F_ENERGY1F,trace_energy2F,trace_energy1F);
	  plot(DD_ENERGY2B_ENERGY1B,trace_energy2B,trace_energy1B);
	  
	  cout << "=====================DOUBLE TRACES=== 1F:" << endl;
	  stringstream ss;
	  ss << trace_energy1F << " 2F: " <<trace_energy2F << " 1B: " << trace_energy1B << " 2B: " << trace_energy2B << " xpos: " << xPosition << " ypos: " << yPosition << endl;
	  Notebook::get()->report(ss.str());

	  for(vector<int>::iterator it = xTrace.begin();it != xTrace.end();it++){	   
	    plot(DD_DOUBLETRACE_FRONT_WITHOUT_MWPC,it-xTrace.begin(),traceNum,*it);
	  }
	  for(vector<int>::iterator it = yTrace.begin();it != yTrace.end();it++){	   
	    plot(DD_DOUBLETRACE_BACK_WITHOUT_MWPC,it-yTrace.begin(),traceNum,*it);
	  }
	  traceNum++;
	}
      }
      
	
	bool hasBeam=false;
	bool hasVeto=false;
	bool hasEscape=false;
	
	double tCalEnergy=xEnergy;
	double tCalRes=10e-9;
	double tCaltime=((time-mwpcTime)*Globals::get()->clockInSeconds());
	DetectorDriver* driver = DetectorDriver::get();
	time_t theTime;
	theTime = driver -> GetWallTime(time);
	//JAEAEvent jaeaevent = JAEAEvent(tCalEnergy, time, mwpc, mwpcTime,
	//				    hasBeam, hasVeto, hasEscape, unknown_jaea);
	// try to input energy information in jaea to gate cs decay and see the distribution of implantation
	JAEAEvent jaeaevent = JAEAEvent(xEnergy, yEnergy,time, mwpc, mwpcTime,
					hasBeam, hasVeto, hasEscape, unknown_jaea);

	
	pickEventType(jaeaevent);
	if (jaeaevent.get_type()== alpha_jaea) {
	  const unsigned int NumGranularities = 8;
	  // time resolution in seconds per bin
	  const double timeResolution[NumGranularities] = 
	    {10e-9, 100e-9, 400e-9, 1e-6, 100e-6, 1e-3, 10e-3, 100e-3};
	  for (unsigned int i = 0; i < NumGranularities; i++) {
	    //    double timeBin = (event.get_time()
	    //			      *Globals::get()->clockInSeconds());
	    double timeBin = (jaeaevent.get_time()*Globals::get()->clockInSeconds()/timeResolution[i]);
	    double energyBin = jaeaevent.get_energy();
	    
	    //cout << "GRANX "<< energyBin << " "<<  timeBin << endl;
	    //plot(DD_ENERGY_DECAY_TIME_GRANX + i,2*energyBin,timeBin);
	  }
	}
	
	correlator_.add_event(jaeaevent, xPosition, yPosition, histo,hasVETO,hasNaI,hasPin);
	
    } // xy event match

  
  
  if (xEvents.size() > 0) {
    ChanEvent* maxFront =
      event.GetSummary("dssd_front_jaea:dssd_front_jaea")->GetMaxEvent(true);
    StripEvent evf(maxFront->GetCalEnergy(), 
		   maxFront->GetTime(),
		   40-maxFront-> GetChanID().GetLocation(),
		   maxFront->IsSaturated(),
		   maxFront->GetTrace());
    //frontPos    = 40-evf.pos;
    // frontPos    = evf.
    frontEnergy = evf.E;
    frontTime   = evf.t;
  }else{frontEnergy=0;}
  if (yEvents.size() > 0) {
    ChanEvent* maxBack = 
      event.GetSummary("dssd_back_jaea:dssd_back_jaea")->GetMaxEvent(true);
    StripEvent evb(maxBack->GetCalEnergy(), 
		   maxBack->GetTime(),
		   40-maxBack->GetChanID().GetLocation(),
		   maxBack->IsSaturated(),
		   maxBack->GetTrace());
    //backPos     = 40-evb.pos;
      backEnergy  = evb.E;
  }else{backEnergy=0;}


  Correlator::EEventType type;
  
  // cout <<  corr.GetCondition() << endl;

  if(hasFront && hasBack){    
    if (frontEnergy > cutoffEnergy && backEnergy > cutoffEnergy) {
      type = Correlator::IMPLANT_EVENT;
    } else if (frontEnergy <= cutoffEnergy && backEnergy <= cutoffEnergy &&
	       !hasMcp) {
      type = Correlator::DECAY_EVENT;
    } else {
      type = Correlator::UNKNOWN_TYPE;
    }
    corr.CorrelateOld(event,type,frontPos,backPos,frontTime);
    //cout << type << " " << (type==Correlator::DECAY_EVENT) << endl;
  } else if (hasFront) {
    if (frontEnergy > cutoffEnergy) {
      type = Correlator::IMPLANT_EVENT;
    } else type = Correlator::DECAY_EVENT; 
	
  } else if (hasBack) {
    if (backEnergy > cutoffEnergy) {
       type = Correlator::IMPLANT_EVENT; 
    } else type = Correlator::DECAY_EVENT;  
  } else type = Correlator::UNKNOWN_TYPE ;
  
  
  // plot stuff
  if(type == Correlator::IMPLANT_EVENT){
    if (hasFront)
      {
	plot(DD_IMPLANT_FRONT_ENERGY__POSITION, frontEnergy, frontPos); 
      }
    if (hasBack)
      plot(DD_IMPLANT_BACK_ENERGY__POSITION, backEnergy, backPos);
    if (hasFront && hasBack)
      plot(DD_IMPLANT_POSITION, frontPos, backPos); 
  } else if (type==Correlator::DECAY_EVENT) {
    if (hasFront)
      {
	plot(D_DECAY_FRONT_ENERGY, frontEnergy);
	plot(DD_DECAY_FRONT_ENERGY__POSITION, frontEnergy, frontPos);
	
	if(hasVETO){
	  plot(D_DECAY_FRONT_ENERGY_VETO, frontEnergy);
	  plot(DD_DECAY_FRONT_ENERGY__POSITION_VETO, frontEnergy, frontPos);
	}else if(!hasVETO){
	  plot(D_DECAY_FRONT_ENERGY_NOVETO, frontEnergy);
	  plot(DD_DECAY_FRONT_ENERGY__POSITION_NOVETO, frontEnergy, frontPos);
	}
	if(hasNaI){
	  plot(D_DECAY_FRONT_ENERGY_NAI, frontEnergy);
	}else if(!hasNaI){
	  plot(D_DECAY_FRONT_ENERGY_NONAI, frontEnergy);
	}
	if(hasPin){
	  plot(D_DECAY_FRONT_ENERGY_PIN, frontEnergy);
	}else if(!hasPin){
	  plot(D_DECAY_FRONT_ENERGY_PIN, frontEnergy);
	}
	
	

      }
    if (hasBack)
      {
	plot(D_DECAY_BACK_ENERGY, backEnergy);
	plot(DD_DECAY_BACK_ENERGY__POSITION, backEnergy, backPos);
	
	if(hasVETO){
	  plot(D_DECAY_BACK_ENERGY_VETO, backEnergy);
	  plot(DD_DECAY_BACK_ENERGY__POSITION_VETO, backEnergy, backPos);
	}else if(!hasVETO){
	  plot(D_DECAY_BACK_ENERGY_NOVETO, backEnergy);
	  plot(DD_DECAY_BACK_ENERGY__POSITION_NOVETO, backEnergy, backPos);
	}
	if(hasNaI){
	  plot(D_DECAY_BACK_ENERGY_NAI, backEnergy);
	}else if(!hasNaI){
	  plot(D_DECAY_BACK_ENERGY_NONAI, backEnergy);
	}
	if(hasPin){
	  plot(D_DECAY_BACK_ENERGY_PIN, backEnergy);
	}else if(!hasPin){
	  plot(D_DECAY_BACK_ENERGY_NOPIN, backEnergy);
	}

	
	

      }
    if (hasFront && hasBack)
      plot(DD_DECAY_POSITION, frontPos, backPos);
  }
  
  
  EndProcess();

    return true;

}

bool Dssd4JAEAProcessor::pickEventType(JAEAEvent& event) {
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

  double cutoffEnergy=6500;
  
  int condition = 0;
  if (event.get_beam()){ 
    condition += 1;
  }
  if (event.get_mwpc() > 0) {
    condition += 2;
  }
  if (event.get_veto()) {
    condition += 4;
  }
  if (condition == 0) {
    double energy = event.get_energy();
      if (energy < cutoffEnergy)
	event.set_type(alpha_jaea);
      else 
	event.set_type(fission_jaea);
    }
  else if (condition == 1) {
    event.set_type(unknown_jaea);
  } 
    else if (condition == 2 || 
             condition == 4 ||
             condition == 6) {
        event.set_type(heavyIon_jaea);
    } 
    else if (condition == 3 && event.get_energy() > recoilEnergyCut_) {
        event.set_type(heavyIon_jaea);
    }
    else if (condition == 5 || condition == 7) {
        event.set_type(lightIon_jaea);
    }
    else
        event.set_type(unknown_jaea);
    

  


    return true;
}
