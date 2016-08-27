/*! \file Correlator.hpp
 *  \brief Header file for correlation
 *
 *  Class which handles implant/decay correlations
 */

#ifndef __CORRELATOR_PROCESSOR_HPP_
#define __CORRELATOR_PROCESSOR_HPP_

#include <utility>
#include <vector>

#include <cmath>

#include "Plots.hpp"
#include "DammPlotIds.hpp"
#include "Globals.hpp"
#include "param.h"

// forward declarations
class LogicProcessor;
class RawEvent;

struct ImplantData
{
    double time;    ///< time of an implant
    double dtime;   ///< time elapsed since previous implant
    bool implanted; ///< previous implant flag
    double tacValue; ///< a TAC value
    bool flagged; ///< previous event in pixel is flagged as interesting

    ImplantData() {
	Clear();
    }
    void Clear(void) { 
	time = dtime = tacValue = NAN; 
	flagged = implanted = false;
    }
};

struct DecayData
{
    double time;    ///< time of a decay
    double dtime;   ///< time elapsed since previous implant
  
    DecayData() {
	Clear();
    }
    void Clear(void) {
	time = dtime = NAN;
    }
};

struct ListData
{
    double time;
    double energy;
    double offTime;
    
    unsigned long clockCount;
    unsigned char logicBits[4];

    ListData(double t, double e, LogicProcessor *lp = NULL);
};


struct EventInfo
{
    /// types of events passed to the correlator
    enum EEventTypes {IMPLANT_EVENT, ALPHA_EVENT, BETA_EVENT, FISSION_EVENT, 
		      PROTON_EVENT, DECAY_EVENT, PROJECTILE_EVENT, GAMMA_EVENT,
		      UNKNOWN_EVENT};

    EEventTypes type; ///< event type
    double time;     ///< timestamp of event
    double dtime;    ///< time since implant [pixie units]
    double energy;   ///< energy of event
    double energyBox; ///< energy depositied into the box
    double offTime;  ///< length of time beam has been off
    double foilTime; ///< time difference to foil event
    double tof;      ///< time of flight for an implant
    double position; ///< calculated strip position
    short  boxMult;  ///< numebr of box hits
    short  boxMax;   ///< location of maximum energy in box
    short  impMult;  ///< number of implant hits
    short  mcpMult;  ///< number of mcp hits
    short  generation; ///< generation number (0 = implant)
    bool   flagged;  ///< flagged of interest
    bool   hasTof;   ///< has time of flight data
    bool   hasVeto;  ///< veto detector has been hit
    bool   beamOn;   ///< beam is on target
    bool   pileUp;   ///< trace is piled-up

    unsigned long clockCount;
    unsigned char logicBits[dammIds::logic::MAX_LOGIC+1];

    EventInfo();
    EventInfo(double t, double e, LogicProcessor *lp);
};



class CorrelationList : public std::vector<EventInfo>
{
private:
    bool flagged;
public:
    CorrelationList();
    double GetDecayTime(void) const;
    double GetImplantTime(void) const;
    void Flag(void);
    bool IsFlagged(void) const;
    // overide the vector clear function so that the flag is also removed
    void clear(void);
    void PrintDecayList(void) const;
private:
};



/*!
  \brief correlate decays with previous implants
  
  The class controls the correlations of decays with previous implants.  There
  are arrays of size MAX_STRIP %x MAX_STRIP that store implants and decays.
  When an event has been identified as either an implant or decay, its
  information is placed in the appropriate array based on its pixel location.
  If a decay was identified, it is correlated with a previous implant.  The
  correlator checks to make sure that the time between implants is
  sufficiently long and that the correlation time has not been exceeded
  before correlating an implant with a decay.
*/

class Correlator
{
public:

    /// correlator condition based on the given events
    enum EConditions {INVALID_LOCATION     = 4,
		      VALID_IMPLANT        = 12,
		      VALID_DECAY          = 16,
		      BACK_TO_BACK_IMPLANT = 32,
		      DECAY_TOO_LATE       = 48,
		      IMPLANT_TOO_SOON     = 52,
		      UNKNOWN_CONDITION    = 100};
  
  enum EEventType{IMPLANT_EVENT,DECAY_EVENT,UNKNOWN_TYPE};
    Correlator();
    virtual ~Correlator();
  
  void DeclarePlots(void);
  void Init(RawEvent &rawev);
  //from old code
  void Correlate(EventInfo &event, 
		 unsigned int fch, unsigned int bch);
  void CorrelateOld(RawEvent &rawev, EEventType type, unsigned int fch, 
		    unsigned int bch, double time,double energy=0);    
  void CorrelateAll(EventInfo &event); 
  void CorrelateAllX(EventInfo &event, unsigned int bch);
  void CorrelateAllY(EventInfo &event, unsigned int fch);
  void PrintDecayList(unsigned int fch, unsigned int bch) const;
  
  double GetDecayTime(void) const;
  //double GetDecayTime(int fch, int bch) const;
  //double GetImplantTime(void) const;
  //double GetImplantTime(int fch, int bch) const;
  
  // double GetDecayTime(void) const {
  //   return lastDecay->dtime;
  // }
  double GetDecayTime(int fch, int bch) const
  {
    return 0;
    //return decaylistold[fch][bch].GetDecayTime();
  }
  double GetImplantTime(void) const {
      return lastImplant->time;
  }
  double GetImplantTime(int fch, int bch) const {
      return implant[fch][bch].time;
  }
  /*
   void Flag(int fch, int bch) {
      implant[fch][bch].flagged = true;
  }
  */
 
  
  void Flag(int fch, int bch);
  bool IsFlagged(int fch, int bch);
  EConditions GetCondition(void) const {
    return condition;
  }
 
private:
  typedef std::vector<ListData> corrlist_t;
  
  Plots histo;
    virtual void plot(int dammId, double val1, double val2 = -1, double val3 = -1, const char* name="h") {
        histo.Plot(dammId, val1, val2, val3, name);
    }
    virtual void DeclareHistogram1D(int dammId, int xSize, const char* title) {
        histo.DeclareHistogram1D(dammId, xSize, title);
    }
    virtual void DeclareHistogram2D(int dammId, int xSize, int ySize, const char* title) {
        histo.DeclareHistogram2D(dammId, xSize, ySize, title);
    }
        
    static const size_t arraySize = 40; /**< Size of the 2D array to hold the decay lists */

    // in units of pixie clocks
    static const double minImpTime; /**< The minimum amount of time that must
				       pass before an implant will be considered
				     for correlation */
    static const double corrTime;   /**< The maximum amount of time allowed
				       between a decay and its previous implant
				       for a correlation between the two to occur */
    static const double fastTime;   /**< Times shorter than this are output as a fast decay */
    
  EventInfo   *lastImplant;  ///< last implant processed by correlator
  EventInfo   *lastDecay;    ///< last decay procssed by correlator
  
  ImplantData   *lastImplantold;  ///< last implant processed by correlator
  DecayData   *lastDecayold;    ///< last decay procssed by correlator
  
    EConditions condition;     ///< condition for last processed event
  // CorrelationList decaylist[arraySize][arraySize]; ///< list of event data for a particular pixel since implant 
  
  corrlist_t decaylistold[40][40]; ///< list of event data for a particular pixel since implant
  ImplantData implant[40][40];
  DecayData   decay[40][40];
  LogicProcessor *logicProc; ///< a logic processor from the detector driver
  
  
  
};

#endif // __CORRELATOR_PROCESSOR_HPP_
