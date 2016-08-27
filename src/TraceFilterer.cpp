/** \file  TraceFilterer.cpp
 *  \brief Implements the analysis of traces using trapezoidal filters
 *
 *  This trace plots the filtered traces as well as energy and timing info
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>

#include "DammPlotIds.hpp"
#include "RandomPool.hpp"
#include "Trace.hpp"
#include "TraceFilterer.hpp"
#include "Globals.hpp"

using namespace std;
using namespace dammIds::trace;


 //< TO BE USED WITH MAGIC +40 ENERGY SAMPLE LOCATION
// const double TraceFilterer::energyScaleFactor = 2.198;
// const double TraceFilterer::energyScaleFactor = 2.547; //< multiply the energy filter sums by this to gain match to raw spectra

/** 
 *  A do nothing constructor
 */
TraceFilterer::PulseInfo::PulseInfo()
{
    isFound = false;
}

/**
 *  Constructor so we can build the pulseinfo in place
 */
TraceFilterer::PulseInfo::PulseInfo(Trace::size_type theTime, 
                                    double theEnergy) :
                                    time(theTime), energy(theEnergy)
{
    isFound = true;
}

TraceFilterer::TraceFilterer(double energyScaleFactor,
                             short fast_rise, short fast_gap,
                             short fast_threshold,
                             short energy_rise, short energy_gap,
                             short slow_rise, short slow_gap,
                             short slow_threshold) :
    fastParms(fast_gap, fast_rise), 
    energyParms(energy_gap, energy_rise), 
    thirdParms(slow_gap, slow_rise)
{
    //? this uses some legacy values for third parms, are they appropriate
    name = "Filterer";
    useThirdFilter = false;
    energyScaleFactor_ = energyScaleFactor;

    //Trace::size_type rise, gap;

    // scale thresholds by the length of integration (i.e. rise time)
    fastThreshold = fast_threshold * fastParms.GetRiseSamples();
    slowThreshold = slow_threshold * thirdParms.GetRiseSamples();
}


TraceFilterer::~TraceFilterer()
{
    // do nothing
}

bool TraceFilterer::Init(const string &filterFileName /* = filter.txt */)
{
    const int maxTraceLength = 6400;

    fastFilter.reserve(maxTraceLength);
    energyFilter.reserve(maxTraceLength);
    thirdFilter.reserve(maxTraceLength);
    return true;
}


void TraceFilterer::DeclarePlots(void)
{

    const int energyBins = SE;
    const int energyBins2 = SB;
    const int traceBins = dammIds::trace::traceBins;

    using namespace dammIds::trace::tracefilterer;
    /*
     * Declare plots within the trace object
     */
    Trace sample_trace = Trace(); // initialization; by YX
    unsigned short numTraces = Globals::get()->numTraces();

    sample_trace.DeclareHistogram2D(DD_TRACE, traceBins, numTraces,
                                    "traces data TracePlotter"); // DD_TRACE, 7500, all traces here; by YX
    sample_trace.DeclareHistogram2D(DD_FILTER1, traceBins, numTraces,
                                    "fast filter");
    sample_trace.DeclareHistogram2D(DD_FILTER2, traceBins, numTraces,
                                    "energy filter");
    if (useThirdFilter) {
        sample_trace.DeclareHistogram2D(DD_FILTER3, traceBins, numTraces,
                                        "3rd filter");
    }
    sample_trace.DeclareHistogram2D(DD_REJECTED_TRACE, traceBins, numTraces,
                                    "rejected traces");

    sample_trace.DeclareHistogram1D(D_ENERGY1, energyBins, "E1 from trace"); 

    sample_trace.DeclareHistogram2D(DD_ENERGY__BOARD_FILTER,
                                 energyBins2, energyBins2, 
                                "Board raw energy vs filter energy (/10)"); 
    sample_trace.DeclareHistogram1D(D_RATIO_BOARD_FILTER,
                energyBins2, "Ratio raw energy to filter (%)"); 

    /*
    //---------------- by Yongchi Xiao; 04/10/2015 -------------------
    sample_trace.DeclareHistogram2D(DD_MYTEST1, traceBins, numTraces, "my traces");  // 7509
    sample_trace.DeclareHistogram2D(DD_MYTEST2, traceBins, numTraces, "my fast filter");// 7510
    //----------------------------------------------------------------
    */
}

void TraceFilterer::Analyze(Trace &trace,
			    const string &type, const string &subtype)
{
    using namespace dammIds::trace::tracefilterer;

    if (level >= 5) {
        const size_t baselineBins = 30;
        const double deviationCut = fastThreshold / 4. /
                                    fastParms.GetRiseSamples();

        double trailingBaseline  = trace.DoBaseline(
                                trace.size() - baselineBins - 1, baselineBins);

	//--------------------- by Yongchi Xiao --------------------
	
	/* The following line for plotting spec.7500(all traces)
	 * originally here
	 * Change the position of this line on 04/10/2015;  by Yongchi Xiao
	 * so that "ONLY traces w/ energy in a certain range can have their traces plotted."
	 * Modified on: 04/10/2015;
	 */

	//	static int numTracesMine = 0;

	trace.Plot(DD_TRACE,numTracesAnalyzed++); // 7500, by YX

	//-----------------------------------------------------------

	//------------------------- by Yongchi Xiao; 04/10/2015 ---------------------------------- 
	
        // start at sample 5 because first samples are occasionally corrupted
        trace.DoBaseline(5, baselineBins);       


	// --- by Yongchi Xiao; 03/02/2016 --- //
	// comment out the following part to see more hist on DSSD
	// even though they might be noise
	
        if ( trace.GetValue("sigmaBaseline") > deviationCut ||
            abs(trailingBaseline - trace.GetValue("baseline")) < deviationCut)
        {	    
            // perhaps check trailing baseline deviation
            // from a simple linear fit 
            static int rejectedTraces = 0;
            unsigned short numTraces = Globals::get()->numTraces();
            if (rejectedTraces < numTraces)
                trace.Plot(DD_REJECTED_TRACE, rejectedTraces++);
            EndAnalyze(); // update timing
            return;
        }
	
	// --- --

        fastFilter.clear();
        energyFilter.clear();

        // determine trace filters, these are trapezoidal filters characterized
        //   by a risetime and a gaptime and a range of the filter 
        trace.TrapezoidalFilter(fastFilter, fastParms);
        trace.TrapezoidalFilter(energyFilter, energyParms);

        if (useThirdFilter) {
            thirdFilter.clear();
            trace.TrapezoidalFilter(thirdFilter, thirdParms);
        }
        FindPulse(fastFilter.begin(), fastFilter.end()); // by YX; This function is defined below;

        if (pulse.isFound) {

	  trace.SetValue("filterTime", (int)pulse.time);
	  trace.SetValue("filterEnergy", pulse.energy); // pulse.energy -> filterEnergy; by YX;
	  /* pulse.energy? We cannot tell whether it is a signal on the front/back, 
	   *or is it a first/second signal
	   */

	  //---------------------- by Yongchi Xiao; 04/10/2015 --------------------
	  /*
	  if( (pulse.energy > (354.9-5.51) ) && 
	      (pulse.energy < (354.9+5.51) ) &&
	      (pulse.energy > (353.7-6.55) ) &&
	      (pulse.energy < (353.7+6.55) ) ) // looking for 109Te traces;
	  */  /*
	  if(pulse.energy > 600) //
	    {
	      trace.Plot(DD_MYTEST1, numTracesMine); // 7509
	      energyFilter.ScalePlot(DD_MYTEST2, numTracesMine++, energyParms.GetRiseSamples() ); // 7510
	      }  // record traces corresponding to 8764keV signals; */
	  //-----------------------------------------------------------------------
        }

        // now plot some stuff
        fastFilter.ScalePlot(DD_FILTER1, numTracesAnalyzed, 
                    fastParms.GetRiseSamples() );

	// --- for test only --- //
	// by Yongchi; 02/20/2016
	//	cout << "\n" << numTracesAnalyzed << ", " << energyParms.GetRiseSamples()
	//   << endl;

	// --- //
        energyFilter.ScalePlot(DD_FILTER2, numTracesAnalyzed,
			       energyParms.GetRiseSamples() ); /* 7502
								* this is in consistency with my own SlowFilter(); by YX;
								* It is JUST what I need? 
								*/

        if (useThirdFilter) {
            thirdFilter.ScalePlot(DD_FILTER3, numTracesAnalyzed,
                    thirdParms.GetRiseSamples() );
        }
        trace.plot(D_ENERGY1, pulse.energy);

    } // sufficient analysis level

    EndAnalyze(trace);
}

const TraceFilterer::PulseInfo& TraceFilterer::FindPulse(Trace::iterator begin, Trace::iterator end)
{
    // class to see if fast filter is above threshold
    static binder2nd< greater<Trace::value_type> > crossesThreshold
	(greater<Trace::value_type>(), fastThreshold);

    //cout << "fast threshold: " << fastThreshold << endl; // by Yongchi Xiao
    Trace::size_type sample;
    Trace::size_type presample;
    pulse.isFound = false;

    while (begin < end) {
        begin = find_if(begin, end, crossesThreshold);
        if (begin == end) {
            break;
        }
                          
        pulse.time = begin - fastFilter.begin();	
        pulse.isFound = true;
        presample = pulse.time - fastParms.GetRiseSamples();

        // sample the slow filter in the middle of its size
        if (useThirdFilter) {
            sample = pulse.time + (thirdParms.GetSize() - fastParms.GetSize()) / 2;
            if (sample >= thirdFilter.size() ||
            thirdFilter[sample] < slowThreshold) {
                begin++; 
                continue;
            }
        }
        //? some investigation needed here for good resolution
        // add a presample location
        // sample = pulse.time + (energyParms.GetSize() - fastParms.GetSize()) / 2;

	//--------- by Yongchi Xiao; 05/06/2015 --------
        //sample = pulse.time + 20;
	sample = pulse.time + 20; // change the sampling point; by Yongchi Xiao
	//cout << "pulse time: " << pulse.time << endl;
	//cout << "sample:     " << sample << endl;
	//----------------------------------------------
        
        RandomPool* randoms = RandomPool::get();
        if (sample < energyFilter.size()) {
	  pulse.energy = energyFilter[sample] + randoms->Get();	    
            // subtract an energy filter baseline
	  //cout << "---presample--- " << presample << endl; // by Yongchi Xiao; 05/06/2015

	  /*
            if (presample >= 0) {
	      pulse.energy -= energyFilter[presample];
            }
	  */
	  // commented out by Yongchi Xiao; 05/06/2015


            // scale to the integration time
            pulse.energy /= energyParms.GetRiseSamples();
            pulse.energy *= energyScaleFactor_;	    
	    //cout << "pulse energy: " << pulse.energy << endl; // by Yongchi Xiao;05/06/2015

        } else 
            pulse.energy = NAN;

        break;
    }

    return pulse;
}
