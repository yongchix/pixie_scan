/*    the detector driver class.

      The main analysis program.  A complete event is create in PixieStd
      passed into this class.  See manual for further details.

      SNL - 7-2-07
      SNL - 7-12-07
            Add root analysis. If the ROOT program has been
            detected on the computer system the and the
            makefile has the useroot flag declared, ROOT
            analysis will be included.
      DTM - Oct. '09
            Significant structural/cosmetic changes. Event processing is
	    now primarily handled by individual event processors which
	    handle their own DetectorDrivers

      SVP - Oct. '10
            Added the VandleProcessor for use with VANDLE.
            Added the PulserProcessor for use with Pulsers.
            Added the WaveformProcessor to determine ps time resolutions.
      KM  - Dec. '12, Jan. '13 
            Huge changes due to switching to XML configuration file.
            See git commits comments.
*/

/*!
  \file DetectorDriver.cpp

  \brief event processing

  In this file are the details for experimental processing of a raw event
  created by ScanList() in PixieStd.cpp.  Event processing includes things
  which do not change from experiment to experiment (such as energy
  calibration and raw parameter plotting) and things that do (differences
  between MTC and RMS experiment, for example).
*/

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>

#include "pugixml.hpp"

#include "DammPlotIds.hpp"
#include "DetectorDriver.hpp"
#include "DetectorLibrary.hpp"
#include "Exceptions.hpp"
#include "RandomPool.hpp"
#include "RawEvent.hpp"
#include "TimingInformation.hpp"
#include "TreeCorrelator.hpp"

#include "BeamLogicProcessor.hpp"
#include "BetaScintProcessor.hpp"
#include "Beta4Hen3Processor.hpp"
#include "DssdProcessor.hpp"
#include "Dssd4SHEProcessor.hpp"
#include "Dssd4JAEAProcessor.hpp"
#include "NaIProcessor.hpp"
#include "PINProcessor.hpp"
#include "Hen3Processor.hpp"
#include "GeProcessor.hpp"
#include "GeCalibProcessor.hpp"
#include "Ge4Hen3Processor.hpp"
#include "ImplantSsdProcessor.hpp"
#include "IonChamberProcessor.hpp"
#include "LiquidScintProcessor.hpp"
#include "McpProcessor.hpp"
#include "MtcProcessor.hpp"
#include "NeutronScintProcessor.hpp"
#include "PositionProcessor.hpp"
#include "PulserProcessor.hpp"
#include "SsdProcessor.hpp"
#include "TraceFilterer.hpp"
#include "TriggerLogicProcessor.hpp"
#include "VandleProcessor.hpp"

#include "CfdAnalyzer.hpp"
#include "DoubleTraceAnalyzer.hpp"
#include "FittingAnalyzer.hpp"
#include "TauAnalyzer.hpp"
#include "TraceAnalyzer.hpp"
#include "TraceExtracter.hpp"
#include "WaveformAnalyzer.hpp"
#include "Notebook.hpp"

#ifdef useroot
#include "RootProcessor.hpp"
#include "ScintRoot.hpp"
#include "VandleRoot.hpp"
#endif

using namespace std;
using namespace dammIds::raw;


/*!
  detector driver constructor

  Creates instances of all event processors
*/

using namespace dammIds::raw;


DetectorDriver* DetectorDriver::instance = NULL;

/** Instance is created upon first call */
DetectorDriver* DetectorDriver::get() {
    if (!instance) {
        instance = new DetectorDriver();
    }
    return instance;
}

DetectorDriver::DetectorDriver() : histo(OFFSET, RANGE, "DetectorDriver") 
{
    Messenger m;
    try {
        m.start("Loading Processors");
        LoadProcessors(m);
    } catch (GeneralException &e) {
        // Any exception in registering plots in Processors 
        // and possible other exceptions in creating Processors
        // will be intercepted here
        m.fail();
        cout << "Exception caught at DetectorDriver::DetectorDriver" << endl;
        cout << "\t" << e.what() << endl;
        exit(EXIT_FAILURE);
    } catch (GeneralWarning &w) {
        cout << "Warning found at DetectorDriver::DetectorDriver" << endl;
        cout << "\t" << w.what() << endl;
    }
    m.done();
}

/*!
  detector driver deconstructor

  frees memory for all event processors
 */
DetectorDriver::~DetectorDriver()
{
    for (vector<EventProcessor *>::iterator it = vecProcess.begin();
	 it != vecProcess.end(); it++) {
	delete *it;
    }

    vecProcess.clear();

    for (vector<TraceAnalyzer *>::iterator it = vecAnalyzer.begin();
	 it != vecAnalyzer.end(); it++) {
	delete *it;
    }

    vecAnalyzer.clear();
}

void DetectorDriver::LoadProcessors(Messenger& m) {
    pugi::xml_document doc;

    pugi::xml_parse_result result = doc.load_file("Config.xml");
    if (!result) {
        stringstream ss;
        ss << "DetectorDriver: error parsing file Config.xml";
        ss << " : " << result.description();
        throw IOException(ss.str());
    }

    /** Force to create Detector Library before loading processors.
     * It is not needed but the output looks nicer ;)
     */
    DetectorLibrary* modChan = DetectorLibrary::get();

    pugi::xml_node driver = doc.child("Configuration").child("DetectorDriver");
    for (pugi::xml_node processor = driver.child("Processor"); processor;
         processor = processor.next_sibling("Processor")) {
        string name = processor.attribute("name").value();

        m.detail("Loading " + name);
        if (name == "BeamLogicProcessor") {
            vecProcess.push_back(new BeamLogicProcessor()); 
        }
        else if (name == "BetaScintProcessor" || name == "Beta4Hen3Processor") {
            double gamma_beta_limit = 
                processor.attribute("gamma_beta_limit").as_double(-1);
            if (gamma_beta_limit == -1) {
                gamma_beta_limit = 200e-9;
                m.warning("Using default gamme_beta_limit = 200e-9", 1);
            }

            double energy_contraction = 
                processor.attribute("energy_contraction").as_double(-1);
            if (energy_contraction == -1) {
                energy_contraction = 1;
                m.warning("Using default energy contraction = 1", 1);
            }

            if (name == "BetaScintProcessor") 
                vecProcess.push_back(
                        new BetaScintProcessor(gamma_beta_limit,
                                               energy_contraction));
            else if (name == "Beta4Hen3Processor")
                vecProcess.push_back(
                        new Beta4Hen3Processor(gamma_beta_limit,
                                               energy_contraction));
        } else if (name == "DssdProcessor") {
            vecProcess.push_back(new DssdProcessor());
        } else if (name == "Dssd4SHEProcessor") {
            double front_back_correlation_time = 
               processor.attribute("front_back_correlation_time").as_double(-1);
            if (front_back_correlation_time == -1) {
                front_back_correlation_time = 300e-9;
                m.warning("Using default front_back_correlation_time = 300e-9",
                          1);
            }
            double delta_energy = 
               processor.attribute("delta_energy").as_double(-1);
            if (delta_energy == -1) {
                delta_energy = 300;
                m.warning("Using default delta_energy = 300", 1);
            }
	    double recoil_energy_cut = 
               processor.attribute("recoil_energy_cut").as_double(-1);
            if (recoil_energy_cut == -1) {
                recoil_energy_cut = 6000;
                m.warning("Using default recoil_energy_cut = 6000", 1);
            }
            double high_energy_cut = 
                processor.attribute("high_energy_cut").as_double(-1);
            if (high_energy_cut == -1) {
                high_energy_cut = 15000.0;
                m.warning("Using default high_energy_cut = 15000", 1);
            }
            double low_energy_cut = 
                processor.attribute("low_energy_cut").as_double(-1);
            if (low_energy_cut == -1) {
                low_energy_cut = 8000.0;
                m.warning("Using default low_energy_cut = 8000", 1);
            }
            double fission_energy_cut = 
                processor.attribute("fission_energy_cut").as_double(-1);
            if (low_energy_cut == -1) {
                low_energy_cut = 50000.0;
                m.warning("Using default fission_energy_cut = 50000", 1);
            }

            stringstream ss;
            unsigned num_front_strips =  
                modChan->GetLocations("dssd_front", "dssd_front").size();
            ss << "Found " << num_front_strips << " front strips (Y)";
            m.detail(ss.str(), 1);
            unsigned num_back_strips =  
                modChan->GetLocations("dssd_back", "dssd_back").size();
            ss.str("");
            ss << "Found " << num_back_strips << " back strips (X)";
            m.detail(ss.str(), 1);

            vecProcess.push_back(new
                    Dssd4SHEProcessor(front_back_correlation_time,
                                      delta_energy,
				      recoil_energy_cut,
                                      high_energy_cut,
                                      low_energy_cut,
                                      fission_energy_cut,
                                      num_back_strips,
                                      num_front_strips));
        }
	else if (name == "Dssd4JAEAProcessor") {
	  double front_back_correlation_time = 
	    processor.attribute("front_back_correlation_time").as_double(-1);
	  if (front_back_correlation_time == -1) {
	    front_back_correlation_time = 300e-9;
	    m.warning("Using default front_back_correlation_time = 300e-9",
		      1);
	  }
	  double delta_energy = 
	    processor.attribute("delta_energy").as_double(-1);
	  if (delta_energy == -1) {
	    delta_energy = 300;
	    m.warning("Using default delta_energy = 300", 1);
	  }
	  double recoil_energy_cut = 
	      processor.attribute("recoil_energy_cut").as_double(-1);
            if (recoil_energy_cut == -1) {
	      recoil_energy_cut = 6000;
	      m.warning("Using default recoil_energy_cut = 6000", 1);
            }
            double high_energy_cut = 
                processor.attribute("high_energy_cut").as_double(-1);
            if (high_energy_cut == -1) {
                high_energy_cut = 15000.0;
                m.warning("Using default high_energy_cut = 15000", 1);
            }
            double low_energy_cut = 
                processor.attribute("low_energy_cut").as_double(-1);
            if (low_energy_cut == -1) {
                low_energy_cut = 8000.0;
                m.warning("Using default low_energy_cut = 8000", 1);
            }
            double fission_energy_cut = 
                processor.attribute("fission_energy_cut").as_double(-1);
            if (fission_energy_cut == -1) {
                fission_energy_cut = 50000.0;
                m.warning("Using default fission_energy_cut = 50000", 1);
            }
			stringstream ss;
            unsigned num_front_strips =  
                modChan->GetLocations("dssd_front_jaea", "dssd_front_jaea").size();
            ss << "Found " << num_front_strips << " front strips (Y)";
            m.detail(ss.str(), 1);
            unsigned num_back_strips =  
                modChan->GetLocations("dssd_back_jaea", "dssd_back_jaea").size();
            ss.str("");
            ss << "Found " << num_back_strips << " back strips (X)";
            m.detail(ss.str(), 1);
			double correlationMatrixWin = 
				processor.attribute("correlationMatrixWin").as_double(-1);
			if (correlationMatrixWin == -1) {
				correlationMatrixWin = 30e-3;
				m.warning("Using default correlationMatrixWin = 30e-3", 1);
			}
			double betaWin = 
				processor.attribute("betaWin").as_double(-1);
			if (betaWin == -1) {
				betaWin = 30e-6;
				m.warning("Using default betaWin = 30e-6", 1); 
			}
            vecProcess.push_back(new
                    Dssd4JAEAProcessor(front_back_correlation_time,
									   delta_energy,
									   recoil_energy_cut,
									   high_energy_cut,
									   low_energy_cut,
									   fission_energy_cut,
									   num_back_strips,
									   num_front_strips, 
									   correlationMatrixWin, 
									   betaWin));
        }
	else if (name == "GeProcessor" || name == "Ge4Hen3Processor") {
            double gamma_threshold = 
                processor.attribute("gamma_threshold").as_double(-1);
            if (gamma_threshold == -1) {
                gamma_threshold = 1.0;
                m.warning("Using default gamma_threshold = 1.0", 1);
            }
            double low_ratio =
                processor.attribute("low_ratio").as_double(-1);
            if (low_ratio == -1) {
                low_ratio = 1.0;
                m.warning("Using default low_ratio = 1.0", 1);
            }
            double high_ratio =
                processor.attribute("high_ratio").as_double(3);
            if (high_ratio == -1) {
                high_ratio = 3.0;
                m.warning("Using default high_ratio = 3.0", 1);
            }
            double sub_event = 
                processor.attribute("sub_event").as_double(-1);
            if (sub_event == -1) {
                sub_event = 100e-9;
                m.warning("Using default sub_event = 100e-9", 1);
            }
            double gamma_beta_limit = 
                processor.attribute("gamma_beta_limit").as_double(-1);
            if (gamma_beta_limit == -1) {
                gamma_beta_limit = 200e-9;
                m.warning("Using default gamme_beta_limit = 200e-9", 1);
            }
            double gamma_gamma_limit = 
                processor.attribute("gamma_gamma_limit").as_double(-1);
            if (gamma_gamma_limit == -1) {
                gamma_gamma_limit = 200e-9;
                m.warning("Using default gamma_gamma_limit = 200e-9", 1);
            }
            double cycle_gate1_min =
                processor.attribute("cycle_gate1_min").as_double(-1);
            if (cycle_gate1_min == -1) {
                cycle_gate1_min = 0.0;
                m.warning("Using default cycle_gate1_min = 0.0", 1);
            }
            double cycle_gate1_max =
                processor.attribute("cycle_gate1_max").as_double(-1);
            if (cycle_gate1_max == -1) {
                cycle_gate1_max = 0.0;
                m.warning("Using default cycle_gate1_max = 0.0", 1);
            }
            double cycle_gate2_min =
                processor.attribute("cycle_gate2_min").as_double(-1);
            if (cycle_gate2_min == -1) {
                cycle_gate2_min = 0.0;
                m.warning("Using default cycle_gate2_min = 0.0", 1);
            }
            double cycle_gate2_max =
                processor.attribute("cycle_gate2_max").as_double(-1);
            if (cycle_gate2_max == -1) {
                cycle_gate2_max = 0.0;
                m.warning("Using default cycle_gate2_max = 0.0", 1);
            }
            if (name == "GeProcessor") {
                vecProcess.push_back(new GeProcessor(gamma_threshold,
                            low_ratio, high_ratio, sub_event,
                            gamma_beta_limit, gamma_gamma_limit,
                            cycle_gate1_min, cycle_gate1_max,
                            cycle_gate2_min, cycle_gate2_max)); 
            } else if (name == "Ge4Hen3Processor") {
	      vecProcess.push_back(new Ge4Hen3Processor(gamma_threshold,
							  low_ratio, high_ratio, sub_event,
							  gamma_beta_limit, gamma_gamma_limit, 
							  cycle_gate1_min, cycle_gate1_max,
							  cycle_gate2_min, cycle_gate2_max)); 
		
            }
        } else if (name == "GeCalibProcessor") {
            double gamma_threshold = 
	      processor.attribute("gamma_threshold").as_double(1);
            double low_ratio =
	      processor.attribute("low_ratio").as_double(1);
            double high_ratio =
	      processor.attribute("high_ratio").as_double(3);
            vecProcess.push_back(new GeCalibProcessor(gamma_threshold,
						      low_ratio, high_ratio)); 
        } else if (name == "Hen3Processor") {
	  vecProcess.push_back(new Hen3Processor()); 
        } else if (name == "ImplantSsdProcessor") {
	  vecProcess.push_back(new ImplantSsdProcessor()); 
        } else if (name == "IonChamberProcessor") {
	  vecProcess.push_back(new IonChamberProcessor()); 
        } else if (name == "LiquidScintProcessor") {
	  vecProcess.push_back(new LiquidScintProcessor());
        } else if (name == "LogicProcessor") {
	  vecProcess.push_back(new LogicProcessor()); 
        } else if (name == "McpProcessor") {
	  vecProcess.push_back(new McpProcessor()); 
	} else if (name == "NaIProcessor") {
	  vecProcess.push_back(new NaIProcessor());
	} else if (name == "PINProcessor"){
	  vecProcess.push_back(new PINProcessor());
	} else if (name == "MtcProcessor") {
	  /** Default value for as_bool() is false */
	  bool double_stop = 
	    processor.attribute("double_stop").as_bool();
	  bool double_start = 
	    processor.attribute("double_start").as_bool();
	  vecProcess.push_back(new MtcProcessor(double_stop, double_start));
        } else if (name == "NeutronScintProcessor") {
	  vecProcess.push_back(new NeutronScintProcessor());
	} else if (name == "PulserProcessor") {
	  vecProcess.push_back(new PulserProcessor());
	} else if (name == "SsdProcessor") {
	  vecProcess.push_back(new SsdProcessor());
	} else if (name == "TriggerLogicProcessor") {
	  vecProcess.push_back(new TriggerLogicProcessor());
	} else if (name == "VandleProcessor") {
	  vecProcess.push_back(new VandleProcessor());
	}
#ifdef useroot
        else if (name == "RootProcessor") {
            vecProcess.push_back(new RootProcessor("tree.root", "tree"));
        } else if (name == "ScintROOT") {
            vecProcess.push_back(new ScintROOT());
        } else if (name == "VandleROOT") {
            vecProcess.push_back(new VandleROOT());
        }
#endif
        else {
            stringstream ss;
            ss << "DetectorDriver: unknown processor type" << name;
            throw GeneralException(ss.str());
        }
        stringstream ss;
        for (pugi::xml_attribute_iterator ait = processor.attributes_begin();
             ait != processor.attributes_end(); ++ait) {
            ss.str("");
            ss << ait->name();
            if (ss.str().compare("name") != 0) {
                ss << " = " << ait->value();
                m.detail(ss.str(), 1);
            }
        }
}

    for (pugi::xml_node analyzer = driver.child("Analyzer"); analyzer;
         analyzer = analyzer.next_sibling("Analyzer")) {
        string name = analyzer.attribute("name").value();

        m.detail("Loading " + name);
        if (name == "TraceFilterer" || 
            name == "DoubleTraceAnalyzer") {

            double gain_match = analyzer.attribute("gain_match").as_double(-1);
            if (gain_match == -1) {
                gain_match = 1.0;
                m.warning("Using gain_match = 1.0", 1);
            }
            int fast_rise = analyzer.attribute("fast_rise").as_int(-1);
            if (fast_rise == -1) {
                fast_rise = 10;
                m.warning("Using fast_rise = 10", 1);
            }
            int fast_gap = analyzer.attribute("fast_gap").as_int(-1);
            if (fast_gap == -1) {
                fast_gap = 10;
                m.warning("Using fast_gap = 10", 1);
            }
            int fast_threshold = 
                analyzer.attribute("fast_threshold").as_int(-1);
            if (fast_threshold == -1) {
                fast_threshold = 50;
                m.warning("Using fast_threshold = 50", 1);
            }
            int energy_rise = analyzer.attribute("energy_rise").as_int(-1);
            if (energy_rise == -1) {
                energy_rise = 50;
                m.warning("Using energy_rise = 50", 1);
            }
            int energy_gap = analyzer.attribute("energy_gap").as_int(-1);
            if (energy_gap == -1) {
                energy_gap = 50;
                m.warning("Using energy_gap = 50", 1);
            }
            int slow_rise = analyzer.attribute("slow_rise").as_int(-1);
            if (slow_rise == -1) {
                slow_rise = 20;
                m.warning("Using slow_rise = 20", 1);
            }
            int slow_gap = analyzer.attribute("slow_gap").as_int(-1);
            if (slow_gap == -1) {
                slow_gap = 20;
                m.warning("Using slow_gap = 20", 1);
            }
            int slow_threshold = 
                analyzer.attribute("slow_threshold").as_int(-1);
            if (slow_threshold == -1) {
                slow_threshold = 10;
                m.warning("Using slow_threshold = 10", 1);
            }

            if (name == "TraceFilterer") 
                vecAnalyzer.push_back(new TraceFilterer(
                            gain_match,
                            fast_rise, fast_gap, fast_threshold, 
                            energy_rise, energy_gap,
                            slow_rise, slow_gap, slow_threshold));
            else if (name == "DoubleTraceAnalyzer")
                vecAnalyzer.push_back(new DoubleTraceAnalyzer(
                            gain_match,
                            fast_rise, fast_gap, fast_threshold, 
                            energy_rise, energy_gap,
                            slow_rise, slow_gap, slow_threshold));

        } else if (name == "TauAnalyzer") {
            vecAnalyzer.push_back(new TauAnalyzer());
        } else if (name == "TraceExtracter") {
            string type = analyzer.attribute("type").as_string();
            string subtype = analyzer.attribute("subtype").as_string();

            vecAnalyzer.push_back(new TraceExtracter(type, subtype));
        }
        else if (name == "WaveformAnalyzer") {
            vecAnalyzer.push_back(new WaveformAnalyzer());
        }
        else if (name == "FittingAnalyzer") {
            vecAnalyzer.push_back(new FittingAnalyzer());
        }
        else if (name == "CfdAnalyzer") {
            vecAnalyzer.push_back(new CfdAnalyzer());
        }
        else {
            stringstream ss;
            ss << "DetectorDriver: unknown analyzer type" << name;
            throw GeneralException(ss.str());
        }
        for (pugi::xml_attribute_iterator ait = analyzer.attributes_begin();
             ait != analyzer.attributes_end(); ++ait) {
            stringstream ss;
            ss << ait->name();
            if (ss.str().compare("name") != 0) {
                ss << " = " << ait->value();
                m.detail(ss.str(), 1);
            }
        }
    }
}

/*!
  Called from PixieStd.cpp during initialization.
  The calibration file cal.txt is read using the function ReadCal() and 
  checked to make sure that all channels have a calibration.
*/

int DetectorDriver::Init(RawEvent& rawev)
{
    // initialize the trace analysis routine
    for (vector<TraceAnalyzer *>::iterator it = vecAnalyzer.begin();
	 it != vecAnalyzer.end(); it++) {
	(*it)->Init();
	(*it)->SetLevel(20); //! Plot traces
    }

    // initialize processors in the event processing vector
    for (vector<EventProcessor *>::iterator it = vecProcess.begin();
         it != vecProcess.end(); it++) {
        (*it)->Init(rawev);	
    }

    try {
        ReadCalXml();
        ReadWalkXml();
    } catch (GeneralException &e) {
        // Any exception in reading calibration and walk correction
        // will be intercepted here
        cout << endl;
        cout << "Exception caught at DetectorDriver::Init" << endl;
        cout << "\t" << e.what() << endl;
        Messenger m;
        m.fail();
        exit(EXIT_FAILURE);
    } catch (GeneralWarning &w) {
        cout << "Warning caught at DetectorDriver::Init" << endl;
        cout << "\t" << w.what() << endl;
    }

    TimingInformation readFiles;
    readFiles.ReadTimingConstants();
    readFiles.ReadTimingCalibration();

    rawev.GetCorrelator().Init(rawev);

    return 0;
}

/*!
  \brief controls event processing

  The ProcessEvent() function is called from ScanList() in PixieStd.cpp
  after an event has been constructed. This function is passed the mode
  the analysis is currently in (the options are either "scan" or
  "standaloneroot").  The function checks the thresholds for the individual
  channels in the event and calibrates their energies. 
  The raw and calibrated energies are plotted if the appropriate DAMM spectra
  have been created.  Then experiment specific processing is performed.  
  Currently, both RMS and MTC processing is available.  After all processing
  has occured, appropriate plotting routines are called.
*/
int DetectorDriver::ProcessEvent(RawEvent& rawev){   
    /*
      Begin the event processing looping over all the channels
      that fired in this particular event.
    */
    plot(dammIds::raw::D_NUMBER_OF_EVENTS, dammIds::GENERIC_CHANNEL);
    
    try {
        for (vector<ChanEvent*>::const_iterator it = rawev.GetEventList().begin();
            it != rawev.GetEventList().end(); ++it) {
            string place = (*it)->GetChanID().GetPlaceName();

            // skip empty channel
            if (place == "__-1")
                continue;

            // check threshold and calibrate
            PlotRaw((*it));
            ThreshAndCal((*it), rawev);
            PlotCal((*it));

            // Do not activate places if saturated or pileup
            if ( (*it)->IsSaturated() || (*it)->IsPileup() )
                continue;

            double time = (*it)->GetTime();
            double energy = (*it)->GetCalEnergy();
            int location = (*it)->GetChanID().GetLocation();
	    
  	    
          	DetectorDriver* driver = DetectorDriver::get();
        	time_t theTime;
     	    theTime = driver -> GetWallTime(time);
            if (abs(theTime-1412393220)<=100 )  {
	             stringstream ss;
                 ss << theTime << " " << time << " " << location <<  " " << energy << endl;
                  Messenger m;
                  m.run_message(ss.str());
        Notebook::get()->report(ss.str());
	}
            EventData data(time, energy, location);
            TreeCorrelator::get()->place(place)->activate(data);
        } 
    
        // have each processor in the event processing vector handle the event
        /* First round is preprocessing, where process result must be guaranteed
        * to not to be dependent on results of other Processors. */
        for (vector<EventProcessor*>::iterator iProc = vecProcess.begin();
        iProc != vecProcess.end(); iProc++) {
            if ( (*iProc)->HasEvent() ) {
                (*iProc)->PreProcess(rawev);
            }
        }
        /* In the second round the Process is called, which may depend on other
        * Processors. */
        for (vector<EventProcessor *>::iterator iProc = vecProcess.begin();
        iProc != vecProcess.end(); iProc++) {
            if ( (*iProc)->HasEvent() ) {
                (*iProc)->Process(rawev);
            }
        }
    } catch (GeneralException &e) {
        // Any exception in activation of basic places, PreProcess and Process
        // will be intercepted here
        cout << "Exception caught at DetectorDriver::ProcessEvent" << endl;
        cout << "\t" << e.what() << endl;
        exit(EXIT_FAILURE);
    } catch (GeneralWarning &w) {
      //cout << "Warning caught at DetectorDriver::ProcessEvent" << endl;
      //cout << "\t" << w.what() << endl;
    }

    return 0;   
}

void DetectorDriver::DeclarePlots()
{
    try {
        DetectorLibrary* modChan = DetectorLibrary::get();

        // Declare plots for each channel
	
        DeclareHistogram1D(D_HIT_SPECTRUM, S7, "channel hit spectrum");
        DeclareHistogram1D(D_SUBEVENT_GAP, SE,
                           "time btwn chan-in event,10ns bin");
        DeclareHistogram1D(D_EVENT_LENGTH, SE,
                           "time length of event, 10 ns bin");
        DeclareHistogram1D(D_EVENT_GAP, SE, "time between events, 10 ns bin");
        DeclareHistogram1D(D_EVENT_MULTIPLICITY, S7,
                           "number of channels in event");
        DeclareHistogram1D(D_BUFFER_END_TIME, SE, "length of buffer, 1 ms bin");
        DeclareHistogram2D(DD_RUNTIME_SEC, SE, S6, "run time - s");
        DeclareHistogram2D(DD_DEAD_TIME_CUMUL, SE, S6, "dead time - cumul");
        DeclareHistogram2D(DD_BUFFER_START_TIME, SE, S6, "dead time - 0.1%");
        DeclareHistogram2D(DD_RUNTIME_MSEC, SE, S7, "run time - ms");
        DeclareHistogram1D(D_NUMBER_OF_EVENTS, S4, "event counter");
        DeclareHistogram1D(D_HAS_TRACE, S8, "channels with traces");
	
        DetectorLibrary::size_type maxChan = modChan->size();

        for (DetectorLibrary::size_type i = 0; i < maxChan; i++) {	 
            if (!modChan->HasValue(i)) {
                continue;
            }
            stringstream idstr; 
            
            const Identifier &id = modChan->at(i);

            idstr << "M" << modChan->ModuleFromIndex(i)
                  << " C" << modChan->ChannelFromIndex(i)
                  << " - " << id.GetType()
                  << ":" << id.GetSubtype()
                  << " L" << id.GetLocation();
            DeclareHistogram1D(D_RAW_ENERGY + i, SE,
                               ("RawE " + idstr.str()).c_str() );
            DeclareHistogram1D(D_FILTER_ENERGY  + i, SE,
                               ("FilterE " + idstr.str()).c_str() );
            DeclareHistogram1D(D_SCALAR  + i, SE,
                               ("Scalar " + idstr.str()).c_str() );
            if (Globals::get()->revision() == "A")
                DeclareHistogram1D(D_TIME + i, SE,
                                ("Time " + idstr.str()).c_str() ); 
            DeclareHistogram1D(D_CAL_ENERGY + i, SE,
                               ("CalE " + idstr.str()).c_str() ); 
        }
	/*for (int j=0; j < 12; j++) {
             DeclareHistogram2D(DD_RAW_V_CAL+j,SC,SE,"Raw vs. Calibrated Energy");
	}*/
        // Now declare histograms present in all used analyzers and
        // processors
        for (vector<TraceAnalyzer *>::const_iterator it = vecAnalyzer.begin();
             it != vecAnalyzer.end(); it++) {
            (*it)->DeclarePlots();
        }

        for (vector<EventProcessor *>::const_iterator it = vecProcess.begin();
             it != vecProcess.end(); it++) {
            (*it)->DeclarePlots();
        } 
        
    } catch (exception &e) {
        // Any exception in histogram declaration will be intercepted here
        cout << "Exception caught at DetectorDriver::DeclarePlots" << endl;
        cout << "\t" << e.what() << endl;
        exit(EXIT_FAILURE);
    }
}

// sanity check for all our expectations
void DetectorDriver::SanityCheck(void) const
{
    /** Use Exceptions to throw an exception here if sanity check was 
     * not succesful */
}

/*!
  \brief check threshold and calibrate each channel.

  Check the thresholds and calibrate the energy for each channel using the
  calibrations contained in the calibration vector filled during ReadCal()
*/

int DetectorDriver::ThreshAndCal(ChanEvent *chan, RawEvent& rawev)
{   
    // retrieve information about the channel
    Identifier chanId = chan->GetChanID();
    int id            = chan->GetID();
    string type       = chanId.GetType();
    string subtype    = chanId.GetSubtype();
    bool hasStartTag  = chanId.HasTag("start");
    Trace &trace      = chan->GetTrace();

    RandomPool* randoms = RandomPool::get();

    double energy = 0.0;

    if (type == "ignore" || type == "") {
        return 0;
    }
    /*
      If the channel has a trace get it, analyze it and set the energy.
    */
    if ( !trace.empty() ) {
        plot(D_HAS_TRACE, id);

        for (vector<TraceAnalyzer *>::iterator it = vecAnalyzer.begin();
            it != vecAnalyzer.end(); it++) {	
                (*it)->Analyze(trace, type, subtype);
        }

        if (trace.HasValue("filterEnergy") ) {     
            if (trace.GetValue("filterEnergy") > 0) {
                energy = trace.GetValue("filterEnergy");
                plot(D_FILTER_ENERGY + id, energy);

                /** These plots are used to determine (or check) the
                 * gain_match parameter to match the filter 
                 * and onboard amplitudes
                 */
                using namespace dammIds::trace::tracefilterer;
                //double board_energy = chan->GetEnergy();
                //trace.plot(DD_ENERGY__BOARD_FILTER, 
                //            board_energy / 10.0, energy / 10.0);
                //trace.plot(D_RATIO_BOARD_FILTER,
                //            board_energy / energy * 100.0);

                trace.SetValue("filterEnergyCal",
                    cali.GetCalEnergy(chanId, trace.GetValue("filterEnergy")));
            } else {
                energy = 2;
            }

            /** Calibrate pulses numbered 2 and forth,
             * add filterEnergyXCal to the trace */
            int pulses = trace.GetValue("numPulses");
            for (int i = 1; i < pulses; ++i) {
                stringstream energyName;
                energyName << "filterEnergy" << i + 1;
                stringstream energyCalName;
                energyCalName << "filterEnergy" << i + 1 << "Cal";
                trace.SetValue(energyCalName.str(),
                    cali.GetCalEnergy(chanId, 
                                      trace.GetValue(energyName.str())));
            }
        }

        if (trace.HasValue("calcEnergy") ) {	    
            energy = trace.GetValue("calcEnergy");
            chan->SetEnergy(energy);
        } else if (!trace.HasValue("filterEnergy")) {
            energy = chan->GetEnergy() + randoms->Get();
            energy /= ChanEvent::pixieEnergyContraction;
        }

        if (trace.HasValue("phase") ) {
            double phase = trace.GetValue("phase");
            chan->SetHighResTime( phase * Globals::get()->adcClockInSeconds() + 
                                  chan->GetTrigTime() *
                                  Globals::get()->filterClockInSeconds());
        }

    } else {
        // otherwise, use the Pixie on-board calculated energy
        // add a random number to convert an integer value to a 
        //   uniformly distributed floating point

        energy = chan->GetEnergy() + randoms->Get();
	//cout << "raw --------------------------" << endl
	// << energy << endl;
        energy /= ChanEvent::pixieEnergyContraction; // energy is 4 times smaller now; by YX
	//cout << energy << endl
	//   << "--------------------------" << endl;
    }

    /** Calibrate energy and apply the walk correction. */
    double time = chan->GetTime();
    double walk_correction = walk.GetCorrection(chanId, energy);

    chan->SetCalEnergy(cali.GetCalEnergy(chanId, energy));
    chan->SetCorrectedTime(time - walk_correction);

    /*
      update the detector summary
    */    
    rawev.GetSummary(type)->AddEvent(chan);
    DetectorSummary *summary;
    
    summary = rawev.GetSummary(type + ':' + subtype, false);
    if (summary != NULL)
        summary->AddEvent(chan);

    if(hasStartTag) {
        summary = 
            rawev.GetSummary(type + ':' + subtype + ':' + "start", false);
        if (summary != NULL)
            summary->AddEvent(chan);
    }
    
    return 1;
}

/*!
  Plot the raw energies of each channel into the damm spectrum number assigned
  to it in the map file with an offset as defined in DammPlotIds.hpp
*/
int DetectorDriver::PlotRaw(const ChanEvent *chan)
{
    int id = chan->GetID();
    float energy = chan->GetEnergy() / ChanEvent::pixieEnergyContraction;

    plot(D_RAW_ENERGY + id, energy);

    return 0;
}

/*!
  Plot the calibrated energies of each channel into the damm spectrum number
  assigned to it in the map file with an offset as defined in DammPlotIds.hpp
*/
int DetectorDriver::PlotCal(const ChanEvent *chan)
{
    int id = chan->GetID();
    int Offset=5;
    double energy = chan->GetEnergy() / ChanEvent::pixieEnergyContraction;
    double calEnergy = chan->GetCalEnergy();
    double time = chan->GetTime();
    static double plottime=0; 
    /*if (id <190) {
	stringstream ss;
            //    ss << time << ' ' << time-plottime << ' ' << plottime <<' ';
		ss << energy << ' ' << calEnergy << ' ' << id;
                Messenger m;
                m.run_message(ss.str());	
    }*/
    plot(D_CAL_ENERGY + id, calEnergy);
    /* Check Cal */
    /*if (id < 48 ) {
	Offset = 0;
    } else if (id < 176 && abs(calEnergy-1.77*energy+280)>=360) {
	Offset = 1;
	plottime= time;
	
    } else if (id < 176 ) {
	Offset = 2;
    } else if (id < 190 ) {
	Offset = 3;
    } else if (id < 195 ) {
	Offset = 4;
    } else {
	Offset = 5;
    }
    if (abs(time-plottime)<=1){
	Offset += 6;
    }
    plot(DD_RAW_V_CAL+Offset,energy,calEnergy);*/
    return 0;
}

vector<EventProcessor *> DetectorDriver::GetProcessors(const string& type) const
{
  vector<EventProcessor *> retVec;

  for (vector<EventProcessor *>::const_iterator it = vecProcess.begin();
       it != vecProcess.end(); it++) {
    if ( (*it)->GetTypes().count(type) > 0 )
      retVec.push_back(*it);
  }

  return retVec;
}

EventProcessor* DetectorDriver::GetProcessor(const std::string &name) const {
	for(vector<EventProcessor*> ::const_iterator it = vecProcess.begin();
		it != vecProcess.end();
		it++) {
		if( (*it)->GetName() == name )
			return (*it);
	}
	return (NULL);
}

void DetectorDriver::ReadCalXml() {
    pugi::xml_document doc;

    pugi::xml_parse_result result = doc.load_file("Config.xml");
    if (!result) {
        stringstream ss;
        ss << "DetectorDriver: error parsing file Config.xml";
        ss << " : " << result.description();
        throw GeneralException(ss.str());
    }

    Messenger m;
    m.start("Loading Calibration");

    pugi::xml_node map = doc.child("Configuration").child("Map");

    /* Note that before this reading in of the xml file, it was already
     * processed for the purpose of creating the channels map.
     * Some sanity checks (module and channel number) were done there
     * so they are not repeated here/
     */
    bool verbose = map.attribute("verbose_calibration").as_bool();
    for (pugi::xml_node module = map.child("Module"); module;
         module = module.next_sibling("Module")) {
        int module_number = module.attribute("number").as_int(-1);
        for (pugi::xml_node channel = module.child("Channel"); channel;
             channel = channel.next_sibling("Channel")) {
            int ch_number = channel.attribute("number").as_int(-1);
            Identifier chanID = DetectorLibrary::get()->at(module_number,
                                                           ch_number);
            bool calibrated = false;
            for (pugi::xml_node cal = channel.child("Calibration");
                cal; cal = cal.next_sibling("Calibration")) {
                string model = cal.attribute("model").as_string("None");
                double min = cal.attribute("min").as_double(0);
                double max = 
                  cal.attribute("max").as_double(numeric_limits<double>::max());

                stringstream pars(cal.text().as_string());
                vector<double> parameters;
                while (true) {
                    double p;
                    pars >> p;
                    if (pars) 
                        parameters.push_back(p);
                    else
                        break;
                }
                if (verbose) {
                    stringstream ss;
                    ss << "Module " << module_number << ", channel " 
                       << ch_number << ": ";
                    ss << " model-" << model;
                    for (vector<double>::iterator it = parameters.begin();
                         it != parameters.end(); ++it)
                        ss << " " << (*it);
                    m.detail(ss.str(), 1);
                }
                cali.AddChannel(chanID, model, min, max, parameters);
                calibrated = true;
            }
            if (!calibrated && verbose) {
                stringstream ss;
                ss << "Module " << module_number << ", channel " 
                   << ch_number << ": ";
                ss << " non-calibrated";
                m.detail(ss.str(), 1);
            }
        }
    }
    m.done();
}

void DetectorDriver::ReadWalkXml() {
    pugi::xml_document doc;

    pugi::xml_parse_result result = doc.load_file("Config.xml");
    if (!result) {
        stringstream ss;
        ss << "DetectorDriver: error parsing file Config.xml";
        ss << " : " << result.description();
        throw GeneralException(ss.str());
    }

    Messenger m;
    m.start("Loading Walk Corrections");

    pugi::xml_node map = doc.child("Configuration").child("Map");
    /** See comment in the similiar place at ReadCalXml() */
    bool verbose = map.attribute("verbose_walk").as_bool();
    for (pugi::xml_node module = map.child("Module"); module;
         module = module.next_sibling("Module")) {
        int module_number = module.attribute("number").as_int(-1);
        for (pugi::xml_node channel = module.child("Channel"); channel;
             channel = channel.next_sibling("Channel")) {
            int ch_number = channel.attribute("number").as_int(-1);
            Identifier chanID = DetectorLibrary::get()->at(module_number,
                                                           ch_number);
            bool corrected = false;
            for (pugi::xml_node walkcorr = channel.child("WalkCorrection");
                walkcorr; walkcorr = walkcorr.next_sibling("WalkCorrection")) {
                string model = walkcorr.attribute("model").as_string("None");
                double min = walkcorr.attribute("min").as_double(0);
                double max = 
                  walkcorr.attribute("max").as_double(
                                              numeric_limits<double>::max());

                stringstream pars(walkcorr.text().as_string());
                vector<double> parameters;
                while (true) {
                    double p;
                    pars >> p;
                    if (pars) 
                        parameters.push_back(p);
                    else
                        break;
                }
                if (verbose) {
                    stringstream ss;
                    ss << "Module " << module_number 
                       << ", channel " << ch_number << ": ";
                    ss << " model: " << model;
                    for (vector<double>::iterator it = parameters.begin();
                         it != parameters.end(); ++it)
                        ss << " " << (*it);
                    m.detail(ss.str(), 1);
                }
                walk.AddChannel(chanID, model, min, max, parameters);
                corrected = true;
            }
            if (!corrected && verbose) {
                stringstream ss;
                ss << "Module " << module_number << ", channel " 
                << ch_number << ": ";
                ss << " not corrected for walk";
                m.detail(ss.str(), 1);
            }
        }
    }
    m.done();
}


extern "C" void detectorend_()
{
  //cout << "ending, no rootfile " << endl;       
}
