#ifndef __DAMM_PLOTIDS_HPP_
#define __DAMM_PLOTIDS_HPP_ 1

/**
 * Histogram sizes consts
 */
const int S1 = 2, S2 = 4, S3 = 8, S4 = 16, S5 = 32, S6 = 64, S7 = 128,
    S8 = 256, S9 = 512, SA = 1024, SB = 2048, SC = 4096,
    SD = 8192, SE = 16384, SF = 32768;

namespace dammIds {
    const int GENERIC_CHANNEL = 10;
    
    // in Correlator.cpp
    namespace correlator {
        const int OFFSET = 600;
        const int RANGE = 10;
    } 
  
  // in DssdProcessor.cpp
  namespace dssd {
    const int OFFSET = 800;
    const int RANGE = 100;
  } 
  // in DssdProcessor.cpp
  namespace dssd4she {
    const int OFFSET = 1800;
    const int RANGE = 100;
  }
  // --------- for Tokai experiment--------------
  namespace nai {
    const int OFFSET = 930; // by Yongchi Xiao; 06/10/2015
    const int RANGE  =  20;
  }
  namespace pin {
    const int OFFSET = 950;
    const int RANGE  =  30; // by Yongchi Xiao; 05/11/2016
  }
    namespace dssd4jaea {
      const int OFFSET = 700;
      const int RANGE = 100;
      // by Yongchi Xiao: put 11 spectrum of NaI into this region
      /* occupied by DSSD4JAEA: 
	 700 - 710;
	 714 - 720;
	 725, 726
	 741 - 744; 
	 750 - 757;
	 759 - 796;
       */

      // choose: 
      // START FROM:  721;
      // END AT: 736;
    }
  namespace mcp {	
    const int OFFSET = 920;
    const int RANGE = 10;
  }
  //----------------------------------------------
  namespace beta_scint {
    const int OFFSET = 2050;
    const int RANGE = 50;
  } 
  
  namespace neutron_scint {
    const int OFFSET = 2100;
    const int RANGE = 50;
  } 
  
  namespace liquid_scint {
    const int OFFSET = 2150;
        const int RANGE = 50;
    } 

    namespace hen3 {
        const int OFFSET = 2200;
        const int RANGE = 50;
    }

    namespace ge {
        const int OFFSET = 2500;
        const int RANGE = 500;
    } 

    namespace logic {
        const int OFFSET = 1100;
        const int RANGE = 100;
        const int MAX_LOGIC = 10; /*< maximum number of logic signals */
    }
    namespace triggerlogic {
        const int OFFSET = 1200;
        const int RANGE = 100;
    }

    namespace vandle{ //The RANGE can be greatly reduced. -SVP
        const int OFFSET = 1300;
        const int RANGE = 200;
    }

    //in PulserProcessor.cpp 
    namespace pulser{ 
        const int OFFSET = 1500;
        const int RANGE = 20;
    } 
   
    namespace waveformanalyzer{ 
        const int OFFSET = 1520;
        const int RANGE = 20;
    }

    // in SsdProcessor.cpp
    namespace ssd {
        const int OFFSET = 1600;
        const int RANGE = 100;
    }
 
    namespace raw {
        /** Notice offset 3000, so allg
         * the ids are effectively +3000 in the his file*/
        const int OFFSET = 3000;
        const int RANGE = 1600;

        /** Notice that there is a space for 300 channels,
         * one 13 modules crate has 208 channels */
        const int D_RAW_ENERGY = 0;
        const int D_FILTER_ENERGY = 300;
        const int D_SCALAR = 600;
        const int D_TIME = 900;
        const int D_CAL_ENERGY = 1200;

        const int D_HIT_SPECTRUM = 1500;
        const int D_SUBEVENT_GAP = 1501;
        const int D_EVENT_LENGTH = 1502;
        const int D_EVENT_GAP = 1503;
        const int D_EVENT_MULTIPLICITY = 1504;
        const int D_BUFFER_END_TIME = 1505;
        const int DD_RUNTIME_SEC = 1506;
        const int DD_DEAD_TIME_CUMUL = 1507;
        const int DD_BUFFER_START_TIME = 1508;
        const int DD_RUNTIME_MSEC = 1509;
        const int D_NUMBER_OF_EVENTS = 1510;
        const int D_HAS_TRACE = 1511;
//	const int DD_RAW_V_CAL = 1512;
    

    }

    // in ImplantSsdProcessor.cpp
    namespace implantSsd {
        const int OFFSET = 6000;
        const int RANGE = 300;
    }

    // in MtcProcessor.cpp
    namespace mtc {
        const int OFFSET = 6500;
        const int RANGE = 100;
    } 

    //in IonChamberProcessor.cpp 
    namespace ionChamber{ 
        const int OFFSET = 6700;
        const int RANGE = 100;
    } 

    namespace position {
        const int OFFSET = 5000;
        const int RANGE = 600;
    } 

    namespace trace {
        const int OFFSET = 7500;
        const int RANGE = 150;
        const int traceBins = SA;

        namespace tracefilterer {
            const int DD_TRACE = 0;
            const int DD_FILTER1 = 1;
            const int DD_FILTER2 = 2;
            const int DD_FILTER3 = 3;
            const int DD_AVERAGE_TRACE = 4;
            const int DD_REJECTED_TRACE = 5;
            const int DD_ENERGY__BOARD_FILTER = 6;
            const int D_RATIO_BOARD_FILTER = 7;
            const int D_ENERGY1 = 8;

	  /*
	  // ---------- by Yongchi Xiao ----------
	  const int DD_MYTEST1 = 9;
	  const int DD_MYTEST2 =10;
	  //--------------------------------------*/
        }

        namespace doubletraceanalyzer {
	  const int D_ENERGY2 = 16;	
	  const int DD_DOUBLE_TRACE = 20;
	  const int DD_ENERGY2__TDIFF = 21;
	  const int DD_ENERGY2__ENERGY1 = 22;
	  const int DD_TRIPLE_TRACE = 30;
	  const int DD_TRIPLE_TRACE_FILTER1 = 31;
	  const int DD_TRIPLE_TRACE_FILTER2 = 32;
	  const int DD_TRIPLE_TRACE_FILTER3 = 33;
	  // added by SG
	  const int DD_DOUBLE_TRACE_FRONT=35;
	  const int DD_DOUBLE_TRACE_BACK=36;
	}

        namespace waveformanalyzer {
            const int DD_TRACES     = 40;
            const int D_CHISQPERDOF = 41;
            const int D_PHASE       = 42;
            const int DD_AMP        = 43;
        }
        
        // 1D-traces from the extracter
        namespace extracter {
            const int maxSingleTraces = 99;
            const int D_TRACE = 50;
        }
    } 
}


#endif // __DAMM_PLOTIDS_HPP_
