/*! \file SheCorrelator.cpp
 *
 * The SheCorrelator is designed to recontruct dssd events in a SHE experiment
 * and to correlate chains of alphas in dssd pixels
 */

#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>

#include "DammPlotIds.hpp"
#include "Dssd4JAEAProcessor.hpp"
#include "JAEACorrelator.hpp"
#include "DetectorDriver.hpp"
#include "Exceptions.hpp"
#include "Notebook.hpp"


using namespace std;

JAEAEvent::JAEAEvent() {
    energy_ = -1.0;
    energy2_= -1.0;
    time_= -1.0;
    mwpc_= -1;
    mwpcTime_= -1;
    has_beam_= false;
    has_veto_= false;
    has_escape_ = false;
    type_= unknown_jaea;
}

JAEAEvent::JAEAEvent(double energy, double energy2,double time, int mwpc, double mwpcTime, 
                   bool has_beam, bool has_veto, bool has_escape, JAEAEventType type /* = unknown*/)
{
    set_energy(energy);
    set_energy2(energy2);
    set_time(time);
    set_mwpc(mwpc);
    set_mwpcTime(mwpcTime);
    set_beam(has_beam);
    set_veto(has_veto);
    set_escape(has_escape);
    set_type(type);
}

JAEACorrelator::JAEACorrelator(int size_x, int size_y) 
{
    size_x_ = size_x + 1;
    size_y_ = size_y + 1;
    pixels_ = new deque<JAEAEvent>*[size_x_];
    for(int i = 0; i < size_x_; ++i)
        pixels_[i] = new deque<JAEAEvent>[size_y_];

}


JAEACorrelator::~JAEACorrelator() {
    for(int i = 0; i < size_y_; ++i) {
        delete[] pixels_[i];
    }
    delete[] pixels_;
}


bool JAEACorrelator::add_event(JAEAEvent& event, int x, int y, Plots& histo,bool hasVETO,bool hasNaI,bool hasPin) {

  //  cout << size_x_ << " " << size_y_ << endl;
  
  if (x < 0 || x >= size_x_) {
    stringstream ss;
    ss << "Requested event at non-existing X strip " << x << endl;
    throw GeneralWarning(ss.str());
  }
  if (y < 0 || y >= size_y_) {
    stringstream ss;
    ss << "Requested event at non-existing Y strip " << y << endl;
    throw GeneralWarning(ss.str());
  }
  
  if (event.get_type() == heavyIon_jaea){
    flush_chain(x, y, histo, hasVETO, hasNaI, hasPin);
    pixels_[x][y].push_back(event);
    //      cout << "heavyion" << endl;
    }
  if (event.get_type() == alpha_jaea){
    //    flush_chain(x, y, histo, hasVETO, hasNaI, hasPin);
    pixels_[x][y].push_back(event);
    //cout << "alpha" << endl;
  }


  return true;
}


bool JAEACorrelator::flush_chain(int x, int y, Plots& histo,bool hasVETO,bool hasNaI,bool hasPin){


    unsigned chain_size = pixels_[x][y].size();
    
    /** If chain too short just clear it */
    if (chain_size < 1) {
        pixels_[x][y].clear();
        return false;
    }
    
    JAEAEvent first = pixels_[x][y].front();

    /** Conditions for interesing chain:
     *      * starts with heavy ion implantation
     *      * has ion + fission_jaea
     *      * or includes at least two alphas
     */

    /** If it doesn't start with hevayIon, clear and exit**/
    if (first.get_type() != heavyIon_jaea) {
        pixels_[x][y].clear();
        return false;
    } 

    /* If it is greater than 1 element long, check if the last is fission,
     *  if not - clear and exit**/
    //    if (chain_size <= 2 && pixels_[x][y].back().get_type() != fission_jaea) {
    //    pixels_[x][y].clear();
    //    return false;
    // }

    stringstream ss;

    time_t wallTime = DetectorDriver::get()->GetWallTime(first.get_time());
    string humanTime = ctime(&wallTime);
    humanTime.erase(humanTime.find('\n', 0), 1);
    ss << humanTime << "\t X = " << x <<  " Y = " << y << endl;

    int alphas = 0;
    double alphaE[6]={0};
    double alphaE2[6]={0};
    double alphaTime[6]={0};
    double dt[6],dt2[6],dt12[6];
    double BeamE=0, BeamTime=0;
    double mwpcTime = 0;
    static int ctr=0, ctra=0;
    int ittr=0,i=0;
    for (deque<JAEAEvent>::iterator it = pixels_[x][y].begin();
         it != pixels_[x][y].end();
         ++it)
      { 	
	int type = (*it).get_type();
	double energy = (*it).get_energy();
	double energy2 = (*it).get_energy2();
	double time = (*it).get_time()*(Globals::get()->clockInSeconds()); 
	mwpcTime = (*it).get_mwpcTime()*(Globals::get()->clockInSeconds());
	
	if( type == heavyIon_jaea) { // implantation
	  BeamE=energy;
	  BeamTime=time;
	  if(x==24 && y==12){
	    //cout << "====heavyIon " <<BeamE << " " << BeamTime << " " << x << " " <<  y << endl;
	  }
	}    	
	if (type == alpha_jaea) {    // decay
            alphas += 1;
	    if(alphas <= 2) {
		alphaE[alphas]=energy;
		alphaE2[alphas]=energy2;

		alphaTime[alphas] = time;
		//cout << "=====decays "  << alphas << " " << alphaE[alphas] << " " << alphaE2[alphas] << endl;
		if(alphas==2){
		  //		if(x==24 && y==12){
		  //cout << "decays" <<alphas << " " <<  alphaE[alphas] << " " << alphaTime[alphas]<< " " << x << " " << y  << endl;
		}
		
	    }
        }
	
	
	if(BeamTime>0 && alphaTime[1]>0){
	  //cout << "Timing "<< BeamTime << " " << alphaTime[1] << " " << dt << endl;
	  const unsigned int NumGranularities = 8;
	  for (unsigned int i = 0; i < NumGranularities; i++) {
	    const double timeResolution[NumGranularities] = 
	      {10e-9, 100e-9, 400e-9, 1e-6, 100e-6, 1e-3, 10e-3, 100e-3};
	    
	    
	    //100000000 is temporary value
	    dt[i]=100000000*(alphaTime[1]-BeamTime)*(Globals::get()->clockInSeconds()/timeResolution[i]);
	    histo.Plot(dammIds::dssd4jaea::DD_ENERGY_DECAY_TIME_GRANX+i,2*alphaE[1],int(dt[i])); // 750-759
	    
	    if(i == 4
	       && int(dt[i]) >= 10 && int(dt[i]) <=250
	       && alphaE[1] > 100
	       && alphaE2[1] > 100
	       )
	      {
		histo.Plot(8, alphaE[1], x); // 708, front
		histo.Plot(9, alphaE2[1], y); // 709, back
	      } // by Yongchi Xiao; 10/22/2015
	    
	    

	    if(hasVETO){
              histo.Plot(dammIds::dssd4jaea::DD_ENERGY_DECAY_TIME_GRANX_VETO+i,2*alphaE[1],int(dt[i]));
            }
            if(!hasVETO){
              histo.Plot(dammIds::dssd4jaea::DD_ENERGY_DECAY_TIME_GRANX_NOVETO+i,2*alphaE[1],int(dt[i]));
            }
            if(hasNaI){
              histo.Plot(dammIds::dssd4jaea::DD_ENERGY_DECAY_TIME_GRANX_NAI+i,2*alphaE[1],int(dt[i]));
            }
            if(!hasNaI){
              histo.Plot(dammIds::dssd4jaea::DD_ENERGY_DECAY_TIME_GRANX_NONAI+i,2*alphaE[1],int(dt[i]));
            }
            if(hasPin){
              histo.Plot(dammIds::dssd4jaea::DD_ENERGY_DECAY_TIME_GRANX_PIN+i,2*alphaE[1],int(dt[i]));
            }
            if(!hasPin){
              histo.Plot(dammIds::dssd4jaea::DD_ENERGY_DECAY_TIME_GRANX_NOPIN+i,2*alphaE[1],int(dt[i]));
            }

	    




	    if(i==3){
	      // Cs condition
	      if(int(dt[i])<50){
		//cout << "Cs: " << alphaE[1] << " " << y << endl; 
		histo.Plot(dammIds::dssd4jaea::DD_IMPLANT_CSGATE,alphaE[1],y);
	      }
	      // I condition
	      if(int(dt[i])>50 && int(dt[i])<200){
		//cout << "I: "<<alphaE[1] << " " << y << endl; 
		histo.Plot(dammIds::dssd4jaea::DD_IMPLANT_IGATE,alphaE[1],y);
	      }

	      
	    }
	    
	    //cout << "check " << BeamE << " " << int(dt[i]) << endl;
	    if(alphaTime[2]>0){
	      dt2[i]=100000000*(alphaTime[2]-BeamTime)*(Globals::get()->clockInSeconds()/timeResolution[i]);
	      dt12[i]=100000000*(alphaTime[2]-alphaTime[1])*(Globals::get()->clockInSeconds()/timeResolution[i]);
	      //	      histo.Plot(dammIds::dssd4jaea::DD_ENERGY_DECAY_TIME_GRANX+i,2*alphaE[1],int(dt[i]));
	      //cout << "Decay E1: "<< alphaE[1] << " Decay E2: " << alphaE[2]  << " Time " << int(dt12[i]) << endl;
	      
	      histo.Plot(dammIds::dssd4jaea::DD_ENERGY_DECAY_TIME_GRANX_SECOND+i,2*alphaE[1],int(dt12[i]));	     
	      

	      

	      


	      
	    }

	    

	  }
	  
	  // front 788 1st decay and 2nd decay
	  if(alphaE[1]>0 && alphaE[2]>0){
	    histo.Plot(dammIds::dssd4jaea::DD_ENERGY_DECAY12,alphaE[1],alphaE[2]);
	  }
	  // back 789 decay 1st decay and 2nd decay
	  if(alphaE2[1]>0 && alphaE2[2]>0){
	    histo.Plot(dammIds::dssd4jaea::DD_ENERGY_DECAY12_2,alphaE2[1],alphaE2[2]);
	  }

	  
	}
	
	human_event_info((*it), ss, first.get_time());
        ss << endl;

    }

    pixels_[x][y].clear();
    
    return true;
}

// Save event to file
void JAEACorrelator::human_event_info(JAEAEvent& event, stringstream& ss,
                                     double clockStart) {
    string humanType;
    switch (event.get_type()) {
        case alpha_jaea:   humanType = "A";
                        break;
        case heavyIon_jaea: humanType = "I";
                        break;
        case fission_jaea: humanType = "F";
                        break;
        case lightIon_jaea: humanType = "L";
                        break;
        case unknown_jaea: humanType = "U";
                        break;
    }

    ss << fixed 
       << humanType 
       << " " 
       << setprecision(0) << setw(12) << event.get_energy()
       << " " 
       << setprecision(3) << setw(12) 
       << (event.get_time() - clockStart) / 1.0e7  
       << " M" << event.get_mwpc() 
       << "B" << event.get_beam() 
       << "V" << event.get_veto()
       << "E" << event.get_escape();

}
