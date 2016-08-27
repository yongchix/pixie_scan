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
#include "Dssd4SHEProcessor.hpp"
#include "SheCorrelator.hpp"
#include "DetectorDriver.hpp"
#include "Exceptions.hpp"
#include "Notebook.hpp"


using namespace std;

SheEvent::SheEvent() {
    energy_ = -1.0;
    time_= -1.0;
    mwpc_= -1;
    mwpcTime_= -1;
    has_beam_= false;
    has_veto_= false;
    has_escape_ = false;
    type_= unknown;
}

SheEvent::SheEvent(double energy, double time, int mwpc, double mwpcTime, 
                   bool has_beam, bool has_veto, bool has_escape, SheEventType type /* = unknown*/)
{
    set_energy(energy);
    set_time(time);
    set_mwpc(mwpc);
    set_mwpcTime(mwpcTime);
    set_beam(has_beam);
    set_veto(has_veto);
    set_escape(has_escape);
    set_type(type);
}


SheCorrelator::SheCorrelator(int size_x, int size_y) 
{
    size_x_ = size_x + 1;
    size_y_ = size_y + 1;
    pixels_ = new deque<SheEvent>*[size_x_];
    for(int i = 0; i < size_x_; ++i)
        pixels_[i] = new deque<SheEvent>[size_y_];

}


SheCorrelator::~SheCorrelator() {
    for(int i = 0; i < size_y_; ++i) {
        delete[] pixels_[i];
    }
    delete[] pixels_;
}


bool SheCorrelator::add_event(SheEvent& event, int x, int y, Plots& histo) {

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

    if (event.get_type() == heavyIon)
        flush_chain(x, y, histo);

    pixels_[x][y].push_back(event);

    if (event.get_type() == fission)
        flush_chain(x, y, histo);
    
    return true;
}


bool SheCorrelator::flush_chain(int x, int y, Plots& histo){

    unsigned chain_size = pixels_[x][y].size();

    /** If chain too short just clear it */
    if (chain_size < 2) {
        pixels_[x][y].clear();
        return false;
    }

    SheEvent first = pixels_[x][y].front();

    /** Conditions for interesing chain:
     *      * starts with heavy ion implantation
     *      * has ion + fission
     *      * or includes at least two alphas
     */

    /** If it doesn't start with hevayIon, clear and exit**/
    if (first.get_type() != heavyIon) {
        pixels_[x][y].clear();
        return false;
    } /*else {
	VRecoilE=first.get_energy();
	VRecoilTime=first.get_time()*2e-7;
    }*/

    /* If it is greater than 1 element long, check if the last is fission,
     *  if not - clear and exit**/
    if (chain_size <= 2 && pixels_[x][y].back().get_type() != fission) {
        pixels_[x][y].clear();
        return false;
    }

    stringstream ss;

    time_t wallTime = DetectorDriver::get()->GetWallTime(first.get_time());
    string humanTime = ctime(&wallTime);
    humanTime.erase(humanTime.find('\n', 0), 1);
    ss << humanTime << "\t X = " << x <<  " Y = " << y << endl;

    int alphas = 0;
    double alphaE[6]={0};
    double alphaTime[6]={0};
    double fissionE =0, fissionTime =0;
    double VRecoilE=0, VRecoilTime=0;
    double mwpcTime = 0;
    static int ctr=0, ctra=0;
    int ittr=0,i=0;
    for (deque<SheEvent>::iterator it = pixels_[x][y].begin();
         it != pixels_[x][y].end();
         ++it)
    { 	// Process Correlations in the SHE Event. Make a logic for a good event and pass it to the plot routines. 
	int type = (*it).get_type();
	double energy = (*it).get_energy();
	double time = (*it).get_time()*(Globals::get()->clockInSeconds());//units in s *1e-3 gives units in us
	mwpcTime = (*it).get_mwpcTime()*(Globals::get()->clockInSeconds());
	 


	if( type == heavyIon) {
	    VRecoilE=energy;
	    VRecoilTime=time;
	}    	
	if(type == fission) {
	    fissionE = energy;
	    fissionTime = time;
	   // cout << time - mwpcTime << endl;
	}
	if (type == alpha) {
            alphas += 1;
	    if(alphas <= 6) {
		alphaE[alphas]=energy;
		alphaTime[alphas] = time; 
	    }
        }
	//tbd process the unknown events in correlation with the good chains.
	human_event_info((*it), ss, first.get_time());
        ss << endl;

    }

    pixels_[x][y].clear();


    if (alphas >= 2 && (alphaTime[1]-mwpcTime) < 16){
	//cout << "INTERESTING alpha!!" << endl;
//	Notebook::get()->report(ss.str()); 
	ittr =0;
        while(ittr < 2 ) { 
	    histo.Plot(dammIds::dssd4she::DD_CHAIN_NUM_ALPHA,VRecoilE/10,ctra);
	    ittr++; 
	}
	/*i=1;
	while (i<6) {
		cout << alphaE[i] << " " << alphaTime[i] << " " <<VRecoilTime << " " << ctra << endl;
		i++;
	
	}*/

	i=1;
	while (alphaE[i]>0){
	    ittr =0;
	    //cout <<alphaTime[i] << " " << VRecoilTime << endl;
     	    histo.Plot(dammIds::dssd4she::DD_CHAIN_ALPHA_V_ALPHA,alphaE[1]/10,alphaE[i+1]/10);
	    while (ittr <= (alphaTime[i]-mwpcTime)*1e6 && ittr <16000) {
		histo.Plot(dammIds::dssd4she::DD_CHAIN_NUM_ALPHA,alphaE[i]/10,ctra);
		ittr++;
	    }
	    i++;
	}

	ctra++;      
    }
    if (fissionE > 0 && fissionTime-mwpcTime <=16 ) {
	cout << "INTERESTING fission!!"  << endl;
	// Messenger m;
    //    m.run_message(ss.str());	  
	
	Notebook::get()->report(ss.str());
	ittr =0;
        while(ittr < 2 ) { 
	    histo.Plot(dammIds::dssd4she::DD_CHAIN_NUM_FISSION,VRecoilE/10,ctr);
	    ittr++; 
	}
	ittr =0;
	//cout << round(fissionTime -VRecoilTime) << endl;
	while (ittr < (fissionTime - mwpcTime)*1e6 && ittr <16000 ) {
	    histo.Plot(dammIds::dssd4she::DD_CHAIN_NUM_FISSION,fissionE/100,ctr);
	    ittr++;
	}
	ctr++;
	       
	//cout << x << " " << y << " " << VRecoilE << " " << fissionE <<  " " << fissionTime - VRecoilTime << endl;
    }
    
      
	
    /*    ss << ' ' << ctr << ' ' << Energy1;
	//ss <<' ' << Time1 << ' ' << Energy2; 
	//ss << ' ' << Time2;
        Notebook::get()->report(ss.str());

        //m.run_message(oss.str());
    */

    return true;
}

// Save event to file
void SheCorrelator::human_event_info(SheEvent& event, stringstream& ss,
                                     double clockStart) {
    string humanType;
    switch (event.get_type()) {
        case alpha:   humanType = "A";
                        break;
        case heavyIon: humanType = "I";
                        break;
        case fission: humanType = "F";
                        break;
        case lightIon: humanType = "L";
                        break;
        case unknown: humanType = "U";
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

