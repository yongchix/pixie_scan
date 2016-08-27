/** \file McpProcessor.cpp
 * Mcp class takes the signals from the mcp detector
 *   and calculates a 2D position based on the readout
 * SNL 2-2-08, Modified DTM 9-09
 */
#ifdef useroot
#include <TTree.h>
#endif

#include "DammPlotIds.hpp"
#include "McpProcessor.hpp"
#include "RawEvent.hpp"

using std::string;
using std::vector;
using namespace dammIds::mcp;
using namespace std;

namespace dammIds {
    namespace mcp {	
      const int D_POSX   = 1;
      const int D_POSY   = 2;
      const int DD_POSXY = 3;
      const int D_AMP1 = 4;
      const int D_AMP2 = 5;
      const int DD_XE_YE = 6;
      const int D_BEAM=9;
    }
}


void McpProcessor::McpData::Clear(void)
{
  for (size_t i=0; i < nPos; i++)
    raw[i] = 0;
  xpos = ypos = 0.;

  mult = 0;
}

McpProcessor::McpProcessor(void) : EventProcessor(OFFSET, RANGE, "mcp")
{
  associatedTypes.insert("mcp");
}

void McpProcessor::DeclarePlots(void)
{
  
  const int posBins   = S9; //
  const int posBins2  = SE; //1024  
  const int posBins2D = S9; // 512
  const int energyBins = SE;
  const int posBins2D2 = SA;
  // 921
  DeclareHistogram1D(D_POSX, posBins2, "Horiz. - MCP Left");
  // 922
  DeclareHistogram1D(D_POSY, posBins2, "Vert. - MCP Left");
  // 923
  DeclareHistogram2D(DD_POSXY, posBins2D2, posBins2D2*2, "MCP Left 2D");

  // 929
  DeclareHistogram1D(D_BEAM,energyBins,"Beam current");

}


bool McpProcessor::Process(RawEvent &event)
{
  
  if (!EventProcessor::Process(event))
    return false;

  // temperarily commented out by Yongchi Xiao; 05/18/2016
  /*
  double qTotal, qRight, qTop, Cathode;
  
  static const vector<ChanEvent*> &mcpEvents = sumMap["mcp"]->GetList();

  data.Clear();
  
  for (vector<ChanEvent*>::const_iterator it = mcpEvents.begin();
       it != mcpEvents.end(); it++) {
      ChanEvent *chan = *it;

      string subtype   = chan->GetChanID().GetSubtype();
      double calEnergy = chan->GetCalEnergy();
      double mcpTime   = chan->GetTime();
      
      double current=0;
      
      if(subtype == "beam"){
	current = calEnergy;
	//plot(D_BEAM,current);
	//cout << subtype << " " << current << endl;
      }

      if(subtype == "1time") 
	{
	  // do nothing
	  Cathode=calEnergy;
	  //cout << "1time :: " << data.mult << endl; // by Yongchi Xiao;
	} else if (subtype == "1position1")
	{
	  data.raw[0] = calEnergy;
	  data.time[0] = mcpTime;
	  data.mult++;
	  //cout << "1position1 :: " << data.mult << endl;// by Yongchi Xiao
	} else if (subtype == "1position2")
	{
	  data.raw[1] = calEnergy;
	  data.mult++;
	  data.time[1] = mcpTime;
	  //cout << "1position2 :: " << data.mult << endl; // by Yongchi Xiao
	} else if (subtype == "1position3")
	{
	  data.raw[2] = calEnergy;
	  data.mult++;
	  data.time[2] = mcpTime;
	  //cout << "1position3 :: " << data.mult << endl;// by Yongchi Xiao
	} else if (subtype == "1position4")
	{
	  data.raw[3] = calEnergy;
	  data.time[3] = mcpTime;
	  data.mult++;
	  //cout << "1position3 :: " << data.mult << endl;// by Yongchi Xiao
	}
  }
    
  // calculation of position from charge collected on the four corners
  // magic numbers here
  data.raw[0] *= 1.3;
    
  qTotal = data.raw[0] + data.raw[1] + data.raw[2] + data.raw[3];
  qRight = data.raw[3] + data.raw[0];
  // qLeft   = data.raw[2] + data.raw[1];  
  qTop   = data.raw[0] + data.raw[1];
  // qBottom = data.raw[2] + data.raw[3];
    
  data.xpos = (qRight / qTotal) * 512. - 75; //horizontal MCP pos
  data.ypos = (qTop   / qTotal) * 512. - 75; //vertical MCP pos
  // qLeft, qBottom not used
  
  if (data.mult >= 2) 
    {
      using namespace dammIds::mcp;
      
      double tDiff=  data.time[1] - data.time[0]+20;
      
      plot(D_POSX, tDiff); //951   
      plot(D_POSY, Cathode);// 952, Cathode = calEnergy for subtyoe=="ltime"; by YX
      
      // old divided value
      //data.raw[2]=data.raw[2]/5.;
      plot(DD_POSXY,data.raw[2]/5,Cathode);
      
      //plot(DD_POSXY,data.raw[2]/4,Cathode);
      
      // for current monitor
      
      
    }
  */


  
  
  EndProcess();
  //  return (data.mult == 3); 
  return true;


  //cout << "mult: " << data.mult << endl;// by Yongchi Xiao
}

#ifdef useroot
bool McpProcessor::AddBranch(TTree *tree)
{
  if (tree) {
    TBranch *mcpBranch = 
      tree->Branch(name.c_str(), &data, "raw[4]/D:xpos:ypos:mult/I");

    return (mcpBranch != NULL);
  } 

  return false;
}

void McpProcessor::FillBranch(void)
{
  if (!HasEvent())
    data.Clear();  
}
#endif //useroot
