/** \file NaIProcessor.hpp
 * 
 * 
 */

#ifndef __NAIPROCESSOR_HPP_
#define __NAIPROCESSOR_HPP_

#include "EventProcessor.hpp"
#include <iomanip>
#include <vector>
#include <utility>
#include <cmath>

using std::cout;
using std::setprecision;
using std::endl;


class NaIProcessor : public EventProcessor
{
 public:
  NaIProcessor();
  virtual void DeclarePlots(void);
  virtual bool PreProcess(RawEvent &event); // by Yongchi Xiao; 10/09/2015
  virtual bool Process(RawEvent &rEvent);
private:
  struct NaIData {
    void Clear(void);
  } data;
};

//----------- by Yongchi Xiao; 06/17/2015; 06/23/2015 -----------
/*
struct NaiEvent
{
  double energy; 
  int channel;
  double time;

  NaiEvent()
  {
    energy = -1;
    channel = -1;
    time = -1;
  }


  void Output()
  {
    std::cout << std::setprecision(6) << energy*2.266 << "keV;   Ch: " << channel << "; time: "
	      << std::setprecision(12) << time << std::endl; // by Yongchi Xiao; 06/21/2015

  }
  void Clear()
  {
    energy = -1;
    channel = -1;
    time = -1;
  }
  double NaiGetTime()
  {return time;}
  double NaiGetCalEnergy()
  {return energy;}
  double NaiGetPos()
  {return channel;}
};
*/
//--------------------------------------------------

#endif // __NAIPROCESSOR_HPP_
