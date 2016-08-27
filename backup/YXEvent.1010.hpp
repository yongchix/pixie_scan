#include <iostream>
#include <vector>
#include <string>

using namespace std;


class SimpleEv // usable for nai/pin
{
private: 
  double time, energy;
  int channel;
  string type;
  
public: 
  SimpleEv(){time = -1; energy = -1; channel = -1; type = "null";}
  ~SimpleEv(){}
  void Clear(){time = -1; energy = -1; channel = -1; type = "null";}
  void AssignValue(double t, double e, double c, string ty)
  {
    time = t;
    energy = e;
    channel = c;
    type = ty;
  }
  string GetType(){return type;}
  double GetTime(){return time;}
  int GetChannel(){return channel;}
  double GetEnergy(){return energy;}
};

class CorrEv
{
private: 
  vector<SimpleEv> corrEvents;
  string corrType;
  bool hasPin;
  
public: 
  CorrEv(){corrType = "null"; hasPin = false; }
  ~CorrEv(){}
  void AddEvent(SimpleEv s)
  {
    if(s.GetType().compare("nai") == 0)
      {
	if(corrType.compare("null") == 0)
	  corrType = s.GetType();
	else
	  corrType += s.GetType();
	corrEvents.push_back(s);
      }
  }
  void TriggerPin(){hasPin = true;}
  void Clear()
  {
    corrEvents.clear();
    corrType = "null";
    hasPin = false;
  }
  string GetCorrType(){ return corrType;}
  bool GetPinInfo(){ return hasPin;}
  double GetTime(){ return corrEvents.at(0).GetTime();}
}; 

//CorrEv naiPair; 
