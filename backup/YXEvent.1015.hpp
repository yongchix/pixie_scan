#include <iostream>
#include <vector>
#include <string>

using namespace std;


class SimpleEv // usable for nai/pin, one-sided DSSD signal
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
  int auxNai;
    
public: 
  CorrEv(){corrType = "null"; auxNai = 0;}
  ~CorrEv(){}
  void AddEvent(SimpleEv s)
  {
    if(corrType.compare("null") == 0)
      corrType = s.GetType();
    else
      corrType += s.GetType();
    corrEvents.push_back(s);
  }
  void AddAuxNai(int n){auxNai += n;}
  void Clear()
  {
    corrEvents.clear();
    corrType = "null";
    auxNai = 0;
  }
  string GetCorrType(){ return corrType;}
  string GetIndivType(int i){ return corrEvents.at(i).GetType();}
  double GetIndivTime(int i){ return corrEvents.at(i).GetTime();}
  double GetIndivEnergy(int i){return corrEvents.at(i).GetEnergy();}
  int GetAuxNai(){ return auxNai;}
}; 

//CorrEv naiPair; 
