#include <iostream>
#include <vector>
#include <string>

using namespace std;


struct SimpleEvent // usable for nai/pin, one-sided DSSD signal
{
  double time, energy;
  int channel;
  string type;
  
  SimpleEvent(){time = -1; energy = -1; channel = -1; type = "null";}
  ~SimpleEvent(){}
  void Clear(){time = -1; energy = -1; channel = -1; type = "null";}
  void AssignValue(double t, double e, double c, string ty)
  {
    time = t;
    energy = e;
    channel = c;
    type = ty;
  }
  // 03/31/2016
  bool operator<(const SimpleEvent &se) const{
    return time < se.time;
  }
};

class CorrFlag
{
private: 
  double time;
  bool correlated;
public:
  CorrFlag(){
    time = -1;
    correlated = false;
  }
  ~CorrFlag(){}
  void Clear(){
    time = -1;
    correlated = false;
  }
  void Mark(double time_, bool corr_){
    time = time_;
    correlated = corr_;
  }
  bool CheckCorr(){
    return correlated;
  }
  double GetTime(){
    return time;
  }

}; 


