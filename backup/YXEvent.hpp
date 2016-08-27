#include <iostream>
#include <vector>
#include <string>

using namespace std;

namespace YXEvent
{
  struct SimpleEv // usable for nai/pin
  {
    double time, energy;
    int channel;
    string type;
    SimpleEv(){time = -1; energy = -1; channel = -1; type = "null";}
    void Clear(){time = -1; energy = -1; channel = -1; type = "null";}
    void AssignValue(double t, double e, double c, string ty)
    {
      time = t;
      energy = e;
      channel = c;
      type = ty;
    }
    ~SimpleEv(){}
  };
    
  struct CorrEv
  {
    double time1, time2;
    double energy1, energy2;
    int channel1, channel2;
    string corrType;
    CorrEv(){time1 = -1; time2 = -1; energy1 = -1; energy2 = -1; channel1 = -1; channel2 = -1; corrType = "null";}
    void AssignValue(SimpleEv s1, SimpleEv s2)
    {
      time1 = s1.time;
      time2 = s2.time;
      energy1 = s1.energy;
      energy2 = s2.energy;
      channel1 = s1.channel;
      channel2 = s2.channel;
      corrType = s1.type + s2.type;
    }
    void Clear(){time1 = -1; time2 = -1; energy1 = -1; energy2 = -1; channel1 = -1; channel2 = -1; corrType = "null";}
    ~CorrEv(){}
  };
}
