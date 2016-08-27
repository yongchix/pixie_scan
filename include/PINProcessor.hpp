/** \file NaIProcessor.hpp
 * 
 * 
 */

#ifndef __PINPROCESSOR_HPP_
#define __PINPROCESSOR_HPP_

#include <iomanip>

#include "EventProcessor.hpp"
#include "ChanEvent.hpp"

using namespace std;

class PINProcessor : public EventProcessor
{
public:
	PINProcessor();
	virtual void DeclarePlots(void);
	virtual bool Process(RawEvent &rEvent);
	std::vector<ChanEvent*> GetFrontPinEvents() {
		return pinFrontEvents;
	}
	std::vector<ChanEvent*> GetBackPinEvents() {
		return pinBackEvents;
	}
	std::vector<ChanEvent*> GetPinEvents() {
		return pinEvents;
	}
private:
	struct PINData {
		void Clear(void);
	} data;
	std::vector<ChanEvent*> pinBackEvents, pinFrontEvents, pinEvents;
};

#endif // __PINPROCESSOR_HPP_
