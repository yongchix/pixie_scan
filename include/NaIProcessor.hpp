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
using std::endl;

class NaIProcessor : public EventProcessor
{
public:
	NaIProcessor();
	virtual void DeclarePlots(void);
	virtual bool PreProcess(RawEvent &event); // by Yongchi Xiao; 10/09/2015, for method3 only
	virtual bool Process(RawEvent &rEvent);
private:
	struct NaIData {
		void Clear(void);
	} data;
};

#endif // __NAIPROCESSOR_HPP_
