/** \file CfdAnalysis.hpp
 * 
 * Class to use Cfd on traces
 */
#ifndef __CFDANALYZER_HPP_
#define __CFDANALYZER_HPP_

#include "TimingInformation.hpp"
#include "Trace.hpp"
#include "TraceAnalyzer.hpp"

class CfdAnalyzer : public TraceAnalyzer
{
 public:
    CfdAnalyzer();
    virtual void DeclarePlots(void) const;
    virtual void Analyze(Trace &, const std::string &, 
			 const std::string &);
    virtual ~CfdAnalyzer() {};
};

#endif // __CFDANALYZER_HPP_
