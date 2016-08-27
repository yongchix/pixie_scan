/** \file doc_for_doxygen.h
  \brief Dummy file for Doxygen structure
  
  This file contains the structure of the doxygen generated
  main page and subpages.  The documentation is lengthy and
  I wanted to separate it from a specific source file for
  easier readability
*

/** \mainpage Pixie16 Analysis 

  \author Sean Liddick
  \author David Miller
  \version code - 2.0, manual - beta v2

  This manual documents the workings of the Pixie16 analysis code.
  The manual is divided into two parts: a Users Manual and a Reference
  Manual.
  
  \section useman User Manual

  - An \subpage philosophy of the pixie16 code
  - A \subpage primer for those who have never seen C++
  - A \subpage directory description
  - A \subpage quickstart guide
  - An \subpage faq
  - Some \subpage codewarn about the code
  
  \section refman Reference Manual

  - An \subpage introduction to the pixie16 analysis
  - A list of \subpage dettype used in the analysis
  - A list of \subpage objects (both detector and others) used in the analysis
  - Experiment Independent Event Processing
       -# \subpage %PixieStd
       -# \subpage DetectorDriver
  - Experiment Dependent Event Processing
       -# \subpage MTC
       -# \subpage RMS
  - A list of \subpage changes in the most recent version

  \section devel Developer Manual (under construction)
  - A list of \subpage flags used in the Makefile
  - A method to \subpage add
  - Future \subpage improvement ideas
  - Some \subpage names for sanity
*/

/*********************************************************************/

/** \page introduction Introduction
  There are two main files which control the pixie16 analysis flow,
  PixieStd.cpp and DetectorDriver.cpp both of which are located
  in the src directory.  The functions in PixieStd.cpp reconstruct a
  complete data spill from the pixie16 modules and, for each spill, 
  create a list of channels that triggered.  The list is then sorted
  according to time and those channels which occur close to each other
  in time are grouped together into events.  The event is then passed
  into DetectorDriver.cpp for processing.

  DetectorDriver.cpp receives each event and calibrates the energy for
  each channel in the event.  The raw and calibrated energies are plotted
  if the appropriate damm spectra have been created in DeclareHistogram.cpp.
  Lastly, experiment specific analysis is performed. Experiment specific
  analysis for generic MTC and RMS experiments remains unimplemented in the
  current version.
    
*/

/*********************************************************************/

/** \page PixieStd Event Creation - PixieStd
  This page describes PixieStd.cpp used for event creation
  and should not have to be altered on a regular basis.  PixieStd 
  includes five important functions which control event creation:
  - hissub_() - The interface between SCAN and the C++ analysis. 
  Reassemble a complete pixie16 spill of data.
  - InitMap() - Analysis initialization
  - New map file info - An alternative way to define the analysis  
  - hissub_sec() - extract channel information from the raw pixie16
  data spill and sort the triggered channels based on time.
  - ScanList() - Group channels into events and send to detector_driver
  for processing
  - HistoStats() - plot various permanent spectra such as the run
  time and time between events.

  \section hissub
  
  PixieStd is the interface between the HRIBF scan program and the
  Pixie16 C++ analysis.  hissub_() is called from scan
  with two variables.  The first, ibuf, is a pointer to the array
  where the raw data is stored.  The second, nhw, describes the
  amount of raw data that is stored in the array.

  In a typical experiment with the new pixie16 readout (default),
  all Pixie16 modules are read out when one module has reached its
  maximum number of events or FIFO threshold.  This is programmed
  during experimental setup.  This spill of data is divided into
  smaller chunks of data for transmission across the network.  hissub_()
  collects the chunks of data to reassemble the complete spill.

  \anchor pixiestdinit
  \section Initialization
  
  Following the reconstruction of the first spill, and before
  continuing the analysis, the analysis is initialized with InitMap()

  \subsection InitMap
  When the first spill has been reconstructed and before analysis
  begins, InitMap() is called to initialize all routines and variables
  used in the scanning process.  First the DetectorLibrary
  class is queried to determine the valid detector types that can be used
  in the analysis.  This is accomplished with the following line of code
  in InitMap()
  \code
  // file pixie_std.cpp - InitMap()
  const set<string> &knownDets = driver.GetKnownDetectors()
  \endcode
  The detector types that are currently valid for analysis are
  dssd_front, dssd_back, ge, idssd_front, timeclass, si, position,
  scint, mcp, mtc, generic, ssd, vandle, pulser, logic, and ion_chamber.  
  The detector types position and idssd_front have to date been 
  exclusively used for fragmentation experiments at MSU (and have no
  implementation in the current version of code).

  Next the file map.txt is opened to determine a unique detector
  that is connected to each Pixie16 module.  Please see
  \subpage mapfile "map.txt" for more information.  For every
  line in the map.txt file the associated parameters are read into
  an identifier variable called id and added to the variable
  modChan (a vector containing a channel identifier for each
  line in the map.txt) which contains information about the 16 channels
  for each Pixie module.

  If the entry is to be ignored (by setting the detector type to
  "ignore"), a dummy entry is inserted into the modChan variable. Channels
  which do not appear in the map are also ignored based on the default
  detector type "" in the Identifier constructor. The position of each 
  channel identifier in modChan is retrieved dynamically based on a
  combination of the channel and module number, nominally
  \verbatim
  module number * 16 + channel number.
  \endverbatim
  Thus the identifier corresponding to module #2 and channel #10
  (id = 42) must be stored in modChan element 42 to allow the retrieval
  to work properly. The file map.txt no longer requires sequential ordering
  of channels nor the presence of all channels for a particular module.

  After reading in map.txt, the individual trace and event processors 
  contained in the DetectorDriver are initialized.
  
  For each Pixie16 channel that has not been set to "ignore," the
  detector type is checked against the known detectors. If the detector
  type is not listed among the known detectors, the program will 
  terminate to provide the user an opportunity to fix the problem.

  \section New map file info
  Alternatively, the detectors can be defined in a new map file format
  through the "map2.txt" file. Each line contains information in order about
  the detectors
  module number, channel number, type, subtype, location
  followed by arbitrary integer tags identified by a string
  e.g. traceanalysis=10 virtual start width=100 
  If no value is given, 1 is assumed.
  
  Module and channel fields can contain ranges (2-4) or asterisks to indicate
  all unfilled channels. ('e' and 'o' can also be provided to limit to even
  or odd values.) Single channel specifiers will be processed first.
  
  Any field after type can be left unspecified:
    If subtype is blank or "--", the subtype is assumed to be the same as type
    If location is blank, numbering will begin as low as available (0 start)
    Tags are often blank.
  The code is smart enough to recognize (based on whether the field is a
  number or a string) what you mean except in the case when location and
  subtype would be left blank. Here use the "--" shorthand for the subtype.
  Locations for a range are given sequentially looping over channels and then
  over modules.

  \section hissub_sec
  
  Once the spill is reassembled, hissub_sec() scans the
  reconstructed buffer to extract individual triggered channels using
  the ReadBuffData() function. The individual channels are inserted
  into the list of events and once the spill has been completely scanned 
  the list is sorted based on the event time of each individual channel.

  \section ScanList
  
  After the list of channels that triggered in the current spill has
  been sorted by time, ScanList() is called.  This function scans the
  list of channels and groups channels with similar times togther into
  events. Starting from the begining of the list and continuing to the end,
  an individual channel event time is compared with the previous channel
  event time to determine if they occur within the time period defined by
  the eventWidth variable (in units of 10 ns) and 
  - if yes - the two channels are considered to belong to the same
  event and the current channel is added to the list of channels
  in rawEvent.
  - if no - the previous rawEvent is considered complete and is
  sent to the DetectorDriver for processing.  Once processing
  is finished, the rawEvent variable is clared and the current
  channel is placed inside it to begin a new event.

*/

/*********************************************************************/

/** \page DetectorDriver Event Processing - DetectorDriver
  This page describes the aspects of DetectorDriver.cpp used for the RawEvent
  which remain consistent between experiments and includes five different
  main functions
  - Init() - initialize the detector driver
  - ProcessEvent() - processes a complete event passed from ScanList()
  - ThreshAndCal() - check threshold and calibrate the raw channel energies
  - PlotRaw() - plot raw experimental energies
  - PlotCal() - plot calibrated energies

  \section Initialization
  As mentioned in PixieStd.cpp \ref pixiestdinit "initialization"
  section, the DetectorDriver initialization routine
  (DetectorDriver::Init()) is called with the detector types that have
  been identified from the reading of the map.txt. This collection of 
  used types is also provided for the raw event initialization which
  creates a DetectorSummary for each detector type and places it into
  the map of detector summaries sumMap using the following code:
  \code
  //file RawEvent.cpp 
  bool DetectorDrive::Init(const set<string> &usedTypes,
                           const set<string> &usedSubtypes) {
  DetectorSummary ds;
  ds.Zero();
  
  for(set<string>::const_iterator it = usedTypes.begin();
     it != usedTypes.end(); i++) {
     ds.SetName(*it); 
     sumMap.insert(make_pair(*it,ds));
  }
  \endcode

  The detector and analysis specific objects are initialized through
  their respective routines.
  
  \subsection ReadCal
  The calibrations are read in during the detector driver
  initialization (DetectorDriver::ReadCal()) from the file cal.txt. 
  Please see \subpage calfile "cal.txt" for detailed information.

  For each line in cal.txt, the calibration is read into the calibration
  vector.  Each calibration lines is compared to match to a particular
  channel identifier as defined in the file map.txt. If a channel
  in map.txt does not exist in cal.txt, the program will
  issues an error and terminate. All channels not set to "ignore"
  must have a valid calibration.  The information from the map.txt
  and cal.txt are combined and output is printed to the screen
  for visual inspection.

  \section ProcessEvent
  ProcessEvent() accesses the global RawEvent and processes it.  
  The analysis loops over all channels contained within the RawEvent
  and checks the threshold, calibrates the energy, plots both the raw and
  calibrated energies (if the appropriate damm spectra have been
  created), and performs a trace analysis if a trace is present.

  \subsection ThreshAndCal
  ThreshAndCal() receives a pointer to the current channel event.
  The channel is calibrated using the element at "module number * 16 +
  channel number" in the calibration vector.  After calibration
  the DetectorSummary for this detector type is updated.
  
  \subsubsection Trace
  If a trace is included with the channel object then the trace
  analysis is invoked through the TraceAnalyzer
  
  \subsection PlotRaw
  The raw energy from the channel is plotted into the appropriate
  damm spectra as long as it was created in DeclareHistogram.cpp
  
  \subsection PlotCal
  The calibrated energy from the channel is plotted into the
  appropriate damm spectrum if it was created in DeclareHistogram.cpp
  
  \warning There is no default calibration
  \warning Provide a calibration line for every channel not
  set to "ignore" in the map.txt file.  Failure to do so will
  cause the program to terminate so the user can fix the problem.
*/

/*********************************************************************/

/** \page MTC MTC
  This page describes the experiment specific analysis that is
  conducted when running a beta decay experiment with the MTC

  \section detreq Detector Types
  For an MTC experiment the following detector types and sub types
  are used.
  - timeclass - MTC1 - this detector type corresponds to an MTC signal.
  It is in a group with other detectors which provide time information.
  There is no significance in the energy of the detected MTC signal
  only its presence.
  - timeclass - MTC2 - this detector type and sub type are optional
  but typically represent a delayed MTC signal so that at least one
  of the MTC channels should register the tape movement.
  - scint - beta1 and beta2 - beta detectors around the beam pipe
  for suppresion of background.
  - ge - clover - a set of clover detectors for gamma ray detection.
  - mcp - mcp - an optinal MCP detector use for timing decays with
  respect to an ion implantation.

  \section logproc Event Processing
  The MTC, MCP, and beta detectors are checked to see if any of them
  fired in the event.  If any of the three detectors fired the time
  of the respective channel is retrieved.

  \warning Currently no correlation is implemented for the MTC and
    beta scintillators, the former process is described below and should
    work with the current code in a similar way.

  If the MTC fired the correlation array is updated.  Front and back
  channels are set to one to use only one pixel of the correlation
  array and all decay-implant times are taken with respect to the MTC
  tape time.  When the MCP is being used the condition for an implant
  into the correlation array should be either the MTC or MCP.

  If the beta detectors fired a decay event is sent to the correlation
  array.

  \section Plotting
  Correlation time versus gamma-ray energy should be plotted in 
  GeProcessor::Process().

  
*/

/*********************************************************************/

/** \page RMS RMS 
  This page describes the experiment specific analysis that is
  conducted when running a particle decay spectroscopy experiment
  at the back of the RMS.  

  \section detreq Detector Types
  For an RMS experiment the following detector types and sub types
  are used.
  - dssd_front - dssd_front - the front of the DSSD.
  - dssd_back  - dssd_back  - the back of the DSSD.
  - MCP - 1position1 - the first channel of the position sensitive MCP.
  - MCP - 1position2 - the second channel of the position sensitive MCP.
  - MCP - 1position3 - the third channel of the position sensitive MCP.
  - MCP - 1position4 - the fourth channel of the position sensitive MCP.
  - MCP - 1time - the timing signal of the position sensitive MCP.
  - ge - clover - Clover detectors if used.

  \section logproc Event Processing
  The multiplicities on the front and back of the DSSD are retrieved
  as are the maximum energies deposited in both the front and back and
  strip numbers where it occured.  If both the energy deposited in the
  front and back is above a user defined threshold 
  (variable DssdProcessor::cutoffEnergy), then the event is treated as an 
  implant.  If both front and back energies are below the cutoff,
  the event is a decay.  If one is a above and the other below the event
  is called an "unknown" type.  The event is then correlated.  If only
  the front or the back fires then the energy deposited into the DSSD
  is checked against the cutoff energy but the event is not correlated.
  If the mcp fired the position is determined from the four mcp signals.

  \section Plotting
  All dssd plots are filled using DssdProcessor::Process() when either 
  the front or back fires. Ge signals are plotted in GeProcessor::Process()
  and the 2D mcp position spectra is plotted in McpProcessor::Process()

*/

/*********************************************************************/

/** \page mapfile mapfile
  New page for map description

*/

/*********************************************************************/

/** \page calfile calfile
  New page for cal description
*/

/*********************************************************************/

/** \page dettype Detector Types
  This page describes all detector types that are possible
  to use in the current analysis.  These are the only
  detector types that are valid for use in the map.txt file.

  \section dssd_front
  This detector type is used for the front of the DSSD.  The front
  and back of the DSSD must be kept separate to independently determine
  the strip with the maximum energy for front and back.  Available 
  sub types include:
  - dssd_front
  This is processed by DssdProcessor.

  \section dssd_back
  This detector type is used for the back of the DSSD.  Available
  sub types include:
  - dssd_back
  This is processed by DssdProcessor.  

  \section ge
  This detector type is used for all germanium detectors.  The
  clover detectors are the furthest developed at present.  Available
  sub types include:
  - clover
  - sega
  This is processed by GeProcessor.
  
  \section si
  This detector type is for all silicon detectors that are not a
  DSSD.  This type should be used for the silicon box. No sub types
  for this type yet.
  
  \section mcp
  This detector type is for all mcp detectors.  Available sub types
  include:
  - 1position1 - MCP position signal 1.
  - 1position2 - MCP position signal 2.
  - 1position3 - MCP position signal 3.
  - 1position4 - MCP position signal 4.
  - 1time - MCP timing signal.
   This is processed by McpProcessor.

  \section timeclass
  This detector type is for timing detectors such as TACs and MTCs. 
  Available sub types include
  - MTC (processed by MtcProcessor)
  - **

  \section scint
  This detector type is for scintillator detectors and could be
  used for NaI, beta scintillators, ...  Available sub types include:
  - beta
  - neutr
  This is processed by ScintProcessor.

  \section idssd_front
  This detector type is for the front of the dssd for implants only.
  This type is intended for use in fragmentation experiments where
  the implant and decay signals use two separate electronics chains.
  Available sub types include:
  - idssd_front

  \section position
  This detector type is for a position signal at a fragmentation experiment.
  This type is intended for correcting a time-of-flight signal based
  on the position signal at an intermediate focus.  Available sub types 
  include:
  - i2spos
  
*/

/*********************************************************************/

/** \page philosophy Overview
  The rewrite of the analysis codes was undertaken for a variety of
  reasons including:
  - Provide a consistent framework for event processing.
  - Take advantage of C++ features.
  - Enable the inclusion of ROOT as a post experimental processing tool.
  - Provide a route for adding new detector and or features that
    do not require alterations significant portions of the code.

  \section analysisover Analysis Overview
  In the Pixie16 analysis the basic unit of information is the channel.
  For each channel that triggers the acquisition the channel
  identification, energy, time, and trace (if applicable) are grouped
  together and sorted according to the times at which they occur (based
  on the Pixie16 100MHz clock).  The time sorted
  list of channels is then scanned and channels that fall within a
  certain time window are grouped together into events.  The events
  are subsequently processed in two stages.  The first stage operates channel
  by channel within the event and performs such operations as
  calibrations and raw parameter plotting and should not be altered.
  The second level of processing is the experiment specific event
  processing and includes such activities as correlations between
  implant events and decay events.  This section section should be
  altered according to experimental need. 

  A schematic data analysis flow chart is shown in the *** file
  located in the manual directory of the pixie16 analysis code.  

  \section evcreat Event Creation 
  For each buffer that is collected (from either a ldf file or an 
  online source) the SCAN routine invokes the hissub_() function.
  The hissub_(), hissub_sec(), and ReadBuffData() functions serve 
  to ultimately reconstruct a spill of Pixie16 data (i.e. a read of 
  all Pixie16 modules) and sort the channels that triggered according 
  to their times.

  ScanList() is called on the time sorted event list to group channels
  with similar times into events, store these channels in the variable
  rawevent and pass them on for further processing.

  \section experinde Experiment Independent Processing
  In the function DetectorDriver::ProcessEvent(), the RawEvent is processed. 
  Each channel is checked against its threshold value read in during
  initialization from cal.txt and calibrated using ThreshAndCal().
  The function ThreshAndCal() also keeps track of the multiplicities
  of different detector types and the maximum energy deposited into each
  detector type.  The raw and calibrated energies are then plotted using
  the functions PlotRaw() and PlotCal() using the DAMM spectra numbers
  read in from map.txt. This is the end of general event processing and
  should not be altered between experiments.

  \section experde Experiment Dependent Processing
  The next step in the event processing is the experiment specific tasks
  and these will vary for different setups but include such
  tasks as determining correlations between decays and implants.

  \section plott Plotting
  Once all event processing has finished various plotting routines are
  invoked for creating the DAMM spectra we all know and love.
  
*/

/*********************************************************************/

/** \page primer Primer
  In this section I will attempt to provide a rough description
  of the aspects of C++ that are being used in the analysis code
  that someone who has never seen C++ before may not recognize.
  Described here are
  - classes
  - pointers
  - vectors
  - and maps

  \section Classes Classes
  The C++ code is structured around classes.  Classes are very similar
  to structures in C, and are objects into which data is placed or
  actions are performed.  For example the class RawEvent is the basic
  object into which all event data will be placed in the analysis.
  Another example is the class Correlator which performs a decay-implant
  correlation using a RawEvent object.  In the analysis, all classes
  consist of two different files, a '.h' file defines the class,
  its data members, and the functions that allow a user to 
  manipulate the data members.
  The following code shows an example for the class called ChanEvent
  which contains the information for each channel that triggered.
  \code
// file RawEvent.h - class ChanEvent \\
class ChanEvent {
    private:
    double energy;             ///< Raw channel energy 
    double calEnergy;          ///< Calibrated channel energy, calibration ...
    double calTime;            ///< Calibrated time, currently unused 
    vector<double> traceInfo;  ///< Values from trace analysis functions 
    vector<int> trace;         ///< Channel trace if present 
    unsigned long trigTime;    ///< The channel trigger time, trigger time and ...
    unsigned long eventTimeLo; ///< Lower 32 bits of pixie16 event time 
    unsigned long eventTimeHi; ///< Upper 32 bits of pixie16 event time 
    unsigned long runTime0;    ///< Lower bits of run time 
    unsigned long runTime1;    ///< Upper bits of run time 
    unsigned long runTime2;    ///< Higher bits of run time 

    double time;               ///< Raw channel time, 64 bit from pixie16 channel event time 
    int    modNum;             ///< Module number 
    int    chanNum;            ///< Channel number 

    void ZeroNums(void);       ///< Zero members which do not have constructors associated with them 
    
    // make the front end responsible for reading the data able to set the channel data directly
    friend int ReadBuffData(unsigned int *, unsigned long *, vector<ChanEvent *> &);
 public:
    static const double pixieEnergyContraction; ///< energies from pixie16 are contracted by this number


    void SetEnergy(double a)    {energy = a;}    ///< Set the raw energy in case we want...
    void SetCalEnergy(double a) {calEnergy = a;} ///< Set the calibrated energy 
    void SetTime(double a)      {time = a;}      ///< Set the raw time 
    void SetCalTime(double a)   {calTime = a;}   ///< Set the calibrated time 
    void AddTraceInfo(double a) {traceInfo.push_back(a);} ///< Add one value to the traceinfo 

    double GetEnergy() const      {return energy;}      ///< Get the raw energy 
    double GetCalEnergy() const   {return calEnergy;}   ///< Get the calibrated energy 
    double GetTime() const        {return time;}        ///< Get the raw time
    double GetCalTime() const     {return calTime;}    ///< Get the calibrated time 
    const vector<int> &GetTraceRef() const {return trace;} ///< Get a reference to the trace 

    unsigned long GetTrigTime() const    
        {return trigTime;}    ///< Return the channel trigger time 
    unsigned long GetEventTimeLo() const
        {return eventTimeLo;} ///< Return the lower 32 bits of event time 
    unsigned long GetEventTimeHi() const
        {return eventTimeHi;} ///< Return the upper 32 bits of event time 
    unsigned long GetRunTime0() const
        {return runTime0;}    ///< Return the lower bits of run time 
    unsigned long GetRunTime1() const
        {return runTime1;}    ///< Return the middle bits of run time 
    unsigned long GetRunTime2() const
        {return runTime2;}    ///< Return the higher bits of run time 

    const Identifier& GetChanID() const; ///< Get the channel identifier
    int GetID() const;                   ///< Get the channel id defined as pixie module # * 16 + channel number 
    double GetTraceInfo(unsigned int a) const; ///< Get a specific value from the traceinfo

    ChanEvent();
    void ZeroVar();
};

  \endcode
  In the ChanEvent example, all the data members are defined as
  private and functions that act on the data members are public.
  This means that only the functions belonging to this class such
  as GetEnergy() can access the private data members.  This prevents
  the variables from being inadvertently altered by a different
  portion of the analysis code.  I've tried to keep all data members
  private but haven't in all cases just for simplicity (for example
  see the rawevent class in the rawevent.h file).

  The arguments that must be passed to a function are defined 
    between the parenthesis after the functions name.  Thus, the
  function SetEnergy() must be provided a double-precision numerical value
  in order to work properly and would look like the following:
  \code
  anEvent.SetEnergy(3.1416);
  \endcode 
  The other important aspect of the functions are their return values, which is
  defined before the function name.  For example, the function GetEnergy()
  returns the value contained in the variable energy which is a double.  
  Therefore to correctly retrieve the value you need to use a double.
  \code
  double tempEnergy;
  tempEnergy = anEvent.GetEnergy();
  \endcode
  Arguments to functions and results therof are not limited to basic data types 
  (int, double, char, ...) but can include other objects as well.
  In the function, GetChanID() function an Identifier object is returned.
  When using any of the functions make sure you know what type of argument
  is required. 

  Functions that are defined in the ".h" file of the class are said
  to be defined "inline".  Most of the functions in ChanEvent are 
  defined inline. One exception is ZeroVar() which is defined
  in the corresponding *.cpp file; a portion of which is shown below:
  \code
  // file rawevent.cpp - ZeroVar()
  void ChanEvent::ZeroVar(){
    ZeroNums();

    // clear objects
    trace.clear();
    traceInfo.clear();
  } 
  \endcode
  The only exception to this *.h and *.cpp pattern is PixieStd.cpp
  No class is defined for PixieStd because the main function
  hissub_ is called from the FORTRAN scan program.  A class could
  probably be constructed but I didn't want to mess with it.

  \section Pointers Pointers
  Pointers are both an incredibly powerful feature of C++ and the easiest
  way to screw up the code.  The two types of pointers in C++ are a
  pointer by value and a pointer by reference.  A pointer by value is
  declared by an '*' in the C++ code and points to the value contained
  in a variable.  A pointer by reference is
  declared by an '&' in the C++ code and points to the memory location 
  where a variable is stored. An example usage is:
  \code
  ChanEvent eventList[100];
  // ...
  ChanEvent *chan;
  chan = &eventList[1];
  \endcode
  In this section of code a pointer by value 'chan' is created that that
  will point to a ChanEvent object. The second line puts the memory
  location of eventList[1] (a ChanEvent object) into the pointer
  chan.  In this way it is possible to perform actions on chan and have
  it affect the values contained in the variable eventList[1].

  Functions also rely heavily on pointers.  By default when a C++
  function is passed a variable a local copy is created inside the
  function.  This means that if a value is passed to a function in
  this manner the function can use the value and alter it but after
  the function is finished the local copy of the value disappears, 
  and the original is not affected.  In general, it is far more useful
  to have a function act on and permanently alter the original value.
  This is where pointers are useful.  A pointer by reference or value can be
  passed to a function and actions on the pointer will affect the
  original data.  The passing of pointers is also much less computationally
  intensive speeding up the data analysis.  For example, a function
  declared to receive a pointer by reference would look like:
  \code
  int DetectorDriver::ThreshAndCal(ChanEvent *chan);
  \endcode
  Where a reference to the variable passed to the function is received.

  \section Vector Vector
  The vector class is part of the STL (standard template library) that
  is used in the analysis.  A vector is essentially a dynamically
  expandable array.  A vector variable is defined as
  \code
  vector<int> intVec;
  \endcode
  where intVec is declared as a vector that contains integer values.
  The two common functions that are used to act on vectors include the
  push_back and clear functions.
  \code
  intVec.push_back(2);
  intVec.push_back(5);
  intVec.clear();
  \endcode
  The push_back function inserts a value at the end of the vector.
  In this example the first element of the vector intVec
  will have a value of 2 and the second element is 5.  The clear()
  function removes all entries in a vector and is extensively used
  in the code to zero the various vectors that are used.  The values
  of a vector can be retrieved in two different ways.  First, you
  can ask for a specific element using the same syntax as an array, thus
  \code
  cout << intVec[1] << endl;
  \endcode
  will print out a value of 5.

  \warning C++ index numbering starts at 0 not at 1 as in FORTRAN.
  Let me repeat since this causes many problems if forgotten.
  C++ index numbering starts at 0 not at 1 as in FORTRAN.

  The second way to retrieve a value from a vector is through an iterator.
  This is a pointer that can be pointed any location in the array to
  retrieve the value at that location.  The following section of code
  would loop over the intVec from beginning to end and print the value
  of the element.  Note the incrementing of the iterator iv using iv++
  to go to the next element.
  \code
  vector<int>::iterator iv;
  iv = intVec.begin()
  while(iv != intVec.end())
  {
     cout << *iv << endl;
     iv++;
  }
  \endcode

  \section Map Map
  The map is the other STL feature used in the code.  The map
  allows an association between an unique key and a value.
  It is defined as
  \code
  map<string,int> stringIntMap;
  \endcode
  In this example the key is a string and the associated value is
  an int.  To create a map you have to make pairs of the keys and values
  as in the following example.
  \code
  stringIntMap.insert(make_pair("dssd",2));
  stringIntMap.insert(make_pair("ge",15));
  \endcode
  In this example the key "dssd" is associated with an integer of 2
  and the key "ge" is associated with an integer of 15.  The
  associated values can be retrieved either by a direct access
  using the key or using of the find function.
  \code
  int temp;
  temp = stringIntMap["dssd"];
  
  int temp1;
  temp1 = stringIntMap.find("ge");
  \endcode
  where temp and temp1 will return 2 and 15 respectively.  While 
  associating strings with integers is not very useful, a map of
  strings associated with objects is a powerful method to associate
  an unknown number of detector types with data objects as is done
  in the pixie16 analysis.
  
  \warning   I (SNL) am not a C++ expert! It is possible that in the analysis code
  things are done inefficiently.  Also, the descriptions I have given
  above may not be completely accurate but instead reflect a lack of
  understanding on my part.  Keep that in mind when reading the manual
  and the code.  If you think things could be done differently for
  clarity or processing speed please let me know.

  \warning I (DTM) might be a C++ expert.
*/

/*********************************************************************/

/** \page directory Directory Structure

  The pixie16 analysis is contained in the directory scan_c++. You
  should see at a minimum the following set of files and directories.
  \verbatim
  cal.txt  Doxyfile  html  include  latex  Makefile  manual  map.txt  scan  src
  \endverbatim

  The files map.txt and cal.txt are used to initialize the pixie16
  analysis and are described in the next section.  The directories src/ and
  include/ contain the *.cpp and *.h files for the various classes
  that are used in the analysis.
  Scan files such as drrsub.f are in the scan/ directory. 
  The html/ and latex/ directories contain the html and latex
  ouput produced by doxygen.  Links to the pixie16manual in both pdf and 
  html format are in the manual directory. 
  The Doxyfile is the configuration file used
  with doxygen to create this manual.  
*/

/*********************************************************************/

/** \page quickstart Quickstart

  \section Overview
  This quickstart guide will quickly describe how to use the current pixie16
  analysis.  The guide assumes you have no desire to know how the
  analysis performs various tasks but only want to know how to make
  the program run and do something useful.  If you are unfamiliar with
  C++ please refer to the C++ primer for a description of a few
  C++ features used in the code that have no FORTRAN equivalents.
  This chapter will use a sample ldf file taken from an RMS experiment
  with a 40x40 DSSD and one MCP.  This tutorial assumes you have a
  working knowledge of the DAMM program for viewing *.his files.

  \section Complilation
  I assume that ROOT has not been installed on your computer and you do not
  want to compile with ROOT capabilities.  In this case make sure that
  declarations of USEROOT and STANDALONEROOT are commented out in 
  the Makefile.
  \verbatim
  # use a "#" on the beginning of a line to comment it out of the Makefile
  #USEROOT = 1
  #STANDALONEROOT = 1
  \endverbatim
  
  The next flag that must be set correctly is whether your data was
  taken with the new pixie16 readout.  This readout was developed
  for the high data rates at the back of the RMS.  If you have DSSD
  RMS data taken after March 20th 2008 then it is likely you have used the
  new pixie16 readout scheme.  If so, then make sure the -Dnewreadout
  flag is set in the following line of the Makefile.
  \code
  G++FLAGS  = -O3 -Wall -Wno-deprecated -fPIC $(CINCLUDEDIRS) -Dnewreadout
  \endcode
  If you have older DSSD RMS data or beta decay spectroscopy data from
  LeRIBSS (or similar setup) then you have most likely not used the 
  new readout scheme and the -Dnewreadout flag should be commented out.
   \code
  G++FLAGS  = -O3 -Wall -Wno-deprecated -fPIC $(CINCLUDEDIRS)# -Dnewreadout
  \endcode

  Also, make sure that the appropriate compilers are inserted into
  the Makefile.  The Makefile assumes that g++ (C++), gcc (C), and
  g77 (FORTRAN) compilers are present on your computer.  If your
  compilers have different names please alter the definitions inside
  the Makefile (at the location shown below) to match your system.
  \verbatim
  #--------- define compilers
  G77 = g77
  GCC = gcc
  G++ = g++
  \endverbatim
  Additionally, you must point the make file to the appropriate
  location where you have installed the static scanor libraries.
  The ones needed are scanorlib.a, orphlib.a, acqlib.a, and ipclib.a.
  Alter locations of these static libraries by changing the definitions
  of DIRA2 and DIRB in the Makefile
  \verbatim
  DIRA2=/usr/hhirf/g77
  DIRB= /usr/acq2/lib

  LIBS = $(DIRA2)/scanorlib.a \
         $(DIRA2)/orphlib.a \
         $(DIRB)/acqlib.a  $(DIRB)/ipclib.a
  \endverbatim

  If you do not have these libraries then please visit the HRIBF
  computing web page to download the appropriate packages.

  Compile the code with the command
  \verbatim
  make
  \endverbatim

  \section inputprep Input Preparation
  In order to successfully run the pixie16 analysis two initialization
  files are required, map.txt and cal.txt.  map.txt describes the
  mapping of pixie16 input channels to specific detectors and cal.txt
  sets the calibration for each detector.  Both files are included in
  the tutorial directory and should be copied to the main scan_c++
  directory to run this tutorial.  Both files are described in
  more detail below.

  \subsection maptxt Map.txt
  The file map.txt contains a total of 7 parameters for each pixie16
  channel.  The 7 parameters are necessary for completely specifying
  the detector plugged into that channel.  The 7 parameters are:
  - module number
  - channel number
  - damm id - this is the spectrum number that the raw energies from this channel will be plotted. 
  - detector type - this is used to logically group detectors.  For example, all front channels on the DSSD have a detector type of dssd_front.
  - detector subtype - this is useful in the program if you need the data from a specific channel you can address it by name.  For example, the first channel of a position sensitive MCP would be 1position1.
  - detector location - this is the physical location of the detector.  For the dssd this is the strip number.
  - trace - this variable is obsolete.

  Comments in the map.txt file are always preceded by "//".  An example
  of the map.txt file for use at the RMS with a DSSD and one MCP would
  look like this:
  \verbatim
  //example map.txt file
  // mod id  chan id  damm    det_type   det_subtype  det_loc  trace
          0        0   600  dssd_front    dssd_front        0      1
          0        1   601  dssd_front    dssd_front        1      1
          0        2   602  dssd_front    dssd_front        2      1
          0        3   603  dssd_front    dssd_front        3      1
          ... 
          4        0   900   dssd_back     dssd_back        0      1
          4        1   901   dssd_back     dssd_back        1      1
          ... 
          5        0   927         mcp    1position1        0      1
          5        1   927         mcp    1position2        0      1
          5        2   927         mcp    1position3        0      1
          5        3   927         mcp    1position4        0      1
          ...
  \endverbatim
  If a channel is not used in the experiment or is not required in the
  analysis the detector type and detector subtype should be set to
  "ignore" in the map.txt line corresponding to that module and
  channel number.  To ignore the 0 and 1 channel on module 2 the
  map.txt file would look like:
  \verbatim
  //example map.txt file ignoring module 2 channels 0 and 1
  //mod id   chan id  damm    det_type   det_subtype  det_loc  trace
         0         0   600  dssd_front    dssd_front        0      1
         0         1   601  dssd_front    dssd_front        1      1
         0         2   602  dssd_front    dssd_front        2      1
         0         3   603  dssd_front    dssd_front        3      1
         ...
         2         0   900  ignore            ignore        0      1
         2         1   901  ignore            ignore        1      1
         ...
  \endverbatim
  
  \subsection caltxt Cal.txt
  The other necessary file for running the program is cal.txt.  As the
  name suggests this file contains the calibration parameters for the
  detectors used in the analysis. For each channel that was defined in
  map.txt that is not set to "ignore" a line must be present in cal.txt
  and include:
  - detector location - this is the physical location of the detector.
  For the dssd this is the strip number
  - detector type - this is used to logically group detectors.  For example,
  all front channels on the DSSD have a detector type of dssd_front.
  - detector subtype - this is useful in the program if you need the data
  from a specific channel you can address it by name.  For example, the
  first channel of a position sensitive MCP would be 1position1.
  - num cal - this is the number of calibrations that will be used for this
  particular channel and can be any positive number.
  - poly order - this is the polynomial order of all calibrations for this
  channel.  1 is linear, 2 is quadratic ...
  - a set of (1+ polynomial order +1) values - the first value of this set
  is the lower limit of the current calibration.  The remaining "polynomial
  order +1" parameters are the coefficients for the  polynomial equation
  from lowest order to highest.  For example, a calibration with one 2nd order
  polynomial would have 4 values.  A calibration line with two 4th order
  polynomials would have two sets of 6 calibration parameters. 

  There is currently no default calibration parameters for detectors.
  Only channels set to "ignore" have default parameters.  The upper
  limit of the calibration (exclusive) is either (1) the lower limit
  of the next calibration or (2) the maximum value set by MAX_PAR in the
  'param.h' file.  Comments in the cal.txt file are always preceded by "//".
  An example of the cal.txt file for one 2nd order polynomial looks like this:
  \verbatim
  //example cal.txt file for 1 - 2nd order polynomial
  //det    det    det   num   poly    low   cal val  cal val  cal val
  //loc   type    sub   cal  order  range     inter    slope     quad
      0   dssd   dssd     1      2    0.0     0.0        1.0      0.0
         front  front
      1   dssd   dssd     1      2    0.0     0.0        1.0      0.0
         front  front
  ...
  \endverbatim

  An example of a calibration file for two 2nd order polynomials is
  shown below.  The first calibration ranges from 0.0 to 16000.
  The second calibration has a range from 16000 to the maximum set by
  MAX_PAR in the 'param.h' file.
  \verbatim
  //example cal.txt file for 2 - 2nd order polynomials
  //det    det    det   num   poly    low    cal   cal   cal   low    cal   cal   cal
  //loc   type    sub   cal  order  range  inter slope  quad        inter slope  quad
      0   dssd   dssd     2      2    0.0    0.0   1.0   0.0 16000    0.0   1.0   0.0
         front  front  
      1   dssd   dssd     2      2    0.0    0.0   1.0   0.0 16000    0.0   1.0   0.0
         front  front
  ...
  \endverbatim

  \section eventcreat Event Creation
  Events are created by grouping together channels which triggered at
  similar times.  This time window is controlled by the variable "eventWidth"
  in PixieStd.cpp in ScanList(). A list of channels sorted by time from 
  lowest to highest is scanned in this function from beginning to end. 
   If two successive channels in the list are within "eventWidth" time,
   they are grouped together in an event.  The variable is in units of
   pixie16 clock ticks and must be multiplied by 10 to obtain the
   event width in ns.  Please note that the "eventWidth" time
   window is only applied to successive events.  Thus it is possible
   (depending on the total trigger rate) to have events that are longer
   than the specified time window.

  \section eventret Retrieving Information from the Event
  After an event has been created it is sent for processing.  The first
  step of processing is to calibrate the individual channels and summarize
  the information contained in the event.  For each detector type that is
  present in the event an object called DetectorSummary is created.
  This object holds detector related information such as the location of
  the maximum energy deposition and the multiplicity of the detector among 
  other things.  For example, the following command will retrieve
  the maximum energy that was deposited in the front of the DSSD in this 
  event.
  \code
  int maxe = revt.GetSummary("dssd_front")->GetMaxEvent()->GetEnergy();
  \endcode
  where revt is the name of the variable holding the raw event information
  and maxe will contain the maximum energy deposited in the front of the
  detector.  Similarly, the strip number where the maximum energy was
  deposited is retrieved with the following command: 
  \code
  int stripno = revt.GetSummary("dssd_front")->GetMaxEvent()->GetLocation();
  \endcode
  Therefore a plot of strip number versus energy for the front of the DSSD
  can be obtained by plotting the variable stripno against maxe (using the
  syntax described later in the plotting section).
  
  To retrieve the multiplicity associated with the front of the DSSD
  the following command is used
  \code
  int fmult = revt.GetSummary("dssd_front")->GetMult();
  \endcode

  If you already have a pointer to the detector summary of interest then
  the multiplicity could be retrieved using the pointer itself as follows
  \code
  int fmult = pdet_s->GetMult();   or
  int fmult = pdet_s.GetMult();
  \endcode
  where pdet_s is a pointer to the detector_summary and the choice between
  a "->" or "." depends on whether pdet_s is a pointer by value or
  pointer by reference, respectively.  The reference manual can provide 
  a list of all commands to retrieve information from the detector_summary
  or the rawevent.

  \section Plotting Plotting
  All plotting is controlled through the "plot" function defined in 
  DeclareHistogram.cpp.  This function is a C++ wrapper that has been
  placed around the DAMM count1cc and set2cc subroutines.  This allows for
  the code to be easily changed between damm and ROOT histograms.  For
  those using DAMM to view the output of the analysis all plots are created
  in the "drrsub.f" file located in the scan directory.
  For example, to plot a one dimensional histogram of energies into a DAMM
  spectrum numbered 101 the command would be 
  \code
  plot(101,energy);
  \endcode
  You can send any type of numerical value to the plot function.  The
  variable is cast into an integer before being passed to the DAMM
  plotting functions.
  A two dimensional histogram is plotted with the command
  \code
  plot(101,energy,stripno);
  \endcode
  and a three dimensional histogram (plotting a trace for example) uses
  the command
  \code
  plot(101,xval,row,traceval);
  \endcode

  \section programexec Program Execution
  After compilation the executable pixie_ldf_c will be present.  Run
  the pixie_ldf_c program as you would any other scan program as follows
  \verbatim
  ./pixie_ldf_c hisname
  \endverbatim
  Where "hisname" is the name of the damm histogram that will be created.

  For this exercise you will use the tutorial.ldf file located in the
  LDF directory. At the scanor prompt load in the appropriate ldf file
  as follows
  \verbatim
  SCANOR->file tutorial/tutorial.ldf
  \endverbatim
  Next start the analysis with the following command
  \verbatim
  SCANOR->go
  \endverbatim

  After starting a variety of output should be printed to the screen.
  The first lists all the detectors that are being used in the analysis.
  For the tutorial file the dssd_front, dssd_back, and mcp detector
  types are used in
  the analysis.  Next is a list of all module and channel number
  combinations and their respective calibrations.  As a sanity check the
  total number of channels that the code has identified is given and should
  be equal to the number of %calibration lines.  Lastly, when
  new detectors are added to the code the number of possible detectors is
  given from two places in the code and the numbers should be identical.

  After completion of the analysis, end the SCANOR program
  \verbatim
  SCANOR->end
  \endverbatim

  Start damm and open your histogram.  You can display 741, 742, 743, and
  744 to see the strip versus energy profile of the dssd for implants and
  decays on the front and back of the dssd.  Spectrum 923 provides a 2D plot
  of the position sensitive MCP.    

  These spectra are filled following the experiment specific event processing
  in the DetectorDriver (defined in DetectorDriver.cpp) by invoking
  the appropriate plotting routines. The 741, 742, 743, and 744 spectra are
  all filled using the DssdProcessor object and mcp 923 spectrum is filled
  by the McpProcessor object.  The plotting routines of these objects
  are invoked from the DetectorDriver in each EventProcessor's
  individual Process() functions.

  To create a decay multiplicity spectrum for the front of the DSSD you need
  to first create the spectrum and fill it in DssdProcessor.cpp.
  Create the spectrum in DeclarePlots() by inserting the following line
  \code
  DeclareHistogram1D(1743, 64, "DSSD decay front mult", 2, 64, 0, 63);
  \endcode
  or more simply (using assumed parameters):
  \code
  DeclareHistogram1D(1743, 64, "DSSD decay front mult");
  \endcode
  Fill the decay front multiplicity spectrum by going inside the the
  "if(type == decay)" and "if(side == front)" condition and inserting
  the following line after the filling of spectrum 743.
  \code
  plot(1743,frontSummary->GetMult());
  \endcode

  Save the two files and recompile the pixie16 scan program and rerun the
  analysis on tutorial.ldf.  Open the resulting damm file and you should
  have a new multiplicity spectrum.

  Next create a decay multiplicity 2 gated version of spectrum 743.  First
  create the new spectrum in DeclarePlots() by inserting the following line
  \code
  DeclareHistogram2D(2743, 8192, 64, "DSSD Strip v E-DF, mult2",
                     2, 8192, 0, 8191, 64, 0, 63);
  \endcode
  or once again:
  \code
  DeclareHistogram2D(2743, 8192, 64, "DSSD Strip v E-DF, mult2");
  \endcode

  Next insert the filling command inside the ''if(type == decay)`` and
  ''if(side == fbcoin)'' condition.
  \code
  if(frontSummary->GetMult()==2)
    plot(2743,frontSummary->GetMaxEvent()->GetEnergy(),
              frontSummary->GetMaxEvent()->GetLocation());
  \endcode

  Save, recompile and rerun the analysis and you should now have a
  multiplicity 2 gated 743 spectrum in spectrum 2743. 

  \section commonerr Common Errors
  To familiarize yourself with some of the more important workings of
  the code we will now perform a few common tasks that will result in
  the failure of the analysis code.

  \subsection mapprob Problems with map.txt
  The easiest mistake is misspelling the detector name.  Open the file
  map.txt and change the entry for module 0 channel 15 to the following
  \verbatim
  //example incorrect map.txt file
  // mod id  chan id  damm    det_type   det_subtype  det_loc  trace
          ...
          0       15   616    dssd_fro    dssd_front       16      1
          ...
  \endverbatim
  This type of error will be caught by the scan program and the following
  warning will be issued.
  \verbatim
  The detector called ``dssd_fro'' read in
  from the file map.txt is unknown to this program!
  This is a fatal error.  Program execution halted.
  If you believe this detector should exist please
  edit the 'GetKnownDetectors' function inside the
  'DetectorDriver.cpp' file.

  The currently known detectors include:
  dssd_front, dssd_back, ge, timeclass, si,
  position, ifront_dssd, scint, mcp
  \endverbatim

  This warning states that you are trying to use a detector (dssd_fro) for
  which no default behavior has been assigned.  This error is the likely
  result of a misspelling and the valid detector names are listed as well.
  If you need to add a new detector type please refer to ***.  Note: the
  program does not check to make sure that the detector subtype is spelled
  correctly.  Therefore if you have code that makes use of this information
  (the MCP is one example) then make sure the spelling is consistent 
  between map.txt and the analysis code.

  \subsection calprob Problems with cal.txt
  
  The next common error is not having the correct number of lines in the
  calibration file.  Open cal.txt and change the calibration line for
  dssd_front strip number 16 to the following
  \verbatim
  //example cal.txt file for 1 - 2nd order polynomial
  //det    det    det   num   poly    low   cal val  cal val  cal val
  //loc   type    sub   cal  order  range     inter    slope     quad
     15   dssd   dssd     1      2    0.0     0.0        1.0      0.0
         front  front
  \endverbatim
  Now you should have two lines referring to dssd_front strip number 15.
  Trying to run with this cal.txt produces the following error.
  \verbatim
  The entry in map.txt with det_type = dssd_front
  det_subtype = dssd_front
  and det_location = 16
  was not found in the cal.txt file, please fix this
  Program terminating
  \endverbatim
  All entries in the map.txt file should have a corresponding entry in the
  cal.txt file.
    
*/

/*********************************************************************/

/** \page faq FAQ

  This FAQ will be filled based on suggestions accumulated during
  the manual beta phase.

*/

/*********************************************************************/

/** \page codewarn Warnings

  This section provides a list of warnings that appear throughout
  both the users manual and reference manual.  If you can not 
  determine why something has either failed or produced incorrect
  results look at the list of warnings to see if it has been mentioned.

  \warning All channels defined in map.txt must have a corresponding
  calibration.
  \warning C++ index numbering starts at 0 not at 1 as in FORTRAN.
  Let me repeat since this causes many problems if forgotten.
  C++ indexing numbering starts at 0 not at 1 as in FORTRAN.
  \warning If you are not using ROOT then the flags USEROOT and
  STANDALONEROOT should be commented out in the Makefile.
  \warning If you have RMS DSSD data taken after March 20th 2008 then
  most likely the new pixie16 readout was used and the flag -Dnewreadout
  should be set in the Makefile.

*/

/*********************************************************************/

/** \page objects Objects
  This page will describe the various objects that are used in the
  analysis.  These objects are both for detectors and processing.
  All objects are declared, initialized and used in the
  DetectorDriver class.

  \section detobj Detector Objects

  \subsection dssd DssdProcessor
  The DssdProcessor object is intended for use with the two DSSD detector
  types dssd_front and dssd_back.  The DssdProcessor object is currently
  only used to plot DSSD spectra including the DSSD strip versus energy 
  spectra (741, 742, 743, 742) and DSSD hit pattern spectra (725, 726). 
  Both dssd_front and dssd_back detector types feed into this object for
  plotting.

  \subsection event EventProcessor
  The EventProcessor is the generic processor from which all the others are
  derived. This allows the common process timing and setting up of the
  associated types and DetectorSummaries in a common framework.

  \subsection gep GeProcessor
  The GeProcessor object is intended for use with any germanium detector.  
  All detectors with a detector type of Ge feed into this object.  The
  GeProcessor object is developed for clover detectors and its two main
  functions plot germanium spectra and perform addback.

  \subsection mcpp McpProcessor
  The McpProcessor object is intended for use with all mcp detectors.  All
  detectors with a detector type of mcp feed into this object.
  The McpProcessor object calculates the 2D position of an event based on 
  the collected charge from the four corners of the detector.

  \subsection mtc MtcProcessor
  The MtcProcessor object is intended for analysis with detectors
  of type timeclass and subtype mtc (alternatively type mtc will work).
  No correlation of implants and decays is currently done, however this
  could be easily implemented through the Correlator.

  \subsection root RootProcessor
  The RootProcessor object is designed to output a ROOT tree for further 
  analysis. In particular, it calls the AddBranch() and FillTree() methods
  of all the EventProcessors included in the present detector driver. Note
  that this requires compilation using the USEROOT flag. Also this is a
  significant departure from the previous ROOT philosophy. In particular,
  the data is separated by its true function in the experiment and significant
  pre-processing of the data is done into what we would want to access about
  the data. Raw channel events are not dumped, but this behaviour probably
  will be again implemented soon.

  \subsection scint ScintProcessor
  The ScintProcessor object is intended for use with an scintillator object 
  and is fed from any detector type labeled as scint. This class was
  particularly discombobulated in the previous version of the code. At the 
  moment, it only plots things for scintillators with subclass neutr.

  \subsection spline SplineProcessor
  \warning In progress
  
  The SplineProcessor is designed to derive sub-sample timing resolution
  based on the pulse shape of the signal. Several methods are being
  investigated and this class will develop accordingly.

  \section proobj Event Processing Objects

  \subsection detdriver DetectorDriver
  The DetectorDriver object controls event processing.  A raw event
  is passed to it from ScanList() and all channels in the raw 
  event are checked against their threshold and calibrated. 
  Then the experimental processing is started.

  \subsection correlator Correlator
  The Correlator object controls the correlation of decays with previous
  implants.  Two 2D arrays of size MAX_STRIP by MAX_STRIP are created
  for both implants and decays.
  Decays are correlated to implants if the time between implants was
  sufficiently large and the time between decays and implants is less
  than the correlation time.  This object also controls the plotting
  of time versus energy plots such as the 750s.

  \subsection trace TraceAnalyzer
  The TraceAnalyzer object implements a quick online trace analysis 
  routine based on trapezoidal filters to identify double pulses.
  The object also plots the trace spectra.

  \section other Other objects
  \subsection stats StatsAccumulator
  This is a utility object for finding means and standard deviations
  from an STL object using the STL accumulate() algorithm.
*/

/*********************************************************************/

/** \page flags flags
  There are currently only three compilation flags that are important
  for the pixie16 analysis.  These include:
  - newreadout - This flag signifies that the new readout scheme is
  being use to transmit data from the pixie16 modules across the
  network.  On any experiment after Apr 2008 this is the most likely
  running mode so compile with this option included.  Without this 
  option the code will compile into the old readout scheme for 
  backwards compatability.
  - USEROOT - To use this flag, you must have installed the ROOT 
  analysis package on your system and have the ROOTSYS, PATH,
  and LD_LIBRARY_PATH environment variables set correctly (according
  to the ROOT installation). This  flag compiles in the ability to
  output a .root file depending on a set of conditions in the
  RootProcessor but otherwise does not affect the functioning of the code.
  - REVD - This sets the scan to work on revision D modules, this only
  effects ReadBuffData() and some functions in PixieStd.cpp

  One flag present in the old analysis will likely be implemented in the future:

  - STANDALONEROOT - The ROOT analysis package mus be installed 
  correctly to use this option.  With this flag the program no 
  longer runs with scan and damm but instead takes an input
  ROOT file (likely created with the USEROOT option above) and produces
  another output .root file.  The STANDALONEROOT flag must be used
  with the USEROOT flag.


*/

/*********************************************************************/

/** \page add add a detector
  This page will describe the mechanism for adding a new detector
  to the analysis.  This should only be done rarely as many detectors
  should fit into one of the already defined \ref dettype.

  To add a new detector first a class must be constructed which
  defines its behavior by creating a .h and .cpp file. For example,
  the following section defines a new detector class called "new".
  \code
  //file NewProcessor.h - adding a new detector to analysis

  #ifndef __NEWPROCESSOR_H_
  #define __NEWPROCESSOR_H_

  #include "EventProcessor.h"

  // add detector
  class NewProcessor : public EventProcessor 
  {
  private:
     //list of private variables
  public:
     NewProcessor();
     virtual void Process(RawEvent &event);
     virtual void DeclarePlots(void) const;
     
     // if you want to add root analysis...
#ifdef USEROOT
     virtual bool AddBranch(TTree *tree);
     virtual void FillBranch(void);
#endif
  };

  #endif // __NEWPROCESSOR_H_

  //file NewProcessor.cpp - adding a new detector to analysis
  NewProcessor::NewProcessor() : EventProcessor()
  {
     name = "new";
     associatedTypes.insert("new");
     //add any special instructions for object creation
  }
  
  void NewProcessor::DeclarePlots() const
  {
     // add DeclareHistogram?D statements
     // related descriptions of the ID numbers go in damm_plotids.h
  }

  bool NewProcessor::Process(RawEvent &evt)
  {
     // the following lines are necessary for the timing of the process
     //   to be implemented, it also ignores the process if none of the
     //   associated types are part of this event
     if (!EventProcessor::Process(evt))
       return false;
     //define plotting routines.

     EndProcess(); //update the processing time
     return true;
  }

  #ifdef USEROOT
  bool NewProcessor::AddBranch(TTree *tree) 
  {
    // create the branch in the tree
    return true;
  }

  void NewProcessor::FillTree(void)
  {
    // update the variables associated with this processor's branch
    // also, since ROOT trees don't allow unfilled values,
    //  one should clear the values for detectors with no associated
    //  events which is often all that is needed since much of the
    //  relevant processing appears in Process() 
  }
  #endif // USEROOT
  \endcode

  Once a class exists for your new detector you must create an
  instance of it in the DetectorDriver.  Include the proper file in 
  DetectorDriver.cpp
  \code
  #include "NewProcessor.h"
  \endcode
  and initialize it with the other processors in the constructor of 
  DetectorDriver.
  \code
  vecProcess.push_back(new NewProcessor());
  \endcode
  
  Now the only thing that remains is to add it to the list of known detectors
  in the GetKnownDetectors() method, and all the rest should be taken care 
  of for you.
*/

/** \page changes Changes
    
    \section general General changes
    - Consistency in the naming conventions of variables/classes included in the code
    - Names made more descriptive when necessary
    - Prevent as much copying of information between classes for optimization
    - Old ROOT-style functionality currently unimplemented
    - C++ streams used instead of stdio functions as much as possible
    - Comments improved -- simple coding statements within comments removed
    - Unused variables removed, or commented out when the extension of the program to include their eventual use is evident
    - Classes are more const correct and therefore users shouldn't be able to mess up things that they shouldn't touch

    \section correlator Correlator.cpp
    - Previously called \c correlator.cpp
    - Header now includes correlation structures (previously in \c raw_event.cpp) with different structures for implants and decays
    - Now only one main correlation function exists
    - DSSD decay spectra plotting now located in DssdProcessor.cpp

    \section histo DeclareHistogram.cpp
    - New set of global functions which take over the behaviour of \c drrsub.f
    - Main plotting functions and generic low-level histogram definitions defined here
    - All used DAMM plot ids defined in damm_plotids.h
    - No DSSD strip spectra automatically defined

    \section driver DetectorDriver.cpp
    - \b Major changes
    - Previously called \c detector_driver.cpp
    - Calculation of calibrated energies now in Calibration class
    - Significant reduction of members of DetectorDriver
    - Simple sanity check function to be run at beginning of analysis
    - Most of analysis is now done through a vector of event processors
    - Raw event is now accessed globally
    - Calibrated spectra are now defaulted to damm id# \f$200+ch\f$
    - Calbrations are now matched to the channel identifiers so only one consistent reading of \c map.txt is done
    - Random number initializations now done here and are scaled appropriately by the contraction factor of the Pixie energies
    - A local copy of the analysis' used detectors is no longer made
    - Detector summaries have been moved to the RawEvent class
    - No "temppoint" handling of detector summaries, instead each EventProcessor will keep its own individual map of detector summaries for the detectors it is interested in

    \section dssd DssdProcessor.cpp
    - Previously called \c dssd_sub.cpp
    - Handles much of the correlation work that was previously done by the DetectorDriver
    - Decays required to be anti-coincident with MCP signals
    - Plotting for decay spectra now done here

    \section event EventProcessor.cpp
    - New class which defines the generic behaviour for the processor which handle events for detectors of different types
    - Includes an interface to ROOT for each individual processor
    
    \section gex GeProcessor.cpp
    - Previously called \c ge.cpp
    - Number of clovers is now determined using the global vector of identifiers not from reading \c map.txt

    \section make Makefile
    - Implicit rules used for most targets
    - REVISIOND flag added to support Pixie Rev. D modules
    - Target "dist" and "distdocs" added to make tarball of code

    \section mcpx McpProcessor.cpp
    - Previously called \c mcp.cpp
    - Now includes a simple branch for "new"-style ROOT analysis

    \section mtc MtcProcessor.cpp
    - Previously called \c timeclass.cpp
    - Handling of old TOF timeclass for NSCL experiments not currently implemented

    \section pixiestd PixieStd.cpp
    - \b Major changes
    - Previously called \c pixie_std.cpp
    - Chunks starting with a -1 delimiter are now skipped over while looking for the beginning of a spill
    - If only the last 5-word chunk is missing from the spill, it is reconstructed
    - The last 5-word chunk is checked to verify that it is the last chunk expected for the spill
    - Real, system, and user time in analysis reported occasionally
    - Random number initialization moved to DetectorDriver.cpp
    - Modules are only expected to be read out cyclically
    - Support up to 14-slot crates
    - Stand-alone ROOT support currently not implemented
    - Events now placed into the list of events using their pointers, so they do not need to be copied
    - Channel identification index numbers are now extracted from the module and channel number in \c map.txt not according to the order of appearance in the file
    - Multiple definitions of the same channel in \c map.txt will now produce an error
    - Platform independce now more obvious through definition of word_t

    \section raw RawEvent.cpp
    - \b Major changes
    - Previously called \c rawevent.cpp
    - Structures for correlation now connected with Correlator.cpp
    - Actual Correlator class now included in the RawEvent
    - Channel IDs are now returned directly from the data imported from \c map.txt -- individual channel events do not contain their own copies
    - ID numbers are calculated (not stored) using module and channel numbers
    - Much less "structification" for channel events -- ReadBuffData() is a friend to set information. In particular, getting a copy of the trace is no longer supported as this can create a large burden on the analysis
    - Detector summaries are now stored with the RawEvent instead of in the DetectorDriver
    - Detector summaries now address their constituent events by pointer (not by arbitrary index in RawEvent)
    - One function call to handle all details regarding associating a channel event to a particular DetectorDriver
    - Less "structification" for the RawEvent -- access to detector summaries is only allowed by pointers to prevent copying of large amounts of channel data
    - Asking for an unexpected DetectorDriver will complain accordingly

    \section buf ReadBuffData.cpp
    - A version for both Rev. A and Rev. D of the Pixie modules
    - Main function now directly takes a pointer to the beginning of the data to be processed of Pixie words of the appropriate size
    - Identifiers are no longer set since they are calculated dynamically in RawEvent.cpp

    \section root RootProcessor.cpp
    - New tact of creating ROOT trees where the data is broken into its constituent detector parts (saved as branches) by the event processors. This is in contrast to the old method of saving raw event data

    \section scint ScintProcessor.cpp
    - Previously called \c scint.cpp
    - If no Ge detectors are present, processor asks for a NULL detectors summary which has to be treated properly
    - No more triple counting for plots
    - Discombobulated in general

    \section spline SplineProcessor.cpp
    - \b "In progress"
    - Previously called \c pulse_spline.cpp
    - Little to no modification yet as this is still under active development

    \section stats StatsAccumulator.cpp
    - New simple class which provides a consistent and concise way of calculating means and standard deviations
    - Designed to work with the STL accumulate algorithm defined in \c numeric
    - Example of use in TraceAnalyzer.cpp

    \section trace TraceAnalyzer.cpp
    - Previously called \c trace.cpp
    - Better implementation of STL methods
    - Filter functions now includes a version to fill the vector in place so an unnecessary copy is not made
    - Processing time information is included
 */

/** \page improvement Improvement
 - More dynamic system for DAMM id definition which separates modular behaviour such that all plots from a particular processor explicitly start with the same two digits; this makes inheritance easier and frees up damm_plotids.h
 - Generic processor and detector type on which the calibrator and trace analyzer can function
 - Optimization and clarification of TraceAnalyzer.cpp
 - Implementation of old standalone root analysis
 - Include QDC and CFD information from channel headers for Revision D
 - Gating class which can be applied to event processors
 - Reading of \c map.txt prior to first buffer so that the correct number of channels and event processor plots can be included in the DAMM file without any excess
 - Expand detector summaries to work on type/subtype combinations
*/

/** \page names Naming conventions
 Some sort of convention -- maybe not the best, but at least it's a choice. Of course, these are solely a suggestion, and exceptions could occur.
 - Classes have names like "ExampleClass"
 - Variables and class members have names like "myVariable"
 - Namespaces like "thisNamespace"
 - Acronyms uncapitalized generally like "DssdProcessor" (exceptions may abound)
 - Global CONSTANT variables can have forms like "CONSTANT_VARIABLE"
 - C++ functions declared as "DoSomething(...)"
 - Fortran functions called as "func_()"
 - Implementation of classes in files with the same name as the class
 - Files with extensions ".cpp", ".h", and ".f"
 - Generic index variables use "i,j,k" or more explicit names such as "mod, ch"; if you need more, consider reorganizing your code
 - Enumeration types like "EColors" with members "BLUE, RED, ORANGE"
 - Typedefs with simple names like "word_t" 
 - Spectrum IDs use namespaces as appropriate and constant names like "D_ENERGY" (1D) or "DD_ENERGY__TIME" (2D, time (y) vs. energy (x)): subject to further consideration
 - Consider using "//!" for comments which should be addressed and "//?" for ideas on how things might be improved
 */
 
