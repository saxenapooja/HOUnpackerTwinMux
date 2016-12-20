/* -*- C++ -*- */
#ifndef HOTPUnpacker_h_included
#define HOTPUnpacker_h_included 1

#include "DataFormats/HcalDigi/interface/HcalCalibDataFrame.h"
#include "DataFormats/HcalDigi/interface/HOTriggerPrimitiveDigi.h"
#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "CondFormats/HcalObjects/interface/HcalElectronicsMap.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveSample.h"
#include "DataFormats/HcalDigi/interface/HcalUnpackerReport.h"
#include <set>

class HOTPUnpacker {
public:

  struct HOCol {
    HOCol();
    std::vector<int>*     etaCont;
    std::vector<int>*     phiCont;
    std::vector<int>*   sampleCont;
    std::vector<int>*      soiCont;
    std::vector<int>* databitsCont;
    
  };
  
  /// for normal data
  HOTPUnpacker(int sourceIdOffset, int beg, int end) : sourceIdOffset_(sourceIdOffset), startSample_(beg), endSample_(end), expectedOrbitMessageTime_(-1), mode_(0) { }
  /// For histograms, no begin and end
  
  HOTPUnpacker(int sourceIdOffset) : sourceIdOffset_(sourceIdOffset), startSample_(-1), endSample_(-1),  expectedOrbitMessageTime_(-1), mode_(0) { }
  void setExpectedOrbitMessageTime(int time) { expectedOrbitMessageTime_=time; }
  void unpack(const FEDRawData& raw, const HcalElectronicsMap& emap, HOCol& conts, HcalUnpackerReport& report, bool silent=false);
  void setMode(int mode) { mode_=mode; }
  
private:
  int sourceIdOffset_; ///< number to subtract from the source id to get the dcc id
  int startSample_; ///< first sample from fed raw data to copy 
  int endSample_; ///< last sample from fed raw data to copy (if present)
  int expectedOrbitMessageTime_; ///< Expected orbit bunch time (needed to evaluate time differences)
  int mode_;
  std::set<HcalElectronicsId> unknownIds_,unknownIdsTrig_; ///< Recorded to limit number of times a log message is generated
};

#endif // HOTPUnpacker_h_included
