#include "EventFilter/L1TXRawToDigi/interface/HOTPUnpacker.h"
#include "EventFilter/HcalRawToDigi/interface/HcalDCCHeader.h"
#include "EventFilter/HcalRawToDigi/interface/HcalDTCHeader.h"
#include "EventFilter/HcalRawToDigi/interface/AMC13Header.h"
#include "EventFilter/HcalRawToDigi/interface/HcalHTRData.h"
#include "DataFormats/HcalDetId/interface/HcalOtherDetId.h"
#include "DataFormats/HcalDigi/interface/HcalQIESample.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "EventFilter/HcalRawToDigi/interface/HcalTTPUnpacker.h"
#include "EventFilter/HcalRawToDigi/plugins/HcalRawToDigi.h"

namespace HOTPUnpacker_impl {
  template <class DigiClass>
  const HcalQIESample* unpack(const HcalQIESample* startPoint, const HcalQIESample* limit, DigiClass& digi, int presamples, const HcalElectronicsId& eid, int startSample, int endSample, int expectedTime, const HcalHTRData& hhd) {
    // set parameters
    digi.setPresamples(presamples);
    digi.setReadoutIds(eid);

    int fiber=startPoint->fiber();
    int fiberchan=startPoint->fiberChan();
    uint32_t zsmask=hhd.zsBunchMask()>>startSample;
    digi.setZSInfo(hhd.isUnsuppressed(),hhd.wasMarkAndPassZS(fiber,fiberchan),zsmask);

    if (expectedTime>=0 && !hhd.isUnsuppressed()) {
      //      std::cout << hhd.getFibOrbMsgBCN(fiber) << " " << expectedTime << std::endl;
      digi.setFiberIdleOffset(hhd.getFibOrbMsgBCN(fiber)-expectedTime);
    }

    // what is my sample number?
    int myFiberChan=startPoint->fiberAndChan();
    int ncurr=0,ntaken=0;
    const HcalQIESample* qie_work=startPoint;
    while (qie_work!=limit && qie_work->fiberAndChan()==myFiberChan) {
      if (ncurr>=startSample && ncurr<=endSample) {
	digi.setSample(ntaken,*qie_work);
	++ntaken;
      }
      ncurr++;
      qie_work++;
    }
    digi.setSize(ntaken);
    return qie_work;
  }


  template <class DigiClass>
  const unsigned short* unpack_compact(const unsigned short* startPoint, const unsigned short* limit, DigiClass& digi, 
				       int presamples, const HcalElectronicsId& eid, int startSample, int endSample, 
				       int expectedTime, const HcalHTRData& hhd) {
    // set parameters
    digi.setPresamples(presamples);
    digi.setReadoutIds(eid);
    int flavor, error_flags, capid0, channelid;

    HcalHTRData::unpack_per_channel_header(*startPoint,flavor,error_flags,capid0,channelid);
    bool isCapRotating=!(error_flags&0x1);
    bool fiberErr=(error_flags&0x2);
    bool dataValid=!(error_flags&0x2);
    int fiberchan=channelid&0x3;
    int fiber=((channelid>>2)&0x7)+1;

    uint32_t zsmask=hhd.zsBunchMask()>>startSample;
    digi.setZSInfo(hhd.isUnsuppressed(),hhd.wasMarkAndPassZS(fiber,fiberchan),zsmask);

    if (expectedTime>=0 && !hhd.isUnsuppressed()) {
      //      std::cout << hhd.getFibOrbMsgBCN(fiber) << " " << expectedTime << std::endl;
      digi.setFiberIdleOffset(hhd.getFibOrbMsgBCN(fiber)-expectedTime);
    }

    // what is my sample number?
    int ncurr=0,ntaken=0;
    const unsigned short* qie_work=startPoint;
    // we branch here between normal (flavor=5) and error mode (flavor=6)
    if (flavor==5) {
      for (qie_work++; qie_work!=limit && !HcalHTRData::is_channel_header(*qie_work); qie_work++) {
	int capidn=(isCapRotating)?((capid0+ncurr)%4):(capid0);
	int capidn1=(isCapRotating)?((capid0+ncurr+1)%4):(capid0);
	// two samples in one...
	HcalQIESample s0((*qie_work)&0x7F,capidn,fiber,fiberchan,dataValid,fiberErr);
	HcalQIESample s1(((*qie_work)>>8)&0x7F,capidn1,fiber,fiberchan,dataValid,fiberErr);
	
	if (ncurr>=startSample && ncurr<=endSample) {
	  digi.setSample(ntaken,s0);
	  ++ntaken;
	}
	ncurr++;
	if (ncurr>=startSample && ncurr<=endSample) {
	  digi.setSample(ntaken,s1);
	  ++ntaken;
	}
	ncurr++;
      }
      digi.setSize(ntaken);
    } else if (flavor==6) {
      for (qie_work++; qie_work!=limit && !HcalHTRData::is_channel_header(*qie_work); qie_work++) {
	if (ncurr>=startSample && ncurr<=endSample) {
	  HcalQIESample sample((*qie_work)&0x7F,((*qie_work)>>8)&0x3,fiber,fiberchan,((*qie_work)>>10)&0x1,((*qie_work)>>11)&0x1);
	  digi.setSample(ntaken,sample);
	  ++ntaken;
	}
	ncurr++;
      }
      digi.setSize(ntaken);
    }
    return qie_work;
  }

}

static inline bool isTPGSOI(const HcalTriggerPrimitiveSample& s) {
  return (s.raw()&0x200)!=0;
}

static int slbChan(uint16_t theSample) { return (theSample>>11)&0x3; }

struct HOUnrolledTP { // parts of an HO trigger primitive, unpacked
  bool valid, checked;
  int ieta, iphi, samples, soi;
  unsigned int databits;
  HOUnrolledTP() {
    valid=false;
    checked=false;
    ieta=0;
    iphi=0;
    samples=0;
    soi=0;
    databits=0;
  }
  void setbit(int i) { databits|=(1<<i); }    
};


void HOTPUnpacker::unpack(const FEDRawData& raw, const HcalElectronicsMap& emap, HOCol& colls, HcalUnpackerReport& report, bool silent) {
  std::cout<<"inside the unpack functiin"<<std::endl;
  if (raw.size()<16) {
    if (!silent) edm::LogWarning("Invalid Data") << "Empty/invalid DCC data, size = " << raw.size();
    return;
  }

  // get the DCC header
  const HcalDCCHeader* dccHeader=(const HcalDCCHeader*)(raw.data());
  const HcalDTCHeader* dtcHeader=(const HcalDTCHeader*)(raw.data());
  bool is_VME_DCC=(dccHeader->getDCCDataFormatVersion()<0x10) || ((mode_&0x1)==0);
  
  int dccid=(is_VME_DCC)?(dccHeader->getSourceId()-sourceIdOffset_):(dtcHeader->getSourceId()-sourceIdOffset_);
  
  std::cout<<"got the DCC header"<<std::endl;
  // walk through the HTR data.  For the uTCA, use spigot=slot+1
  HcalHTRData htr;
  const unsigned short* daq_first, *daq_last, *tp_first, *tp_last;
  const HcalTriggerPrimitiveSample *tp_begin, *tp_end, *tp_work; 
  //  for (int spigot=0; spigot<HcalDCCHeader::SPIGOT_COUNT; spigot++) {
  for (int spigot=0; spigot<100; spigot++) {
    
    std::cout<<"is_VME_DCC"<<is_VME_DCC<<" spigot : "<<spigot<<std::endl;
    if (is_VME_DCC) {
      if (!dccHeader->getSpigotPresent(spigot)) continue;
      
      std::cout<<"spigot presnts"<<std::endl;
      int retval=dccHeader->getSpigotData(spigot,htr,raw.size());
      std::cout<<"retval is:"<<retval<<std::endl;
      if (retval!=0) {
	if (retval==-1) {
	  if (!silent) edm::LogWarning("Invalid Data") << "Invalid HTR data (data beyond payload size) observed on spigot " << spigot << " of DCC with source id " << dccHeader->getSourceId();
	  report.countSpigotFormatError();
	}
	continue;
      }
      
      std::cout<<"zeroth set of checks are done"<<std::endl;
      
      // check
      if (dccHeader->getSpigotCRCError(spigot)) {
	if (!silent) 
	  edm::LogWarning("Invalid Data") << "CRC Error on HTR data observed on spigot " << spigot << " of DCC with source id " << dccHeader->getSourceId();
	report.countSpigotFormatError();
	continue;
      } 
    }
    std::cout<<"first set of checks are done"<<std::endl;

    // check for EE
    if (htr.isEmptyEvent()) {
      report.countEmptyEventSpigot();
    }
    if (htr.isOverflowWarning()) {
      report.countOFWSpigot();
    }
    if (htr.isBusy()) {
      report.countBusySpigot();
    }
    if (!htr.check()) {
      if (!silent) 
	edm::LogWarning("Invalid Data") << "Invalid HTR data observed on spigot " << spigot << " of DCC with source id " << dccHeader->getSourceId();
      report.countSpigotFormatError();
      continue;
    }  
   
    std::cout<<"second  set of checks are done"<<std::endl;
    if (htr.getFirmwareFlavor()>=0x80) {
      if (!silent) edm::LogWarning("HcalUnpackerReport") << "Skipping data on spigot " << spigot << " of DCC with source id " << dccHeader->getSourceId() << " which is of unknown flavor " << htr.getFirmwareFlavor();
      continue;
    }

    // get pointers
    std::cout<<"will getting the pointers"<<std::endl;
    htr.dataPointers(&daq_first,&daq_last,&tp_first,&tp_last);
    unsigned int smid=htr.getSubmodule();
    int htr_tb=smid&0x1;
    int htr_slot=(smid>>1)&0x1F;
    int htr_cr=(smid>>6)&0x1F;
    
    tp_begin=(HcalTriggerPrimitiveSample*)tp_first;
    tp_end=(HcalTriggerPrimitiveSample*)(tp_last+1); // one beyond last..
    
    /// work through the samples
    bool isHOtpg=htr.getFormatVersion()>=3 && htr.getFirmwareFlavor()==0; // HO is flavor zero
  
    /*
      Unpack the HO trigger primitives
    */
    std::cout<<"will enter in the isHOtps loop"<<std::endl;
    if (isHOtpg) {
      HOUnrolledTP unrolled[24];
      for (tp_work=tp_begin; tp_work!=tp_end; tp_work++) {
	if (tp_work->raw()==0xFFFF) continue; // filler word
	int sector=slbChan(tp_work->raw());
	if (sector>2) continue;

	for (int ibit=0; ibit<8; ibit++) {
	  int linear=sector*8+ibit; 
	  if (!unrolled[linear].checked) {
	    unrolled[linear].checked=true;
	    int fiber=(linear/3)+1;
	    int fc=(linear%3);
	    // electronics id (use precision match for HO TP)
	    HcalElectronicsId eid(fc,fiber,spigot,dccid);	
	    eid.setHTR(htr_cr,htr_slot,htr_tb);
	    DetId did=emap.lookup(eid);
	    if (!did.null()) {
	      if (did.det()==DetId::Hcal && ((HcalSubdetector)did.subdetId())==HcalOuter ) {
		HcalDetId hid(did);
		unrolled[linear].valid=true;
		unrolled[linear].ieta=hid.ieta();
		unrolled[linear].iphi=hid.iphi();
	      }
	    } else {
	      report.countUnmappedTPDigi(eid);
	    }
	  }
	  if (unrolled[linear].valid) {
	    if (isTPGSOI(*tp_work)) unrolled[linear].soi=unrolled[linear].samples;
	    if (tp_work->raw()&(1<<ibit)) unrolled[linear].setbit(unrolled[linear].samples);
	    unrolled[linear].samples++;
	  }
	}
      }
      for (int i=0; i<24; i++) {
	std::cout<<"will do the filling, "<<i<<std::endl;
	if (unrolled[i].valid) { 
	  std::cout<<"all info:"<<unrolled[i].ieta<<", "<<unrolled[i].iphi<<","<<unrolled[i].databits<<std::endl;
	  colls.etaCont->push_back(unrolled[i].ieta);
	  colls.phiCont->push_back(unrolled[i].iphi);
	  colls.sampleCont->push_back(unrolled[i].samples);
	  colls.soiCont->push_back(unrolled[i].soi);
	  colls.databitsCont->push_back(unrolled[i].databits);
	}
      }
    }
  }
}
      

HOTPUnpacker::HOCol::HOCol() {
  etaCont=0;
  phiCont=0;
  sampleCont=0;
  soiCont=0;
  databitsCont=0;
}

      

