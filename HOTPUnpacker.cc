#include "EventFilter/L1TXRawToDigi/plugins/HOTPUnpacker.h"
#include "EventFilter/HcalRawToDigi/interface/HcalDCCHeader.h"
#include "EventFilter/HcalRawToDigi/interface/HcalDTCHeader.h"
#include "EventFilter/HcalRawToDigi/interface/HcalHTRData.h"
#include "DataFormats/HcalDetId/interface/HcalOtherDetId.h"
#include "DataFormats/HcalDigi/interface/HcalQIESample.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "EventFilter/HcalRawToDigi/interface/HcalTTPUnpacker.h"

//namespace HOTPUnpacker_impl {


static inline bool isTPGSOI(const HcalTriggerPrimitiveSample& s) {
  return (s.raw()&0x200)!=0;
}


struct HOTPUnpacker::HOUnrolledTP { // parts of an HO trigger primitive, unpacked
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

static int slbChan(uint16_t theSample) { return (theSample>>11)&0x3; }

void HOTPUnpacker::unpack(const FEDRawData& raw, 
			  const HcalElectronicsMap& emap,
			  // std::vector<HOTriggerPrimitiveDigi>& tphoCont,
			  Collections& colls,
			  HcalUnpackerReport& report, 
			  bool silent) {

  if (raw.size()<16) {
    if (!silent) edm::LogWarning("Invalid Data") << "Empty/invalid DCC data, size = " << raw.size();
    return;
  }

  // get the DCC header
  const HcalDCCHeader* dccHeader=(const HcalDCCHeader*)(raw.data());
  const HcalDTCHeader* dtcHeader=(const HcalDTCHeader*)(raw.data());
  bool is_VME_DCC=(dccHeader->getDCCDataFormatVersion()<0x10) || ((mode_&0x1)==0);
  
  int dccid=(is_VME_DCC)?(dccHeader->getSourceId()-sourceIdOffset_):(dtcHeader->getSourceId()-sourceIdOffset_);

  // check the summary status
  
  // walk through the HTR data.  For the uTCA, use spigot=slot+1
  HcalHTRData htr;
  const unsigned short* daq_first, *daq_last, *tp_first, *tp_last;
  //  const HcalQIESample* qie_begin, *qie_end, *qie_work;
  const HcalTriggerPrimitiveSample *tp_begin, *tp_end, *tp_work; 
  for (int spigot=0; spigot<HcalDCCHeader::SPIGOT_COUNT; spigot++) {
    
    if (is_VME_DCC) {
      if (!dccHeader->getSpigotPresent(spigot)) continue;

      int retval=dccHeader->getSpigotData(spigot,htr,raw.size());
      if (retval!=0) {
	if (retval==-1) {
	  if (!silent) edm::LogWarning("Invalid Data") << "Invalid HTR data (data beyond payload size) observed on spigot " << spigot << " of DCC with source id " << dccHeader->getSourceId();
	  report.countSpigotFormatError();
	}
	continue;
      }
      // check
      if (dccHeader->getSpigotCRCError(spigot)) {
	if (!silent) 
	  edm::LogWarning("Invalid Data") << "CRC Error on HTR data observed on spigot " << spigot << " of DCC with source id " << dccHeader->getSourceId();
	report.countSpigotFormatError();
	continue;
      } 
    } else { // is_uTCA (!is_VME_DCC)
      int slot=spigot+1;
      if (slot>HcalDTCHeader::MAXIMUM_SLOT) continue;

      if (!dtcHeader->getSlotPresent(slot)) continue;

      int retval=dtcHeader->getSlotData(slot,htr,raw.size());
      if (retval!=0) {
	if (retval==-1) {
	  if (!silent) edm::LogWarning("Invalid Data") << "Invalid uHTR data (data beyond payload size) observed on slot " << slot << " of DTC with source id " << dtcHeader->getSourceId();
	  report.countSpigotFormatError();
	}
	continue;
      }
      // check
      if (dtcHeader->getSlotCRCError(slot)) {
	if (!silent) 
	  edm::LogWarning("Invalid Data") << "CRC Error on uHTR data observed on slot " << slot << " of DTC with source id " << dtcHeader->getSourceId();
	report.countSpigotFormatError();
	continue;
      } 
    }


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
   
    if (htr.getFirmwareFlavor()>=0x80) {
      if (!silent) edm::LogWarning("HOTPUnpacker") << "Skipping data on spigot " << spigot << " of DCC with source id " << dccHeader->getSourceId() << " which is of unknown flavor " << htr.getFirmwareFlavor();
      continue;
    }

    // calculate "real" number of presamples
    //    int nps=htr.getNPS()-startSample_;
    
    // get pointers
    htr.dataPointers(&daq_first,&daq_last,&tp_first,&tp_last);
    unsigned int smid=htr.getSubmodule();
    int htr_tb=smid&0x1;
    int htr_slot=(smid>>1)&0x1F;
    int htr_cr=(smid>>6)&0x1F;
    
    tp_begin=(HcalTriggerPrimitiveSample*)tp_first;
    tp_end=(HcalTriggerPrimitiveSample*)(tp_last+1); // one beyond last..
    
    /// work through the samples
    //    int currFiberChan=0x3F; // invalid fiber+channel...
    // int ncurr=0;
    // bool valid=false;

    //    bool tpgSOIbitInUse=htr.getFormatVersion()>=3; // version 3 and later
    bool isHOtpg=htr.getFormatVersion()>=3 && htr.getFirmwareFlavor()==0; // HO is flavor zero
    //    int npre=0;

    /*
      Unpack the trigger primitives
    */
    if (isHOtpg) {
      HOUnrolledTP unrolled[24];
      for (tp_work=tp_begin; tp_work!=tp_end; tp_work++) {
	if (tp_work->raw()==0xFFFF) continue; // filler word
	//	int sector=tp_work->slbChan();
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
	if (unrolled[i].valid) 
	  colls.tphoCont->push_back(HOTriggerPrimitiveDigi(
						    unrolled[i].ieta,
						    unrolled[i].iphi,
						    unrolled[i].samples,
						    unrolled[i].soi,
						    unrolled[i].databits));
      }
    }
    
  }
}

HOTPUnpacker::Collections::Collections() {
  tphoCont=0;
}


// void HOTPUnpacker::unpack(const FEDRawData& raw, 
// 			  const HcalElectronicsMap& emap, 
// 			  std::vector<HODataFrame>& container, 
// 			  std::vector<HcalTriggerPrimitiveDigi>& tp) {
//   Collections c;
//   c.hoCont=&container;
//   c.tpCont=&tp;
//   HOTPUnpackerReport r;
//   unpack(raw,emap,c,r);
// }

