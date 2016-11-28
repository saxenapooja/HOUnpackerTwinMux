#ifndef HODIGI_TWINMUX_H
#define HODIGI_TWINMUX_H

#include <boost/cstdint.hpp>
#include <ostream>
//#include "DataFormats/HcalDetId/interface/HcalDetId.h"

/** \class HOTwinMuxDigi
  *  Simple container packer/unpacker for HO in TwinMUX
  *  Trigger Primitive from an HO HTR, via oSLB links, to TwinMUX
  *
  *  \author Saxena, Pooja - DESY
  */

class HOTwinMuxDigi {
 public:
  HOTwinMuxDigi() {theTP_HO=0;}
  HOTwinMuxDigi(uint64_t data) {theTP_HO = data;}
  HOTwinMuxDigi(int ieta, int iphi, int mip, int sector, int wheel, int index, int link, bool validbit );

  /// get raw packed HO 
  uint64_t raw() const {return theTP_HO; }

  /// get the raw ieta value
  int raw_ieta() const {return theTP_HO&0x1F; }

  /// get the sign of ieta (int: +/- 1)
  int ieta_sign() const {return ( (theTP_HO&0x10)?(-1):(+1)); }

  /// get the absolute value of ieta
  int ieta_abs() const {return (theTP_HO&0x000F); }

  /// get the signed ieta value
  int ieta() const {return (ieta_abs() * ieta_sign()); }

  /// get the raw iphi value
  int iphi() const {return (theTP_HO>>5)&0x007F; }

  /// get the mip value
  int mip() const {return (theTP_HO>>16)&0x1; }

  /// get the sector value
  int sector() const {return (theTP_HO>>17)&0xF; }
  
  /// get the raw wheel value
  int raw_wheel() const {return (theTP_HO>>21)&0x7; }

  /// get the sign of wheel (int: +/- 1)                                                                                                                                   
  int wheel_sign() const {return ( ( (theTP_HO>>23)&0x10) ?(-1):(+1)); }

  /// get the absolute value of wheel                                                                                                                                       
  int wheel_abs() const {return (theTP_HO>>21)&0x03; }
  
  /// get the signed wheel value                                                                                                                                          
  int wheel() const {return (wheel_abs() * wheel_sign()); }

  /// get the index
  int index() const {return (theTP_HO>>24)&0x3F; }

  /// get the link value
  int link() const {return (theTP_HO>>30)&0x3; }

  /// get the valid bit
  int validbit() const {return (theTP_HO>>32)&0x1; }

  static const int HO_SECTOR_MAX = 12;

 private:
  uint64_t theTP_HO;
};

std::ostream& operator<<(std::ostream&, const HOTwinMuxDigi&);

#endif


  
