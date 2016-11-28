#include "DataFormats/L1Trigger/interface/HOTwinMuxDigi.h"
#include <cstdio>

HOTwinMuxDigi::HOTwinMuxDigi(int ieta, int iphi, int mip, int sector, int wheel, int index, int link, bool validbit) {
  if((ieta > 15) || (ieta < -15)) printf("HOTwinMuxDigi:  ieta out of range");
  if ((mip < 0) || (mip > 1)) printf("HOTwinMuxDigi: mip value out of range");
  
  theTP_HO=(abs(ieta)&0xF) | ((ieta<0)?(0x10) : (0x00)) |
    ((iphi&0x7F)<<5) |  ((mip&0x1)<<16) |
    ((sector&0xF)<<17) | 
    ((abs(wheel)&0x3)<<21) | ((wheel<0)?(0x800000) : (0x00)) |
    ((index &0x3F)<<24) | ((link&0x3)<<30) |
    ((validbit&0x1)<<32);
}

std::ostream& operator<<(std::ostream& s, const HOTwinMuxDigi& HOtp) {
  s << "HO digi in TwinMUX "<<HOtp.ieta() <<", " << HOtp.iphi() << ", " << HOtp.mip();
  s << "(wh, sec, index, link)" <<HOtp.wheel()<<", " << HOtp.sector() <<", "<< HOtp.index() << ", " << HOtp.link(); 
  s << "validbit is : "<< HOtp.validbit();
  return s;
}    
