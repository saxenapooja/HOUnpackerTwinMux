#include "DataFormats/L1Trigger/interface/HOTwinMuxDigi.h"
#include <cstdio>
#include <iostream>

using namespace std;

HOTwinMuxDigi::HOTwinMuxDigi(int ieta, int iphi, int bx, int mip, int validbit, int wheel, int sector, int index, int link) {

  //HOTwinMuxDigi(int ieta, int iphi, int mip, int sector, int wheel, int index, int link, bool validbit, signed int bx) {
  //  printf("filling the bit values here!!!");
  if((ieta > 15) || (ieta < -15)) printf("HOTwinMuxDigi:  ieta out of range");
  if ((mip < 0) || (mip > 1)) printf("HOTwinMuxDigi: mip value out of range");
 
  //  std::cout<<"ieta, iphi :"<< ieta<<", "<< iphi << std::endl;

  theTP_HO=(abs(ieta)&0xF) | ((ieta<0)?(0x10) : (0x00)) |
    ((iphi&0x7F)<<5) |  
    ((abs(bx)&0x1)<<12) | ((bx<0)?(0x2000) : (0x00)) |
    ((mip&0x1)<<14) | 
    ((validbit&0x1)<<15) |
    ((abs(wheel)&0x3)<<16) | ((wheel<0)?(0x40000) : (0x00)) |
    ((sector&0xF)<<19) | 
    ((index &0x1F)<<23) | ((link&0x3)<<28);
}

std::ostream& operator<<(std::ostream& s, const HOTwinMuxDigi& HOtp) {
  s << "HO digi in TwinMUX "<<HOtp.ieta() <<", " << HOtp.iphi() << ", " << HOtp.mip();
  s << "(wh, sec, index, link)" <<HOtp.wheel()<<", " << HOtp.sector() <<", "<< HOtp.index() << ", " << HOtp.link(); 
  s << "validbit is : "<< HOtp.validbit();
  return s;
}    
