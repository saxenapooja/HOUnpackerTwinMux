//-------------------------------------------------
//
//   Class: DTTM7FEDReader
//
//   L1 DT TwinMux Raw-to-Digi
//
//
//
//   Author :
//   C. F. Bedoya  - CIEMAT
//   G. Codispoti -- INFN Bologna
//   J. Pazzini   -- INFN Padova
//
//--------------------------------------------------

#include "L1TTwinMuxRawToDigi.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "EventFilter/HcalRawToDigi/interface/HcalDCCHeader.h"
#include "EventFilter/HcalRawToDigi/interface/HcalHTRData.h"
#include "EventFilter/HcalRawToDigi/interface/HcalDTCHeader.h" 
#include "FWCore/Framework/interface/ESHandle.h" 
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include <iostream>
#include <fstream>

L1TTwinMuxRawToDigi::L1TTwinMuxRawToDigi(const edm::ParameterSet& pset) :

  debug_( pset.getUntrackedParameter<bool>("debug", false) ), 
  nfeds_(0),
  DTTM7InputTag_( pset.getParameter<edm::InputTag>("DTTM7_FED_Source") ),
  feds_( pset.getUntrackedParameter<std::vector<int> >("feds", std::vector<int>()) ),
  wheels_( pset.getUntrackedParameter<std::vector<int> >("wheels", std::vector<int>())),
  amcsecmap_( pset.getUntrackedParameter<std::vector<long long int> >("amcsecmap", std::vector<long long int>())),
  //  hounpacker_(pset.getUntrackedParameter<int>("HcalFirstFED",int(FEDNumbering::MINHCALFEDID)),pset.getParameter<int>("firstSample"),pset.getParameter<int>("lastSample")),
  hofedUnpackList_(pset.getUntrackedParameter<std::vector<int> >("FEDs", std::vector<int>())),
  hofirstFED_(pset.getUntrackedParameter<int>("HcalFirstFED",FEDNumbering::MINHCALFEDID)),
  silent_(pset.getUntrackedParameter<bool>("silent",true)),
  complainEmptyData_(pset.getUntrackedParameter<bool>("ComplainEmptyData",false))
{
    

  if (hofedUnpackList_.empty()) {
    for (int i=FEDNumbering::MINHCALFEDID; i<=FEDNumbering::MAXHCALFEDID; i++)
      // (int i=724; i<=FEDNumbering::MAXHCALFEDID; i++) // 724-731 are HO FEDs
      hofedUnpackList_.push_back(i);
  }


  produces<L1MuDTChambPhContainer>();
  produces<L1MuDTChambThContainer>();
  produces<HOTrigPrimDigiCollection>();
  produces<HcalUnpackerReport>();

  Raw_token = consumes<FEDRawDataCollection> (DTTM7InputTag_);
 
  nfeds_ = feds_.size();
  
  if ( nfeds_ != wheels_.size() )
    throw cms::Exception("TwinMux_unpacker") << "Configuration file error. Size of \'wheels\' and \'feds\' differs.\n";  

  if ( amcsecmap_.size() != wheels_.size() )
    throw cms::Exception("TwinMux_unpacker") << "Configuration file error. Size of \'wheels\' and \'amcsecmap\' differs.\n";  
    
  for (size_t wh_i = 0; wh_i < amcsecmap_.size(); ++wh_i){
    std::array<short, 12> whmap;      
    for (size_t amc_i = 1; amc_i < 13; ++amc_i ){
      short shift = (12-amc_i)*4;
      whmap[amc_i-1] = ( amcsecmap_[wh_i] >> shift ) & 0xF;
    }
    amcsec_.push_back(whmap);
  }
   
}

L1TTwinMuxRawToDigi::~L1TTwinMuxRawToDigi(){}

void L1TTwinMuxRawToDigi::produce(edm::Event& e, 
                             const edm::EventSetup& c) {

  std::auto_ptr<L1MuDTChambPhContainer> TM7phi_product(new L1MuDTChambPhContainer);
  std::auto_ptr<L1MuDTChambThContainer> TM7the_product(new L1MuDTChambThContainer);

  L1MuDTChambPhContainer::Phi_Container phi_data;
  L1MuDTChambThContainer::The_Container the_data;
  //  HOTriggerPrimitiveDigi ho_data;

  if(dodebug)  std::cout<<"************ Start the produce loop for event: "<< (e.eventAuxiliary().event()) << std::endl;

  // HO
  edm::ESHandle<HcalDbService> pSetup; 
  c.get<HcalDbRecord>().get( pSetup );
  const HcalElectronicsMap*  readoutMap = pSetup->getHcalMapping();

  if(dodebug)  std::cout<<"unpacked the hofedlist, size is :"<< hofedUnpackList_.size() << std::endl;

  std::vector<HOTriggerPrimitiveDigi> hotp;
  Collections colls;
  colls.tphoCont = &hotp;

  if(dodebug) std::cout<<"Entering in the fillRawData loop"<< std::endl;
  if ( !fillRawData(e, phi_data, the_data, colls, readoutMap) ) return;

  TM7phi_product->setContainer(phi_data);
  TM7the_product->setContainer(the_data);

  e.put(TM7phi_product);
  e.put(TM7the_product);

  //HO
  std::auto_ptr<HcalUnpackerReport> report(new HcalUnpackerReport);
  std::auto_ptr<HOTrigPrimDigiCollection> hotp_product(new HOTrigPrimDigiCollection());
  hotp_product->swap_contents(hotp);
  hotp_product->sort();
  e.put(hotp_product);
}


bool L1TTwinMuxRawToDigi::fillRawData( edm::Event& e,
				       L1MuDTChambPhContainer::Phi_Container& phi_data,
				       L1MuDTChambThContainer::The_Container& the_data,
				       Collections& colls,
				       const HcalElectronicsMap *readoutMap) {
  edm::Handle<FEDRawDataCollection> data;
  e.getByToken( Raw_token, data );

  if(dodebug)   std::cout<<"In fillRawData loop, will start analyzing the TMux feds, No.:"<<nfeds_<< std::endl;
  for ( size_t w_i = 0; w_i < nfeds_; ++w_i ) {
    if(dodebug) std::cout<<"processing the FED: "<<feds_[w_i]<< std::endl;
    processFed( feds_[w_i], wheels_[w_i], amcsec_[w_i], data, phi_data, the_data, colls, readoutMap );
  }
  
  return true;
}

int L1TTwinMuxRawToDigi::normBx( int bx_, 
                            int bxCnt_ ){
    
    int bxNorm_ = bx_ - bxCnt_;    
    if ( abs( bxNorm_ ) < 3000 ) return bxNorm_; 
    
    if ( bxNorm_ > 0 ) return bxNorm_ - 3564;
    if ( bxNorm_ < 0 ) return bxNorm_ + 3564;
    
    return -99;
    
}

int L1TTwinMuxRawToDigi::radAngConversion( int radAng_  ) {
    
    if (radAng_>2047) 
        return radAng_-4096;

    return radAng_;
    
}

int L1TTwinMuxRawToDigi::benAngConversion( int benAng_  ) {
    
    if (benAng_>511) 
        return benAng_-1024;

    return benAng_;
    
}


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

static int slbChan(uint16_t theSample) { return (theSample>>11)&0x3; }

                                                                                                                             
static inline bool isTPGSOI(const HcalTriggerPrimitiveSample& s) {
  return (s.raw()&0x200)!=0;
}

void L1TTwinMuxRawToDigi::processFed( int twinMuxFed, 
				      int twinMuxWheel,
				      std::array<short, 12> twinMuxAmcSec,
				      edm::Handle<FEDRawDataCollection> data,
				      L1MuDTChambPhContainer::Phi_Container& phiSegments,
				      L1MuDTChambThContainer::The_Container& theSegments,
				      Collections& colls,
				      const HcalElectronicsMap* readoutMap) {
  
  if(dodebug)  std::cout<<"In processFed loop"<< std::endl;
  /// Container
  std::vector<long> DTTM7WordContainer;

  /// Debug
  std::ofstream logfile;
  if ( debug_ ) {
    std::ostringstream fname;
    fname << "eventDump_" <<  twinMuxFed << ".txt";
    logfile.open( fname.str() );
  }

  /// Header
  if(dodebug)  std::cout<<"checking if there is TMdata?"<<std::endl;
  FEDRawData TM7data = data->FEDData(twinMuxFed); 
  if ( TM7data.size() == 0 ) return;

  if(dodebug) std::cout<<"TM7data has data!"<<std::endl;

  /// Variables
  LineFED_ = TM7data.data();
  int nline  = 0; // counting already include header
  long dataWord = 0;
  int newCRC = 0xFFFF;

  ///--> Header - line 1 [must start with 0x5]
  readline( nline, dataWord );
  calcCRC( dataWord, newCRC );
  
  std::cout<<"dataWord : "<< dataWord<<std::endl;
  std::cout<<"dataWord >> 8: "<< (dataWord >> 8) <<std::endl;
  std::cout<<"dataWord >> 8) & 0xfff: "<< ((dataWord >> 8) & 0xFFF) <<std::endl;


  int TM7fedId = ( dataWord >> 8 ) & 0xFFF;  // positions 8 -> 19
  /*** NOT UNPACKED  
  int bunchCnt = ( dataWord >> 20 ) & 0xFFF;  // positions 20 -> 31
  int eventCnt = ( dataWord >> 32 ) & 0xFFFFFF;  // positions 32 -> 55
  ***/
  int BOEevTy  = ( dataWord >> 60 ) & 0xF;  // positions 60 -> 63

  int linecounter = 0;
  if ( debug_ ) logfile << '[' << ++linecounter << "]\t"
                        << std::hex << dataWord << std::dec << "\t|\t"
                        << "BOEevTy " << BOEevTy << '\t'
                        << "TM7fedId "  << TM7fedId << '\n';

  if ( (BOEevTy != 0x5) || ( TM7fedId != twinMuxFed ) ) {
            
    edm::LogWarning("TwinMux_unpacker") << "Not a TM7 of FED " 
                                        << twinMuxFed << " header "
                                        << std::hex << dataWord;
    return;
    
  }

  ///--> Header - line 2
  readline( nline, dataWord );
  calcCRC( dataWord, newCRC );

  std::map<int, int> AMCsizes;
  /*** NOT UNPACKED  
  int orbit = ( dataWord >> 4 ) & 0xFFFFFFFF;  // positions 4 -> 35
  ***/
  int nAMC = ( dataWord >> 52 ) & 0xF;  // positions 52 -> 55

  if ( debug_ ) logfile << '[' << ++linecounter << "]\t" << std::hex
                        << dataWord << std::dec << "\t|\t"
                        << "nAMC " << nAMC << '\n';

  ///--> AMC - line 3 to 3+nAMC
  for ( int j = 0; j < nAMC; ++j ) {
  
    readline( nline, dataWord ); 
    calcCRC( dataWord, newCRC );
   
    int AMCno = (dataWord >> 16 ) & 0xF;  // positions 16 -> 19
    /*** NOT UNPACKED  
    int TM7boardID = dataWord & 0xFFFF;  // positions 0 -> 15
    int bulkno = (dataWord >> 20 ) & 0xFF;  // positions 20 -> 27
    ***/
    if ( (AMCno < 1) || (AMCno > 12) ) {
        edm::LogWarning("TwinMux_unpacker") << "AMCnumber " << std::dec << AMCno 
                                            << " out of range (1-12)";
        return;
    }    

    AMCsizes[AMCno] = ( dataWord >> 32 ) & 0xFFFFFF;  // positions 32 -> 55

    if ( debug_ ) logfile << '[' << ++linecounter << "]\t"
                          << std::hex << dataWord
                          << std::dec << "\t|\t"
                          << "AMCsizes[" << AMCno << "] "
                          << AMCsizes[AMCno]
                          << std::dec << '\n';
  }

  ///--> Store payloads
  std::map<int,int>::iterator AMCiterator = AMCsizes.begin();
  std::map<int,int>::iterator AMCitend = AMCsizes.end();  
  for ( ; AMCiterator != AMCitend; ++AMCiterator ) {
      
    for ( int k=0; k<AMCiterator->second; ++k) {
        
       readline( nline, dataWord );
       calcCRC( dataWord, newCRC);
       DTTM7WordContainer.push_back( dataWord );
    }
  }  

  ///--> Trailer - line 1
  readline( nline, dataWord );
  calcCRC( dataWord, newCRC);

  ///--> Trailer - line 2 [must start with 0xA]
  readline( nline, dataWord );
  calcCRC( dataWord & 0xFFFFFFFF0000FFFF, newCRC); /// needed not to put crc in crc calc

  ///--> AMC trailer - line 2
  int chkEOE = (dataWord >> 60 ) & 0xF;  // positions 60 -> 63
  int CRC = ( dataWord >> 16 ) & 0xFFFF; // positions 17 ->32
  int evtLgth = ( dataWord >> 32 ) & 0xFFFFFF; // positions 33 ->56

  if ( chkEOE != 0xA ) {
      edm::LogWarning("TwinMux_unpacker") << "AMC block closing line " << std::hex << dataWord 
                                          << std::dec << " does not start with 0xA";
      return;
  }    

  if ( debug_ ) logfile << "\tevtLgth " << std::hex
                        << evtLgth << "\tCRC " << CRC << std::dec << '\n';

  if ( nline != evtLgth ) {
    edm::LogWarning("TwinMux_unpacker") << "Number of words read " << std::dec << nline 
                                        << " and event length " << std::dec << evtLgth 
                                        << " differ ";
    return;
  }

  if ( newCRC != CRC ) {
    edm::LogWarning("TwinMux_unpacker") << "Calculated CRC " << std::hex << newCRC 
                                        << " differs from CRC in trailer " << std::hex << CRC;
    return;
  }

  // --> Analyze event 
  std::vector<long>::iterator DTTM7iterator = DTTM7WordContainer.begin();
  std::vector<long>::iterator DTTM7itend = DTTM7WordContainer.end();

  if(dodebug) std::cout<<"will analyze the DTTM7 word!"<<std::endl;

  int lcounter = 0;
  for ( ; DTTM7iterator != DTTM7itend; ++DTTM7iterator ) {
  
    dataWord  = (*DTTM7iterator);
    int dataLenght = (dataWord & 0xFFFFF);         // positions 0 -> 19
    int bxCounter  = (dataWord >> 20 ) & 0xFFF;    // positions 20 -> 31
    int event      = (dataWord >> 32 ) & 0xFFFFFF; // positions 32 -> 55
    int AMC_ID     = (dataWord >> 56 ) & 0xF;      // positions 56 -> 59
    int control    = (dataWord >> 60 ) & 0xF;      // positions 59 -> 63 
    int wheel      = twinMuxWheel;
    
    if( ( AMC_ID < 1 ) or ( AMC_ID > 12 ) ) {
      edm::LogWarning("TwinMux_unpacker") << "%%%%%% AMC_ID OUT OF RANGE \n"
                                          << " TM7fedId "     << TM7fedId
                                          << " AMC_ID "       << AMC_ID;
      break;
    }
    
    int sector     = twinMuxAmcSec[AMC_ID-1];
    if(dodebug) std::cout<<"sector analysed is "<<sector<<std::endl;   
 
    if( ( sector < 1 ) or ( sector > 12 ) ) {
      if( sector != 15 ) edm::LogWarning("TwinMux_unpacker") << "%%%%%% VALID AMC_ID POINTS TO SECTOR OUT OF RANGE \n"
                                                             << " TM7fedId "     << TM7fedId
                                                             << " AMC_ID "       << AMC_ID
                                                             << " wheel "        << wheel
                                                             << " sector "       << sector;         
      break;
    }

    if ( debug_ ) logfile << '[' << ++lcounter << "]\t"
                          << std::hex << dataWord << std::dec << "\t|\t"
                          << "AMC_ID "     << AMC_ID << '\t' 
                          << "control "    << control   << '\t' 
                          << "event "      << event  << '\t' 
                          << "bxCounter "  << bxCounter  << '\t' 
                          << "dataLenght " << dataLenght << '\n';

    ++DTTM7iterator; // User word empty  /// ==>> increment 2

    
if( DTTM7iterator == DTTM7itend ) {
      edm::LogInfo("TwinMux_unpacker") << "TRAILING WORD AS A PAYLOAD END in FED " 
                                       << std::hex << TM7fedId 
                                       << std::hex << dataWord 
                                       << std::dec<< " [it pos " 
                                       << int(DTTM7iterator - DTTM7itend)  << " ]";
      break;
    }
  
    if(dodebug) std::cout<<"clear the test:01"<<std::endl;    

    dataWord = (*DTTM7iterator);
    int boardID   = (dataWord & 0xFFFF); // positions  0 -> 15
    int orbit     = (dataWord >> 16 ) & 0xFFFF; // positions 15 -> 32
    
    if ( DTTM7iterator == DTTM7itend ) {
      edm::LogWarning("TwinMux_unpacker") << "%%%%%% AMC_ID " << AMC_ID
                                          << " control "      << control
                                          << " event "        << event
                                          << " bxCounter "    << bxCounter
                                          << " size "         << dataLenght
                                          << " orbit "        << orbit
                                          << " board "        << boardID
                                          << " AMCsizes "     << AMCsizes[AMC_ID]
                                          << " it pos "       << int(DTTM7iterator - DTTM7itend);
      break;
    }

    if(dodebug) std::cout<<"clear the test:02"<<std::endl;    

    if (debug_ ) logfile << '[' << ++lcounter << "]\t" 
                         << std::hex << dataWord 
                         << std::dec << "\t|\t" 
                         << " orbit " << orbit
                         << " board " << boardID << '\n';

    int AMCsize = AMCsizes[AMC_ID] - 1; /// do not consider the trailer
    int bxID =  99;
    int bc0  = -99;
    int bxNr = -99;
    
    /// 2 words already read, last removed because trailer with CRC
    for ( int tm7eventsize = 2; tm7eventsize < AMCsize; ++tm7eventsize ) {
  
      ++DTTM7iterator; /// ==>> increment 3   
      if ( DTTM7iterator == DTTM7itend ) {
          
        edm::LogWarning("TwinMux_unpacker") << "UNEXPECTED END OF PAYLOAD INSIDE CHAMBER DESCRIPTION"
                                          << " [it pos " << int(DTTM7iterator - DTTM7itend)  << " ]" ;
        break;
        
      }

      if(dodebug) std::cout<<"clear the test:03"<<std::endl;    
      
      long dataWordSub = (*DTTM7iterator);
      int selector = ( dataWordSub >> 60 ) & 0xF; // positions 60 -> 63

      if(dodebug) std::cout<<"selector is :"<< selector << std::endl;

      if ( selector == 0x4 ) { //TSC word
	
        bxID = ( dataWordSub >> 48 ) & 0xFFF; // positions 48 -> 60
        bc0  = ( dataWordSub >> 22 ) & 0x1; // positions 22 -> 23
        bxNr = normBx(bxID, bxCounter); /// bx normalized to the bxcounter
	
        if ( debug_ ) logfile << '[' << ++lcounter << "]\t" << std::hex 
                              << dataWordSub << std::dec
                              << "\t TSC WORD\t"
                              << "bxID " << bxID << '\t'
                              << "bc0  " << bc0  << '\n';

      }//TSC WORD
    
      else if ( selector == 0x1 ) { //MB1/2 word
	std::cout<<"dataWordSub : "<< dataWordSub << std::endl;

        int mb2_phi =    ( dataWordSub & 0xFFF);        // positions  0 -> 11
        int mb2_phib =   ( dataWordSub >> 12 ) & 0x3FF; // positions 12 -> 21
        int mb2_qual =   ( dataWordSub >> 22 ) & 0x7;   // positions 22 -> 24
        int mb2_ts2tag = ( dataWordSub >> 28 ) & 0x1;   // positions 28
        /*** NOT UNPACKED  
        int mb2_parity = ( dataWordSub >> 29) & 0x1;    // positions 29
        ***/

        int mb1_phi  =   ( dataWordSub >> 30 ) & 0xFFF; // positions 30 -> 41
        int mb1_phib =   ( dataWordSub >> 42 ) & 0x3FF; // positions 42 -> 51
        int mb1_qual =   ( dataWordSub >> 52 ) & 0x7;   // positions 52 -> 54
        int mb1_ts2tag = ( dataWordSub >> 58 ) & 0x1;   // positions 58
        /*** NOT UNPACKED  
        int mb1_parity = ( dataWordSub >> 59 ) &0x1;    // positions 59
        ***/
        
        int mb1_phi_conv  = radAngConversion(mb1_phi);
        int mb1_phib_conv = benAngConversion(mb1_phib);

        int mb2_phi_conv  = radAngConversion(mb2_phi);
        int mb2_phib_conv = benAngConversion(mb2_phib);

        phiSegments.push_back( L1MuDTChambPhDigi( bxNr, wheel, sector-1, 
                                                1, mb1_phi_conv, mb1_phib_conv, 
                                                mb1_qual, mb1_ts2tag, bxCounter ) );
        phiSegments.push_back( L1MuDTChambPhDigi( bxNr, wheel, sector-1, 
                                                2, mb2_phi_conv, mb2_phib_conv, 
                                                mb2_qual, mb2_ts2tag, bxCounter ) );

        if ( debug_ ) logfile << '[' << ++lcounter << "]\t"<< std::hex 
                              << dataWordSub   << std::dec      << "\t|\t"
                              << "mb1_ts2tag " << mb1_ts2tag    << '\t'
                              << "mb1_qual "   << mb1_qual      << '\t'
                              << "mb1_phib "   << mb1_phib_conv << '\t'
                              << "mb1_phi "    << mb1_phi_conv  << '\t'
                              << "mb2_ts2tag " << mb2_ts2tag    << '\t'
                              << "mb2_qual "   << mb2_qual      << '\t'
                              << "mb2_phib "   << mb2_phib_conv << '\t'
                              << "mb2_phi "    << mb2_phi_conv  << '\n';
      }//MB1/2 word
   
      else if ( selector == 0x2 ) {

        int mb4_phi =    ( dataWordSub & 0xFFF);        // positions  0 -> 11
        int mb4_phib =   ( dataWordSub >> 12 ) & 0x3FF; // positions 12 -> 21
        int mb4_qual =   ( dataWordSub >> 22 ) & 0x7;   // positions 22 -> 24
        int mb4_ts2tag = ( dataWordSub >> 28 ) & 0x1;   // positions 28
        /*** NOT UNPACKED  
        int mb4_parity = ( dataWordSub >> 29) & 0x1;    // positions 29
        ***/

        int mb3_phi  =   ( dataWordSub >> 30 ) & 0xFFF; // positions 30 -> 41
        int mb3_phib =   ( dataWordSub >> 42 ) & 0x3FF; // positions 42 -> 51
        int mb3_qual =   ( dataWordSub >> 52 ) & 0x7;   // positions 52 -> 54
        int mb3_ts2tag = ( dataWordSub >> 58 ) & 0x1;   // positions 58
        /*** NOT UNPACKED  
        int mb3_parity = ( dataWordSub >> 59 ) &0x1;    // positions 59
        ***/

        int mb3_phi_conv  = radAngConversion(mb3_phi);
        int mb3_phib_conv = benAngConversion(mb3_phib);

        int mb4_phi_conv  = radAngConversion(mb4_phi);
        int mb4_phib_conv = benAngConversion(mb4_phib);

        phiSegments.push_back( L1MuDTChambPhDigi( bxNr, wheel, sector-1, 
                                                  3, mb3_phi_conv, mb3_phib_conv, 
                                                  mb3_qual, mb3_ts2tag, bxCounter) );
        phiSegments.push_back( L1MuDTChambPhDigi( bxNr, wheel, sector-1, 
                                                  4, mb4_phi_conv, mb4_phib_conv, 
                                                  mb4_qual, mb4_ts2tag, bxCounter) );

        if ( debug_ ) logfile << '[' << ++lcounter << "]\t"<< std::hex
                              << dataWordSub   << std::dec      << "\t|\t"
                              << "mb3_ts2tag " << mb3_ts2tag    << '\t'
                              << "mb3_qual "   << mb3_qual      << '\t'
                              << "mb3_phib "   << mb3_phib_conv << '\t'
                              << "mb3_phi "    << mb3_phi_conv  << '\t'
                              << "mb4_ts2tag " << mb4_ts2tag    << '\t'
                              << "mb4_qual "   << mb4_qual      << '\t'
                              << "mb4_phib "   << mb4_phib_conv << '\t'
                              << "mb4_phi "    << mb4_phi_conv  << '\n';

      }//MB3/4 word
  
      else if ( selector == 0x3 ) { //etha word
       
        int posALL, posBTI[7];
    
        int mb3_eta    = ( dataWordSub & 0xFF );        // positions  0 -> 7
        int mb3_eta_HQ = ( dataWordSub >> 8  ) & 0xFF;  // positions  8 -> 15
        int mb2_eta    = ( dataWordSub >> 16 ) & 0xFF;  // positions 16 -> 23

        int mb2_eta_HQ = ( dataWordSub >> 32 ) & 0xFF;  // positions 32 -> 39
        int mb1_eta    = ( dataWordSub >> 40 ) & 0xFF;  // positions 40 -> 47
        int mb1_eta_HQ = ( dataWordSub >> 48 ) & 0xFF;  // positions 48 -> 55

        if ( debug_ ) logfile << '[' << ++lcounter << "]\t" << std::hex
                              << dataWordSub << std::dec << "\t|\t"
                              << "mb1_eta_HQ " <<  mb1_eta_HQ   << '\t'
                              << "mb1_eta " << mb1_eta  << '\t'
                              << "mb2_eta_HQ " <<  mb2_eta_HQ   << '\t'
                              << "mb2_eta " << mb2_eta  << '\t'
                              << "mb3_eta_HQ " <<  mb3_eta_HQ   << '\t'
                              << "mb3_eta " << mb3_eta  << '\n';
        
        posALL    = mb1_eta_HQ & 0x7F;
        posBTI[0] = mb1_eta_HQ & 0x01;
        posBTI[1] = (mb1_eta_HQ & 0x02)>>1;
        posBTI[2] = (mb1_eta_HQ & 0x04)>>2;
        posBTI[3] = (mb1_eta_HQ & 0x08)>>3;
        posBTI[4] = (mb1_eta_HQ & 0x10)>>4;
        posBTI[5] = (mb1_eta_HQ & 0x20)>>5;
        posBTI[6] = (mb1_eta_HQ & 0x40)>>6;

        if ( posALL ) {
            
          theSegments.push_back( L1MuDTChambThDigi( bxNr, wheel, sector-1, 1, posBTI) );
            
        }
        
        posALL    = mb2_eta_HQ & 0x7F;
        posBTI[0] = mb2_eta_HQ & 0x01;
        posBTI[1] = (mb2_eta_HQ & 0x02)>>1;
        posBTI[2] = (mb2_eta_HQ & 0x04)>>2;
        posBTI[3] = (mb2_eta_HQ & 0x08)>>3;
        posBTI[4] = (mb2_eta_HQ & 0x10)>>4;
        posBTI[5] = (mb2_eta_HQ & 0x20)>>5;
        posBTI[6] = (mb2_eta_HQ & 0x40)>>6;

        if ( posALL ) {
            
          theSegments.push_back( L1MuDTChambThDigi( bxNr, wheel, sector-1, 2, posBTI) );
            
        }

        posALL    = mb3_eta_HQ & 0x7F;
        posBTI[0] = mb3_eta_HQ & 0x01;
        posBTI[1] = (mb3_eta_HQ & 0x02)>>1;
        posBTI[2] = (mb3_eta_HQ & 0x04)>>2;
        posBTI[3] = (mb3_eta_HQ & 0x08)>>3;
        posBTI[4] = (mb3_eta_HQ & 0x10)>>4;
        posBTI[5] = (mb3_eta_HQ & 0x20)>>5;
        posBTI[6] = (mb3_eta_HQ & 0x40)>>6;

        if ( posALL ) {
            
          theSegments.push_back( L1MuDTChambThDigi( bxNr, wheel, sector-1, 3, posBTI) );
            
        }
    
      }//etha word
     
      else if ( selector == 0xB ) { //MB1/2 output word
 
        /*** NOT UNPACKED  
        int mb2_phi =    ( dataWordSub & 0xFFF);        // positions  0 -> 11
        int mb2_phib =   ( dataWordSub >> 12 ) & 0x3FF; // positions 12 -> 21
        int mb2_qual =   ( dataWordSub >> 22 ) & 0x7;   // positions 22 -> 24
        int mb2_q3 =     ( dataWordSub >> 25 ) & 0x1;   // positions 25
        int mb2_q4 =     ( dataWordSub >> 26 ) & 0x1;   // positions 26
        int mb2_ts2tag = ( dataWordSub >> 28 ) & 0x1;   // positions 28
        int mb2_parity = ( dataWordSub >> 29) & 0x1;    // positions 29

        int mb1_phi  =   ( dataWordSub >> 30 ) & 0xFFF; // positions 30 -> 41
        int mb1_phib =   ( dataWordSub >> 42 ) & 0x3FF; // positions 42 -> 51
        int mb1_qual =   ( dataWordSub >> 52 ) & 0x7;   // positions 52 -> 54
        int mb1_q3 =     ( dataWordSub >> 55 ) & 0x1;   // positions 55
        int mb1_q4 =     ( dataWordSub >> 56 ) & 0x1;   // positions 56
        int mb1_ts2tag = ( dataWordSub >> 58 ) & 0x1;   // positions 58
        int mb1_parity = ( dataWordSub >> 59 ) &0x1;    // positions 59

        if ( debug_ ) logfile << '[' << ++lcounter << "]\t"<< std::hex 
                              << dataWordSub   << std::dec      << "\t|\t"
                              << "mb1_ts2tag_out " << mb1_ts2tag    << '\t'
                              << "mb1_qual_out "   << mb1_qual      << '\t'
                              << "mb1_q3_out "     << mb1_q3        << '\t'
                              << "mb1_q4_out "     << mb1_q4        << '\t'
                              << "mb1_phib_out "   << mb1_phib_conv << '\t'
                              << "mb1_phi_out "    << mb1_phi_conv  << '\t'
                              << "mb2_ts2tag_out " << mb2_ts2tag    << '\t'
                              << "mb2_qual_out "   << mb2_qual      << '\t'
                              << "mb2_q3_out "     << mb2_q3        << '\t'
                              << "mb2_q4_out "     << mb2_q4        << '\t'
                              << "mb2_phib_out "   << mb2_phib_conv << '\t'
                              << "mb2_phi_out "    << mb2_phi_conv  << '\n';
        ***/

        if ( debug_ ) logfile << '[' << ++lcounter << "]\t" << std::hex 
                              << dataWordSub << std::dec
                              << "\t MB1/2 OUTPUT WORD\n";

      }//MB1/2 output word
 
 
      else if ( selector == 0xC ) { //MB3/4 output word
 
        /*** NOT UNPACKED  
        int mb4_phi =    ( dataWordSub & 0xFFF);        // positions  0 -> 11
        int mb4_phib =   ( dataWordSub >> 12 ) & 0x3FF; // positions 12 -> 21
        int mb4_qual =   ( dataWordSub >> 22 ) & 0x7;   // positions 22 -> 24
        int mb4_q3 =     ( dataWordSub >> 25 ) & 0x1;   // positions 25
        int mb4_q4 =     ( dataWordSub >> 26 ) & 0x1;   // positions 26
        int mb4_ts2tag = ( dataWordSub >> 28 ) & 0x1;   // positions 28
        int mb4_parity = ( dataWordSub >> 29) & 0x1;    // positions 29

        int mb3_phi  =   ( dataWordSub >> 30 ) & 0xFFF; // positions 30 -> 41
        int mb3_phib =   ( dataWordSub >> 42 ) & 0x3FF; // positions 42 -> 51
        int mb3_qual =   ( dataWordSub >> 52 ) & 0x7;   // positions 52 -> 54
        int mb3_q3 =     ( dataWordSub >> 55 ) & 0x1;   // positions 55
        int mb3_q4 =     ( dataWordSub >> 56 ) & 0x1;   // positions 56
        int mb3_ts2tag = ( dataWordSub >> 58 ) & 0x1;   // positions 58
        int mb3_parity = ( dataWordSub >> 59 ) &0x1;    // positions 59

        if ( debug_ ) logfile << '[' << ++lcounter << "]\t"<< std::hex 
                              << dataWordSub   << std::dec      << "\t|\t"
                              << "mb3_ts2tag_out " << mb3_ts2tag    << '\t'
                              << "mb3_qual_out "   << mb3_qual      << '\t'
                              << "mb3_q3_out "     << mb3_q3        << '\t'
                              << "mb3_q4_out "     << mb3_q4        << '\t'
                              << "mb3_phib_out "   << mb3_phib_conv << '\t'
                              << "mb3_phi_out "    << mb3_phi_conv  << '\t'
                              << "mb4_ts2tag_out " << mb4_ts2tag    << '\t'
                              << "mb4_qual_out "   << mb4_qual      << '\t'
                              << "mb4_q3_out "     << mb4_q3        << '\t'
                              << "mb4_q4_out "     << mb4_q4        << '\t'
                              << "mb4_phib_out "   << mb4_phib_conv << '\t'
                              << "mb4_phi_out "    << mb4_phi_conv  << '\n';
        ***/

        if ( debug_ ) logfile << '[' << ++lcounter << "]\t" << std::hex 
                              << dataWordSub << std::dec
                              << "\t MB3/4 OUTPUT WORD\n";

      }//MB3/4 output word
 
      else if ( selector == 0xD ) { //etha output word
         
        if ( debug_ ) logfile << '[' << ++lcounter << "]\t" << std::hex 
                              << dataWordSub << std::dec
                              << "\t ETHA OUTPUT WORD\n";
    
      }//etha output word
       
      else if ( selector == 0x9 || selector == 0xE ) { //RPC word
          
        edm::LogInfo("TwinMux_unpacker") << "RPC WORD [" << std::dec << tm7eventsize << "] : "
                                         << std::hex << dataWordSub << std::dec
                                         << " it pos " << int(DTTM7iterator - DTTM7itend);
  
        if ( debug_ ) logfile << '[' << ++lcounter << "]\t" << std::hex 
                              << dataWordSub << std::dec
                              << "\t RPC WORD\n";          

      }//RPC word
      
      else if ( selector == 0x6 ) { //HO word
	
        edm::LogInfo("TwinMux_unpacker") << "HO WORD [" << std::dec << tm7eventsize << "] : "
                                         << std::hex << dataWordSub << std::dec
                                         << "1 it pos " << int(DTTM7iterator - DTTM7itend);
	
        if ( debug_ ) logfile << '[' << ++lcounter << "]\t" << std::hex 
                              << dataWordSub << std::dec
                              << "\t HO WORD\n";          


        if(dodebug) std::cout<<"TwinMux_unpacker::HO WORD [" << std::dec << tm7eventsize << "] : "
                                         << std::hex << dataWordSub <<  " it pos " << int(DTTM7iterator - DTTM7itend)<<std::endl;
	
        if ( dodebug ) std::cout << '[' << ++lcounter << "]\t" << std::hex 
                              << dataWordSub << std::dec
                              << "\t HO WORD\n";          
	
	

	
	std::ostringstream ss;
	for (unsigned int i=0; i<hofedUnpackList_.size(); i++)
	  ss << hofedUnpackList_[i] << " ";
	edm::LogInfo("HO") << "HO TP will unpack FEDs ( " << ss.str() << ")";

	if(dodebug) std::cout<<"HO TP will unpack FEDs ( " << ss.str() << ")"<<std::endl;	


	int top_mip =   ( dataWordSub >> 23 );          // positions  0 -> 23
        // int bot_mip =   ( dataWordSub >> 29 ) & 0xFF;   // positions 29 -> 52
        //int top_bx  =   ( dataWordSub >> 24 ) & 0xFF;   // positions 24 -> 27
        //int bot_bx  =   ( dataWordSub >> 53 ) & 0xFF;   // positions 53 -> 56
	
	std::cout<<"mip is : "<< top_mip<<std::endl;

	// unpack all HO FEDs                                       
	HOUnrolledTP unrolled[24];

	for (std::vector<int>::const_iterator i=hofedUnpackList_.begin(); i!=hofedUnpackList_.end(); i++) {
	  const FEDRawData& fed = data->FEDData(*i);
	  if(dodebug)	  std::cout<<"fed.size() "<< fed.size() << std::endl;
	  if (fed.size()==0) {
	    if (complainEmptyData_) {
	      if (!silent_) edm::LogWarning("EmptyData") << "No data for FED " << *i;
	      if(dodebug) std::cout<<"Empty data for fed"<<*i << std::endl;
	      //	      report->addError(*i);
	    }
	  } else if (fed.size() < 8*3) {
	    if (!silent_) edm::LogWarning("EmptyData") << "Tiny data " << fed.size() << " for FED " << *i;
	    //	    report->addError(*i);
	  } else {
	    try {

	      //main code begin here
	      // get the DCC header
	      const HcalDCCHeader* dccHeader=(const HcalDCCHeader*)(fed.data());
	      const HcalDTCHeader* dtcHeader=(const HcalDTCHeader*)(fed.data());
	      bool is_VME_DCC=(dccHeader->getDCCDataFormatVersion()<0x10) || ((mode_&0x1)==0);
	      
	      int dccid=(is_VME_DCC)?(dccHeader->getSourceId()-sourceIdOffset_):(dtcHeader->getSourceId()-sourceIdOffset_);
	      if(dodebug) std::cout<<"********************** dccid is :"<< dccid<< std::endl;

	      // check the summary status
	      
	      // walk through the HTR data.  For the uTCA, use spigot=slot+1
	      HcalHTRData htr;
	      const unsigned short* daq_first, *daq_last, *tp_first, *tp_last;
	      //  const HcalQIESample* qie_begin, *qie_end, *qie_work;
	      const HcalTriggerPrimitiveSample *tp_begin, *tp_end, *tp_work; 
	      for (int spigot=0; spigot<HcalDCCHeader::SPIGOT_COUNT; spigot++) {
		if(dodebug)		std::cout<<"spigot is :"<< spigot << std::endl;

		if (is_VME_DCC) {
		  if (!dccHeader->getSpigotPresent(spigot)) continue;
		  
		  int retval=dccHeader->getSpigotData(spigot,htr,fed.size());
		  if (retval!=0) {
		    if (retval==-1) {
		      if (!silent) edm::LogWarning("Invalid Data") << "Invalid HTR data (data beyond payload size) observed on spigot " << spigot << " of DCC with source id " << dccHeader->getSourceId();
		      //		      report.countSpigotFormatError();
		    }
		    continue;
		  }
		  // check
		  if (dccHeader->getSpigotCRCError(spigot)) {
		    if (!silent) 
		      edm::LogWarning("Invalid Data") << "CRC Error on HTR data observed on spigot " << spigot << " of DCC with source id " << dccHeader->getSourceId();
		    //		    report.countSpigotFormatError();
		    continue;
		  } 
		} else { // is_uTCA (!is_VME_DCC)
		  int slot=spigot+1;
		  if (slot>HcalDTCHeader::MAXIMUM_SLOT) continue;
		  
		  if (!dtcHeader->getSlotPresent(slot)) continue;
		  
		  int retval=dtcHeader->getSlotData(slot,htr,fed.size());
		  if (retval!=0) {
		    if (retval==-1) {
		      if (!silent) edm::LogWarning("Invalid Data") << "Invalid uHTR data (data beyond payload size) observed on slot " << slot << " of DTC with source id " << dtcHeader->getSourceId();
		      //		      report.countSpigotFormatError();
		    }
		    continue;
		  }
		  
		  // check
		  if (dtcHeader->getSlotCRCError(slot)) {
		    if (!silent) 
		      edm::LogWarning("Invalid Data") << "CRC Error on uHTR data observed on slot " << slot << " of DTC with source id " << dtcHeader->getSourceId();
		    //		    report.countSpigotFormatError();
		    continue;
		  } 
		}
		
		
		// check for EE
		if (htr.isEmptyEvent()) {
		  //		  report.countEmptyEventSpigot();
		}
		if (htr.isOverflowWarning()) {
		  //		  report.countOFWSpigot();
		}
		if (htr.isBusy()) {
		  //		  report.countBusySpigot();
		}
		if (!htr.check()) {
		  if (!silent) 
		    edm::LogWarning("Invalid Data") << "Invalid HTR data observed on spigot " << spigot << " of DCC with source id " << dccHeader->getSourceId();
		  //		  report.countSpigotFormatError();
		  continue;
		}  
		
		if (htr.getFirmwareFlavor()>=0x80) {
		  if (!silent) edm::LogWarning("HOTPUnpacker") << "Skipping data on spigot " << spigot << " of DCC with source id " << dccHeader->getSourceId() << " which is of unknown flavor " << htr.getFirmwareFlavor();
		  continue;
		}
		
		// calculate "real" number of presamples
		// int nps=htr.getNPS()-startSample_;
		
		// get pointers
		htr.dataPointers(&daq_first,&daq_last,&tp_first,&tp_last);
		unsigned int smid = htr.getSubmodule();
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
		
		
		/*
		  Unpack the trigger primitives
		*/
		if (isHOtpg) {
		  if(dodebug) std::cout <<"hotp is being obserbed, will unpack"<< std::endl;

		  for (tp_work=tp_begin; tp_work!=tp_end; tp_work++) {
		    if (tp_work->raw()==0xFFFF) continue; // filler word
		    //int sector=tp_work->slbChan();
		    int sector=slbChan(tp_work->raw());
		    if(dodebug)	 std::cout<<"sector for the HO is: " << sector << std::endl;

		    if (sector>2) continue;
		    
		    for (int ibit=0; ibit<8; ibit++) {
		      int linear=sector*8+ibit; 
		      if (!unrolled[linear].checked) {
			unrolled[linear].checked=true;
			int fiber=(linear/3)+1;
			int fc=(linear%3);
			
			// electronics id (use precision match for HO TP)
			//http://cmslxr.fnal.gov/source/DataFormats/HcalDetId/interface/HcalElectronicsIdh.
			HcalElectronicsId eid(fc,fiber,spigot,dccid);	
			eid.setHTR(htr_cr,htr_slot,htr_tb);
			std::cout<<"electronics Id: "<< eid() << std::endl;
	
			//http://cmslxr.fnal.gov/source/CondFormats/HcalObjects/interface/HcalElectronicsMap.h
			DetId did= readoutMap->lookup(eid);
			std::cout<<"sector: "<<sector <<" subdet : "<< ((HcalSubdetector)did.subdetId()) << " htrchan:"<<eid.htrChanId()<<" htrSlot :"<<eid.htrSlot()<< " crateId() : "<<eid.crateId()<< std::endl;
			std::cout<<"detId :"<< did() <<" det: " << did.det()<<" subdetId : "<< did.subdetId() << std::endl;
			
			if (!did.null()) {
			  if (did.det()==DetId::Hcal && ((HcalSubdetector)did.subdetId())==HcalOuter ) {
			    HcalDetId hid(did);
			    unrolled[linear].valid=true;
			    unrolled[linear].ieta=hid.ieta();
			    unrolled[linear].iphi=hid.iphi();
			    if(dodebug)   std::cout<<"eta, phi is :"<< hid.ieta() <<"," << hid.iphi() << std::endl;
			  }
			}
		      }
		      if (unrolled[linear].valid) {
			if (isTPGSOI(*tp_work)) unrolled[linear].soi=unrolled[linear].samples;
			if (tp_work->raw()&(1<<ibit)) unrolled[linear].setbit(unrolled[linear].samples);
			unrolled[linear].samples++;
		      }
		    } // for (int ibit=0; ibit<8; ibit++)
		  } // for (tp_work=tp_begin; tp_work!=tp_end; tp_work++)
		  
		  for (int i=0; i<24; i++) {
		    if (unrolled[i].valid) 
		      colls.tphoCont->push_back(HOTriggerPrimitiveDigi(
								       unrolled[i].ieta,
								       unrolled[i].iphi,
								       unrolled[i].samples,
								       unrolled[i].soi,
								       unrolled[i].databits)
						);
		  }
		}
	      }
	      
	      //    unpack(fed, *readoutMap, colls, *report, silent_);
	      
	      
	      // main code ends here
	      //	      report->addUnpacked(*i);
	    } catch (cms::Exception& e) {
	      if (!silent_) edm::LogWarning("Unpacking error") << e.what();
	      //	      report->addError(*i);
	    } catch (...) {
	      if (!silent_) edm::LogWarning("Unpacking exception");
	      //	      report->addError(*i);
	    }
	  }
	} //  for (std::vector<int>::const_iterator
	
	
	// if ( debug_ ) logfile << '[' << ++lcounter << "]\t"<< std::hex
	// 		      << dataWordSub   << std::dec      << "\t|\t"
	// 		      << "mb3_ts2tag " << mb3_ts2tag    << '\t'
	// 		      << "mb3_qual "   << mb3_qual      << '\t'
	// 		      << "mb3_phib "   << mb3_phib_conv << '\t'
	// 		      << "mb3_phi "    << mb3_phi_conv  << '\t'
	// 		      << "mb4_ts2tag " << mb4_ts2tag    << '\t'
	// 		      << "mb4_qual "   << mb4_qual      << '\t'
	// 		      << "mb4_phib "   << mb4_phib_conv << '\t'
	// 		      << "mb4_phi "    << mb4_phi_conv  << '\n';
	
	
      }//HO word
      
      else if ( selector == 0xF ) { //ERROR word

        edm::LogInfo("TwinMux_unpacker") << "ERROR WORD [" << std::dec << tm7eventsize << "] : "
                                         << std::hex << dataWordSub << std::dec
                                         << " it pos " << int(DTTM7iterator - DTTM7itend);
  
        if ( debug_ ) logfile << '[' << ++lcounter << "]\t" << std::hex 
                              << dataWordSub << std::dec
                              << "\t ERROR WORD\n";
      }//ERROR word

      else { //unkown word

        edm::LogInfo("TwinMux_unpacker") << "UNKNOWN WORD received " << std::hex << dataWordSub 
                                           << " in FED " << std::hex << TM7fedId;

   	    if ( debug_ ) logfile << '[' << ++lcounter << "]\t" << std::hex 
                              << dataWordSub << std::dec
                              << "\t UNKNOWN WORD\n";
      }
  
      if( DTTM7iterator == DTTM7itend ) break;
      
    } //end of loop over AMCsize


    /// Trailer of payload with CRC
    ++DTTM7iterator;

    if( DTTM7iterator == DTTM7itend ) break;

  } // end for-loop container content

  return;
}



void L1TTwinMuxRawToDigi::calcCRC( long word, int & myC ) {

  int myCRC[16], D[64], C[16];

  for ( int i = 0; i < 64; ++i ) { D[i]    = (word >> i) & 0x1; }
  for ( int i = 0; i < 16; ++i ) { C[i]    = (myC>>i)  & 0x1; }

  myCRC[0] = ( D[63] + D[62] + D[61] + D[60] + D[55] + D[54] +
               D[53] + D[52] + D[51] + D[50] + D[49] + D[48] +
               D[47] + D[46] + D[45] + D[43] + D[41] + D[40] +
               D[39] + D[38] + D[37] + D[36] + D[35] + D[34] +
               D[33] + D[32] + D[31] + D[30] + D[27] + D[26] +
               D[25] + D[24] + D[23] + D[22] + D[21] + D[20] +
               D[19] + D[18] + D[17] + D[16] + D[15] + D[13] +
               D[12] + D[11] + D[10] + D[9]  + D[8]  + D[7]  +
               D[6]  + D[5]  + D[4]  + D[3]  + D[2]  + D[1]  +
               D[0]  + C[0]  + C[1]  + C[2]  + C[3]  + C[4]  +
               C[5]  + C[6]  + C[7]  + C[12] + C[13] + C[14] +
               C[15] )%2;

  myCRC[1] = ( D[63] + D[62] + D[61] + D[56] + D[55] + D[54] +
	           D[53] + D[52] + D[51] + D[50] + D[49] + D[48] +
	           D[47] + D[46] + D[44] + D[42] + D[41] + D[40] +
	           D[39] + D[38] + D[37] + D[36] + D[35] + D[34] +
    	       D[33] + D[32] + D[31] + D[28] + D[27] + D[26] +
    	       D[25] + D[24] + D[23] + D[22] + D[21] + D[20] +
    	       D[19] + D[18] + D[17] + D[16] + D[14] + D[13] +
    	       D[12] + D[11] + D[10] + D[9]  + D[8]  + D[7]  +
	           D[6]  + D[5]  + D[4]  + D[3]  + D[2]  + D[1]  +
	           C[0]  + C[1]  + C[2]  + C[3]  + C[4]  + C[5]  +
	           C[6]  + C[7]  + C[8]  + C[13] + C[14] + C[15] )%2;

  myCRC[2] = ( D[61] + D[60] + D[57] + D[56] + D[46] + D[42] +
	           D[31] + D[30] + D[29] + D[28] + D[16] + D[14] +
	           D[1]  + D[0]  + C[8]  + C[9]  + C[12] + C[13] )%2;

  myCRC[3] = ( D[62] + D[61] + D[58] + D[57] + D[47] + D[43] +
	           D[32] + D[31] + D[30] + D[29] + D[17] + D[15] +
	           D[2]  + D[1]  + C[9]  + C[10] + C[13] + C[14] )%2;

  myCRC[4] = ( D[63] + D[62] + D[59] + D[58] + D[48] + D[44] +
    	       D[33] + D[32] + D[31] + D[30] + D[18] + D[16] + 
	           D[3]  + D[2]  + C[0]  + C[10] + C[11] + C[14] +
	           C[15] )%2;

  myCRC[5] = ( D[63] + D[60] + D[59] + D[49] + D[45] + D[34] +
	           D[33] + D[32] + D[31] + D[19] + D[17] + D[4]  +
    	       D[3]  + C[1]  + C[11] + C[12] + C[15] )%2;

  myCRC[6] = ( D[61] + D[60] + D[50] + D[46] + D[35] + D[34] +
	           D[33] + D[32] + D[20] + D[18] + D[5]  + D[4]  +
	           C[2]  + C[12] + C[13] )%2;

  myCRC[7] = ( D[62] + D[61] + D[51] + D[47] + D[36] + D[35] +
    	       D[34] + D[33] + D[21] + D[19] + D[6]  + D[5]  +
	           C[3]  + C[13] + C[14] )%2;

  myCRC[8] = ( D[63] + D[62] + D[52] + D[48] + D[37] + D[36] +
	           D[35] + D[34] + D[22] + D[20] + D[7]  + D[6]  +
    	       C[0]  + C[4]  + C[14] + C[15] )%2;

  myCRC[9] = ( D[63] + D[53] + D[49] + D[38] + D[37] + D[36] +
	           D[35] + D[23] + D[21] + D[8]  + D[7]  + C[1]  +
	           C[5]  + C[15] )%2;

  myCRC[10] = ( D[54] + D[50] + D[39] + D[38] + D[37] + D[36] + 
       		    D[24] + D[22] + D[9]  + D[8]  + C[2]  + C[6] )%2;

  myCRC[11] = ( D[55] + D[51] + D[40] + D[39] + D[38] + D[37] +
		        D[25] + D[23] + D[10] + D[9]  + C[3]  + C[7] )%2;

  myCRC[12] = ( D[56] + D[52] + D[41] + D[40] + D[39] + D[38] +
        		D[26] + D[24] + D[11] + D[10] + C[4]  + C[8] )%2;

  myCRC[13] = ( D[57] + D[53] + D[42] + D[41] + D[40] + D[39] +
		        D[27] + D[25] + D[12] + D[11] + C[5]  + C[9] )%2;

  myCRC[14] = ( D[58] + D[54] + D[43] + D[42] + D[41] + D[40] +
        		D[28] + D[26] + D[13] + D[12] + C[6]  + C[10] )%2;

  myCRC[15] = ( D[63] + D[62] + D[61] + D[60] + D[59] + D[54] +
		        D[53] + D[52] + D[51] + D[50] + D[49] + D[48] + 
	        	D[47] + D[46] + D[45] + D[44] + D[42] + D[40] +
        		D[39] + D[38] + D[37] + D[36] + D[35] + D[34] + 
		        D[33] + D[32] + D[31] + D[30] + D[29] + D[26] +
        		D[25] + D[24] + D[23] + D[22] + D[21] + D[20] + 
        		D[19] + D[18] + D[17] + D[16] + D[15] + D[14] +
        		D[12] + D[11] + D[10] + D[9]  + D[8]  + D[7]  + 
        		D[6]  + D[5]  + D[4]  + D[3]  + D[2]  + D[1]  +
        		D[0]  + C[0]  + C[1]  + C[2]  + C[3]  + C[4]  + 
	        	C[5]  + C[6]  + C[11] + C[12] + C[13] + C[14] +
	        	C[15] )%2;

  int tempC = 0x0;  
  for ( int i = 0; i < 16 ; ++i) { tempC = tempC + ( myCRC[i] << i ); }
  myC = tempC;
  return;
}


L1TTwinMuxRawToDigi::Collections::Collections() {
  tphoCont=0;
}


//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TTwinMuxRawToDigi);
