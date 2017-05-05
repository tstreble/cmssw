///
/// Description: Firmware headers
///
/// Implementation:
///    Concrete firmware implementations
///
/// \author: Jim Brooke - University of Bristol
///

//
//

#ifndef Stage2Layer2DemuxSumsAlgoFirmware_H
#define Stage2Layer2DemuxSumsAlgoFirmware_H

#include "L1Trigger/L1TCalorimeter/interface/Stage2Layer2DemuxSumsAlgo.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloParamsHelper.h"
#include "L1Trigger/L1TCalorimeter/interface/Cordic.h"

namespace l1t {

  // Imp1 is for v1 and v2
  class Stage2Layer2DemuxSumsAlgoFirmwareImp1 : public Stage2Layer2DemuxSumsAlgo {
  public:
    Stage2Layer2DemuxSumsAlgoFirmwareImp1(CaloParamsHelper* params);
    virtual ~Stage2Layer2DemuxSumsAlgoFirmwareImp1();
    virtual void processEvent(const std::vector<l1t::EtSum> & inputSums,
			      std::vector<l1t::EtSum> & outputSums);
    virtual void calibrateEnergySums();
    virtual void resetEnergySums();
  private:

    CaloParamsHelper* params_;

    Cordic cordic_;

    int32_t et_, etem_, metx_, mety_, metxHF_, metyHF_;
    int32_t ht_, mhtx_, mhty_, mhtxHF_, mhtyHF_; 
    int32_t metPhi_, metPhiHF_, mhtPhi_, mhtPhiHF_;
    uint32_t met_, metHF_, mht_, mhtHF_;
    uint32_t mbp0_, mbm0_, mbp1_, mbm1_;
    uint32_t ntow_;

  };

}

#endif
