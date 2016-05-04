///
/// \class l1t::Stage2Layer2EtSumAlgorithmFirmwareImp1
///
/// \author: Jim Brooke
///
/// Description: first iteration of stage 2 jet algo

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "L1Trigger/L1TCalorimeter/interface/Stage2Layer2EtSumAlgorithmFirmware.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"
#include <math.h>


l1t::Stage2Layer2EtSumAlgorithmFirmwareImp1::Stage2Layer2EtSumAlgorithmFirmwareImp1(CaloParamsHelper* params) :
  params_(params)
{

  // Add some LogDebug for these settings

  metTowThresholdHw_ = floor(params_->etSumEtThreshold(0)/params_->towerLsbSum());
  ettTowThresholdHw_ = floor(params_->etSumEtThreshold(2)/params_->towerLsbSum());

  metEtaMax_ = params_->etSumEtaMax(0);
  ettEtaMax_ = params_->etSumEtaMax(2);
}


l1t::Stage2Layer2EtSumAlgorithmFirmwareImp1::~Stage2Layer2EtSumAlgorithmFirmwareImp1() {}

void l1t::Stage2Layer2EtSumAlgorithmFirmwareImp1::processEvent(const std::vector<l1t::CaloTower> & towers,
                                                               std::vector<l1t::EtSum> & etsums) {

  
  // etaSide=1 is positive eta, etaSide=-1 is negative eta
  for (int etaSide=1; etaSide>=-1; etaSide-=2) {

    int32_t ex(0), ey(0), et(0);

    for (unsigned absieta=1; absieta<CaloTools::kHFEnd; absieta++) {

      int ieta = etaSide * absieta;

      // TODO add the eta and Et thresholds

      int32_t ringEx(0), ringEy(0), ringEt(0);

      for (int iphi=1; iphi<=CaloTools::kHBHENrPhi; iphi++) {
      
        l1t::CaloTower tower = l1t::CaloTools::getTower(towers, ieta, iphi);

		if (tower.hwPt()>metTowThresholdHw_ && CaloTools::mpEta(abs(tower.hwEta()))<=metEtaMax_) {
		  
		  // x- and -y coefficients are truncated by after multiplication of Et by trig coefficient.
		  // The trig coefficients themselves take values [-1023,1023] and so were scaled by
		  // 2^10 = 1024, which requires bitwise shift to the right of the final value by 10 bits.
		  // This is accounted for at ouput of demux (see Stage2Layer2DemuxSumsAlgoFirmwareImp1.cc)
		  ringEx += (int32_t) (tower.hwPt() * CaloTools::cos_coeff[iphi - 1] );
		  ringEy += (int32_t) (tower.hwPt() * CaloTools::sin_coeff[iphi - 1] );

		}
		if (tower.hwPt()>ettTowThresholdHw_ && CaloTools::mpEta(abs(tower.hwEta()))<=ettEtaMax_) 
		  ringEt += tower.hwPt();
      }    
      
      ex += ringEx;
      ey += ringEy;
      et += ringEt;
    }

    math::XYZTLorentzVector p4;

    l1t::EtSum etSumTotalEt(p4,l1t::EtSum::EtSumType::kTotalEt,et,0,0,0);
    l1t::EtSum etSumEx(p4,l1t::EtSum::EtSumType::kTotalEtx,ex,0,0,0);
    l1t::EtSum etSumEy(p4,l1t::EtSum::EtSumType::kTotalEty,ey,0,0,0);

    etsums.push_back(etSumTotalEt);
    etsums.push_back(etSumEx);
    etsums.push_back(etSumEy);

  }

}
