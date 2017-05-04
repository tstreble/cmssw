///
/// \class l1t::Stage2Layer2SumsAlgorithmFirmwareImp1
///
/// \author:
///
/// Description:

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "L1Trigger/L1TCalorimeter/interface/Stage2Layer2DemuxSumsAlgoFirmware.h"

#include "L1Trigger/L1TCalorimeter/interface/CaloParamsHelper.h"

#include <vector>
#include <algorithm>


l1t::Stage2Layer2DemuxSumsAlgoFirmwareImp1::Stage2Layer2DemuxSumsAlgoFirmwareImp1(CaloParamsHelper* params) :
  params_(params), cordic_(Cordic(144*16,17,8))  // These are the settings in the hardware - should probably make this configurable
{
}


l1t::Stage2Layer2DemuxSumsAlgoFirmwareImp1::~Stage2Layer2DemuxSumsAlgoFirmwareImp1() {


}


void l1t::Stage2Layer2DemuxSumsAlgoFirmwareImp1::processEvent(const std::vector<l1t::EtSum> & inputSums,
                                                              std::vector<l1t::EtSum> & outputSums) {


  // Add up the x, y and scalar components
  for (std::vector<l1t::EtSum>::const_iterator eSum = inputSums.begin() ; eSum != inputSums.end() ; ++eSum )
    {
      switch (eSum->getType()) {

      case l1t::EtSum::EtSumType::kTotalEt:
        et_ += eSum->hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalEtEm:
        etem_ += eSum->hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalEtx:
        metx_ += eSum->hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalEty:
        mety_ += eSum->hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalHt:
        ht_ += eSum->hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalHtx:
        mhtx_ += eSum->hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalHty:
        mhty_ += eSum->hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalEtxHF:
        metxHF_ += eSum->hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalEtyHF:
        metyHF_ += eSum->hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalHtxHF:
        mhtxHF_ += eSum->hwPt();
        break;
	
      case l1t::EtSum::EtSumType::kTotalHtyHF:
        mhtyHF_ += eSum->hwPt();
        break;

      case l1t::EtSum::EtSumType::kMinBiasHFP0:
	mbp0_ = eSum->hwPt();
	break;

      case l1t::EtSum::EtSumType::kMinBiasHFM0:
	mbm0_ = eSum->hwPt();
	break;

      case l1t::EtSum::EtSumType::kMinBiasHFP1:
	mbp1_ = eSum->hwPt();
	break;

      case l1t::EtSum::EtSumType::kMinBiasHFM1:
	mbm1_ = eSum->hwPt();
	break;

      case l1t::EtSum::EtSumType::kTowerCount:
	ntow_ = eSum->hwPt();
	break;

      default:
        continue; // Should throw an exception or something?
      }
    }


  calibrateEnergySums();
  
  if (et_>0xFFF)   et_   = 0xFFF;
  //if (metx_>0xFFF) metx_ = 0xFFF;
  //if (mety_>0xFFF) mety_ = 0xFFF;
  if (ht_>0xFFF)   ht_   = 0xFFF;
  //if (mhtx_>0xFFF) mhtx_ = 0xFFF;
  //if (mhty_>0xFFF) mhty_ = 0xFFF;
  //if (metxHF_>0xFFF) metxHF_ = 0xFFF;
  //if (metyHF_>0xFFF) metyHF_ = 0xFFF;

  mhtPhi_ = (111 << 4);
  mhtPhiHF_ = (111 << 4); // to match hw value if undefined
  
  // Final MET calculation
  if (metx_ != 0 || mety_ != 0 ) cordic_( metx_ , mety_ , metPhi_ , met_ );
  // sets the met scale back to the original range for output into GT, this corresponds to
  // the previous scaling of sin/cos factors in calculation of metx and mety by 2^10 = 1024
  met_ >>= 10; 

  // Final METHF calculation
  if (metxHF_ != 0 || metyHF_ != 0 ) cordic_( metxHF_ , metyHF_ , metPhiHF_ , metHF_ );
  metHF_ >>= 10;


  // Final MHT calculation
  if (mhtx_ != 0 || mhty_ != 0 ) cordic_( mhtx_ , mhty_ , mhtPhi_ , mht_ );
  // sets the mht scale back to the original range for output into GT, the other 4
  // bits are brought back just before the accumulation of ring sum in MP jet sum algorithm
  mht_ >>= 6; 

  if (mhtxHF_ != 0 || mhtyHF_ != 0 ) cordic_( mhtxHF_ , mhtyHF_ , mhtPhiHF_ , mhtHF_ );
  mhtHF_ >>= 6; 

  // Make final collection
  math::XYZTLorentzVector p4;

  l1t::EtSum etSumTotalEt(p4,l1t::EtSum::EtSumType::kTotalEt,et_,0,0,0);
  l1t::EtSum etSumTotalEtEm(p4,l1t::EtSum::EtSumType::kTotalEtEm,etem_,0,0,0);
  l1t::EtSum etSumMissingEt(p4,l1t::EtSum::EtSumType::kMissingEt,met_,0,metPhi_>>4,0);
  l1t::EtSum etSumMissingEtHF(p4,l1t::EtSum::EtSumType::kMissingEtHF,metHF_,0,metPhiHF_>>4,0);
  l1t::EtSum htSumht(p4,l1t::EtSum::EtSumType::kTotalHt,ht_,0,0,0);
  l1t::EtSum htSumMissingHt(p4,l1t::EtSum::EtSumType::kMissingHt,mht_,0,mhtPhi_>>4,0);
  l1t::EtSum htSumMissingHtHF(p4,l1t::EtSum::EtSumType::kMissingHtHF,mhtHF_,0,mhtPhiHF_>>4,0);
  l1t::EtSum etSumMinBiasHFP0(p4,l1t::EtSum::EtSumType::kMinBiasHFP0,mbp0_,0,0,0);
  l1t::EtSum etSumMinBiasHFM0(p4,l1t::EtSum::EtSumType::kMinBiasHFM0,mbm0_,0,0,0);
  l1t::EtSum etSumMinBiasHFP1(p4,l1t::EtSum::EtSumType::kMinBiasHFP1,mbp1_,0,0,0);
  l1t::EtSum etSumMinBiasHFM1(p4,l1t::EtSum::EtSumType::kMinBiasHFM1,mbm1_,0,0,0);
  l1t::EtSum etSumTowCount(p4,l1t::EtSum::EtSumType::kTowerCount,ntow_,0,0,0);

  l1t::EtSum etSumTotalEtx(p4, l1t::EtSum::EtSumType::kTotalEtx,metx_,0,0,0);
  l1t::EtSum etSumTotalEty(p4, l1t::EtSum::EtSumType::kTotalEty,mety_,0,0,0);

  outputSums.push_back(etSumTotalEt);
  outputSums.push_back(etSumTotalEtEm);
  outputSums.push_back(etSumMinBiasHFP0);
  outputSums.push_back(htSumht);
  outputSums.push_back(etSumMinBiasHFM0);
  outputSums.push_back(etSumMissingEt);
  outputSums.push_back(etSumMinBiasHFP1);
  outputSums.push_back(htSumMissingHt);
  outputSums.push_back(etSumMinBiasHFM1);
  outputSums.push_back(etSumMissingEtHF);
  outputSums.push_back(htSumMissingHtHF);
  outputSums.push_back(etSumTowCount);
  
  //also add MP Etx & Ety for calibration derivation
  outputSums.push_back(etSumTotalEtx);
  outputSums.push_back(etSumTotalEty);

  resetEnergySums();

}



void l1t::Stage2Layer2DemuxSumsAlgoFirmwareImp1::calibrateEnergySums()
{

  if(params_->etSumXCalibrationType() == "LUT"){
    //calibrate Ex
    uint exCalibLUTAddr = metx_ >> 5;
    uint exMult = ( 0x3FF & params_->etSumXCalibrationLUT()->data(exCalibLUTAddr) );
    int8_t exAdd  = ( 0xFF | ( params_->etSumXCalibrationLUT()->data(exCalibLUTAddr) >> 10 ) ); 
    uint exCorr = (  ( metx_ * exMult ) >> 9 ) + exAdd;
    metx_ = exCorr;
    
    //calibrate ExHF
    uint exHFCalibLUTAddr = metxHF_ >> 5;
    uint exHFMult = ( 0x3FF & params_->etSumXCalibrationLUT()->data(exHFCalibLUTAddr) );
    int8_t exHFAdd  = ( 0xFF | ( params_->etSumXCalibrationLUT()->data(exHFCalibLUTAddr) >> 10 ) ); 
    uint exHFCorr = (  ( metxHF_ * exHFMult ) >> 9 )  + exHFAdd;
    metxHF_ = exHFCorr;
    
  } else {
      if(params_->etSumXCalibrationType() != "None" && params_->etSumXCalibrationType() != "none") 
	edm::LogError("l1t|stage 2") << "Invalid etSumX calibration type in calo params. Not applying etSumX calibration to Stage 2 ETT & MET" << std::endl;
      return;
  }
  
  if(params_->etSumYCalibrationType() == "LUT"){
    //calibrate Ey
    uint eyCalibLUTAddr = mety_ >> 5;
    uint eyMult = ( 0x3FF & params_->etSumYCalibrationLUT()->data(eyCalibLUTAddr) );
    int8_t eyAdd  = ( 0xFF | ( params_->etSumYCalibrationLUT()->data(eyCalibLUTAddr) >> 10 ) ); 
    uint eyCorr = (  ( mety_ * eyMult ) >> 9 )  + eyAdd;
    mety_ = eyCorr;
    
    //calibrate EyHF
    uint eyHFCalibLUTAddr = metyHF_ >> 5;
    uint eyHFMult = ( 0x3FF & params_->etSumYCalibrationLUT()->data(eyHFCalibLUTAddr) );
    int8_t eyHFAdd  = ( 0xFF | ( params_->etSumYCalibrationLUT()->data(eyHFCalibLUTAddr) >> 10 ) ); 
    uint eyHFCorr = (  ( metyHF_ * eyHFMult ) >> 9 )  + eyHFAdd;
      metyHF_ = eyHFCorr;
      
  } else {
    if(params_->etSumYCalibrationType() != "None" && params_->etSumYCalibrationType() != "none") 
      edm::LogError("l1t|stage 2") << "Invalid etSumY calibration type in calo params. Not applying etSumY calibration to Stage 2 ETT & MET" << std::endl;
    return;
  }


  if(params_->etSumEttCalibrationType() == "LUT"){
    //calibrate Et
    uint etCalibLUTAddr = et_ >> 5;
    uint etMult = ( 0x3FF & params_->etSumEttCalibrationLUT()->data(etCalibLUTAddr) );
    int8_t etAdd  = ( 0xFF | ( params_->etSumEttCalibrationLUT()->data(etCalibLUTAddr) >> 10 ) ); 
    uint etCorr = (  ( et_ * etMult ) >> 9 )  + etAdd;
    et_ = etCorr;
    
  } else {
    if(params_->etSumEttCalibrationType() != "None" && params_->etSumEttCalibrationType() != "none") 
      edm::LogError("l1t|stage 2") << "Invalid etSumEtt calibration type in calo params. Not applying etSumEtt calibration to Stage 2 ETT & MET" << std::endl;
    return;
  }
  
  if(params_->etSumEcalSumCalibrationType() == "LUT"){
    //calibrate Etem
    uint etemCalibLUTAddr = etem_ >> 5;
    uint etemMult = ( 0x3FF & params_->etSumEcalSumCalibrationLUT()->data(etemCalibLUTAddr) );
    int8_t etemAdd  = ( 0xFF | ( params_->etSumEcalSumCalibrationLUT()->data(etemCalibLUTAddr) >> 10 ) ); 
    uint etemCorr = (  ( etem_ * etemMult ) >> 9 )  + etemAdd;
    etem_ = etemCorr;
    
  } else {
    if(params_->etSumEcalSumCalibrationType() != "None" && params_->etSumEcalSumCalibrationType() != "none") 
      edm::LogError("l1t|stage 2") << "Invalid etSumEcalSum calibration type in calo params. Not applying etSumEcalSum calibration to Stage 2 ETT & MET" << std::endl;
    return;
  }
  
   
}



void l1t::Stage2Layer2DemuxSumsAlgoFirmwareImp1::resetEnergySums()
{
  
  et_ = 0; 
  etem_ = 0; 
  metx_ = 0; 
  mety_ = 0; 
  metxHF_ = 0; 
  metyHF_ = 0;
  ht_ = 0; 
  mhtx_ = 0; 
  mhty_ = 0; 
  mhtxHF_ = 0; 
  mhtyHF_ = 0; 
  metPhi_ = 0; 
  metPhiHF_ = 0; 
  mhtPhi_ = 0; 
  mhtPhiHF_ = 0;
  met_ = 0; 
  metHF_ = 0; 
  mht_ = 0; 
  mhtHF_ = 0;
  mbp0_ = 0; 
  mbm0_ = 0; 
  mbp1_ = 0; 
  mbm1_ = 0;
  ntow_ = 0;
  
}
