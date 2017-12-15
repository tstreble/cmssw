#include "DataFormats/L1THGCal/interface/HGCalTriggerTowerMap.h"

using namespace l1t;


void HGCalTriggerTowerMap::addTriggerCell( l1t::HGCalTriggerCell TC ){
  
  this->setHwPt(this->hwPt()+TC.hwPt());
  if(TC.hwPt()>maxTriggerCell_.hwPt()) maxTriggerCell_ = TC;
  triggerCells_.emplace_back(TC);

}
