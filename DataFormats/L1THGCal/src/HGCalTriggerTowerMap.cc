#include "DataFormats/L1THGCal/interface/HGCalTriggerTowerMap.h"

using namespace l1t;

namespace l1t {
  bool operator < ( const l1t::HGCalTriggerCell& a, const l1t::HGCalTriggerCell & b )
  {
      return a.hwPt() < b.hwPt();
  }
}


void HGCalTriggerTowerMap::processTriggerCells(){

  std::sort (triggerCells_.rbegin(), triggerCells_.rend());
  
  int hwPt = 0; 
  for(auto& TC : triggerCells_) hwPt += TC.hwPt();
  setHwPt(hwPt);

}
