#ifndef DataFormats_L1TCalorimeter_HGCalTriggerTowerMap_h
#define DataFormats_L1TCalorimeter_HGCalTriggerTowerMap_h


#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/L1THGCal/interface/HGCalTriggerCell.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"

namespace l1t 
{

    class HGCalTriggerTowerMap;
    typedef BXVector<HGCalTriggerTowerMap> HGCalTriggerTowerMapBxCollection;

    class HGCalTriggerTowerMap : public HGCalTriggerCell
    {

        public:

            HGCalTriggerTowerMap() {}

            ~HGCalTriggerTowerMap() {}

	    void addTriggerCell( HGCalTriggerCell TC );
	    std::vector<HGCalTriggerCell> triggerCells() { return triggerCells_; }
	    HGCalTriggerCell maxTriggerCell() { return maxTriggerCell_; }
	    
          
        private:
  
	    std::vector<HGCalTriggerCell> triggerCells_;
	    HGCalTriggerCell maxTriggerCell_;

    };

}

#endif
