#ifndef __L1Trigger_L1THGCal_HGCalTriggerCellTowerMapCodec_h__
#define __L1Trigger_L1THGCal_HGCalTriggerCellTowerMapCodec_h__

#include "L1Trigger/L1THGCal/interface/HGCalTriggerFECodecBase.h"
#include "L1Trigger/L1THGCal/interface/fe_codecs/HGCalTriggerCellTowerMapCodecImpl.h"


inline std::ostream& operator<<(std::ostream& o, const HGCalTriggerCellTowerMapDataPayload& data) 
{ 
    for(const auto& dat : data.payload)
    {
        o << "(" << std::hex << dat.detId() 
            << std::dec << " " << dat.hwPt() << ") ";
    }
    o << "\n";
    return o;
}


class HGCalTriggerCellTowerMapCodec : public HGCalTriggerFE::Codec<HGCalTriggerCellTowerMapCodec,HGCalTriggerCellTowerMapDataPayload> 
{
    public:
        typedef HGCalTriggerCellTowerMapDataPayload data_type;

        HGCalTriggerCellTowerMapCodec(const edm::ParameterSet& conf);

        void setDataPayloadImpl(const HGCEEDigiCollection& ee,
                const HGCHEDigiCollection& fh,
                const HGCBHDigiCollection& bh );

        void setDataPayloadImpl(const l1t::HGCFETriggerDigi& digi);

        std::vector<bool> encodeImpl(const data_type&) const ;
        data_type         decodeImpl(const std::vector<bool>&, const uint32_t) const;  

    private:
        HGCalTriggerCellTowerMapCodecImpl codecImpl_;
};

#endif
