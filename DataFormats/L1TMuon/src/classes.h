#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/RefToBase.h"

#include "DataFormats/L1TMuon/interface/MuonCaloSumFwd.h"
#include "DataFormats/L1TMuon/interface/MuonCaloSum.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/EMTFOutput.h"
#include "DataFormats/L1TMuon/interface/EMTFTrack.h"
#include "DataFormats/L1TMuon/interface/EMTFHit.h"

namespace {
  struct dictionary {
    l1t::MuonCaloSumBxCollection caloSum;
    edm::Wrapper<l1t::MuonCaloSumBxCollection> caloSumWrap;
    std::vector<l1t::MuonCaloSum> vCaloSum;

    l1t::RegionalMuonCandBxCollection regCand;
    edm::Wrapper<l1t::RegionalMuonCandBxCollection> regCandWrap;
    std::vector<l1t::RegionalMuonCand> vRegCand;
   
    l1t::EMTFOutputCollection emtfOutput;
    edm::Wrapper<l1t::EMTFOutputCollection> emtfOutputWrap;
   
    l1t::EMTFTrackCollection emtfTrack;
    edm::Wrapper<l1t::EMTFTrackCollection> emtfTrackWrap;
   
    l1t::EMTFHitCollection emtfHit;
    edm::Wrapper<l1t::EMTFHitCollection> emtfHitWrap;
   
  };
}


