#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

//
// class declaration
//

using namespace std;

class PFCandProducer : public edm::EDProducer {
    
public:
    
    explicit PFCandProducer(const edm::ParameterSet &iConfig);
    
    ~PFCandProducer() override {};
    
    
private:

    virtual void produce(edm::Event&, const edm::EventSetup&);

    pair<double,double> computeDCA(const pat::PackedCandidate &pfCand,
				   edm::ESHandle<MagneticField> bFieldHandle,
                                   reco::BeamSpot beamSpot);

    // ----------member data ---------------------------
    
    edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
    edm::EDGetTokenT<edm::View<pat::PackedCandidate>> PFCandSrc_;
    edm::EDGetTokenT<edm::View<pat::IsolatedTrack>> IsoTrackSrc_;

};

PFCandProducer::PFCandProducer(const edm::ParameterSet &iConfig):
beamSpotSrc_( consumes<reco::BeamSpot> ( iConfig.getParameter<edm::InputTag>( "beamSpot" ) ) ),
PFCandSrc_( consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "PFCandCollection" ) ) ),
IsoTrackSrc_( consumes<edm::View<pat::IsolatedTrack>> ( iConfig.getParameter<edm::InputTag>( "IsoTrackCollection" ) ) )
{
    produces<pat::CompositeCandidateCollection>();
}



void PFCandProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
    edm::ESHandle<MagneticField> bFieldHandle;
    edm::Handle<reco::BeamSpot> beamSpotHandle;
    
    iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
    iEvent.getByToken(beamSpotSrc_, beamSpotHandle);
    
    if ( ! beamSpotHandle.isValid() ) {
        edm::LogError("PFCandProducer") << "No beam spot available from EventSetup" ;
    }

    reco::BeamSpot beamSpot = *beamSpotHandle;
    
    edm::Handle<edm::View<pat::PackedCandidate>> pfCandHandle;
    
    iEvent.getByToken(PFCandSrc_, pfCandHandle);
    
    unsigned int pfCandNumber = pfCandHandle->size();

    edm::Handle<edm::View<pat::IsolatedTrack>> isoTrackHandle;

    iEvent.getByToken(IsoTrackSrc_, isoTrackHandle);

    unsigned int isoTrackNumber = isoTrackHandle->size();
    
    // Output collection
    std::unique_ptr<pat::CompositeCandidateCollection> result( new pat::CompositeCandidateCollection );

    for (unsigned int i = 0; i < pfCandNumber; ++i) {            

      const pat::PackedCandidate & pfCand = (*pfCandHandle)[i];
      auto pfCandPtr = edm::Ptr<pat::PackedCandidate>(pfCandHandle,i);

      double DCABS = -1.;
      double DCABSErr = -1.;

      if(pfCand.hasTrackDetails()) {

	pair<double,double> DCA = computeDCA(pfCand,
					     bFieldHandle,
					     beamSpot);
	 
	DCABS = DCA.first;
	DCABSErr = DCA.second;

      }

      double dEdXStrip = -1.;
      double dEdXPixel = -1.;

      for (unsigned int j = 0; j < isoTrackNumber; ++j) {

	const pat::IsolatedTrack & isoTrack = (*isoTrackHandle)[j];
	auto pfCand_fromTrack_Ptr = edm::refToPtr(isoTrack.packedCandRef());

	if(isoTrack.packedCandRef().isNonnull() && pfCand_fromTrack_Ptr==pfCandPtr){
	  dEdXStrip = isoTrack.dEdxStrip();
	  dEdXPixel = isoTrack.dEdxPixel();
	}

      }

      pat::CompositeCandidate pfCandNew;
      pfCandNew.addDaughter( pfCand );
      pfCandNew.addUserFloat("DCASig", DCABS/DCABSErr);
      pfCandNew.addUserFloat("dEdXStrip", dEdXStrip);
      pfCandNew.addUserFloat("dEdXPixel", dEdXPixel);

      result->push_back(pfCandNew);

    }

    iEvent.put(std::move(result));

}




pair<double,double> PFCandProducer::computeDCA(const pat::PackedCandidate &pfCand,
					       edm::ESHandle<MagneticField> bFieldHandle,
					       reco::BeamSpot beamSpot){
    
    const reco::TransientTrack trackTT((*(pfCand.bestTrack())), &(*bFieldHandle));
    
    TrajectoryStateClosestToPoint theDCAXBS = trackTT.trajectoryStateClosestToPoint( GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()) );
    
    double DCABS = theDCAXBS.perigeeParameters().transverseImpactParameter();
    double DCABSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
    
    pair<double,double> DCA = make_pair(DCABS,DCABSErr);
    
    return DCA;
}







DEFINE_FWK_MODULE(PFCandProducer);
