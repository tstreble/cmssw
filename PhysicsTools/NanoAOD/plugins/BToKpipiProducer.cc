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
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include <TLorentzVector.h>
#include "DataFormats/Math/interface/Vector3D.h"
#include <TVector.h>
#include <TMatrix.h>



//
// class declaration
//

using namespace std;

class BToKpipiProducer : public edm::EDProducer {
    
public:
    
    explicit BToKpipiProducer(const edm::ParameterSet &iConfig);
    
    ~BToKpipiProducer() override {};
    
    
private:
    
    virtual void produce(edm::Event&, const edm::EventSetup&);
    
    bool KPiVertexRefitting(const pat::PackedCandidate &kaon,
                            const pat::PackedCandidate &pi,
                            edm::ESHandle<MagneticField> bFieldHandle,
                            RefCountedKinematicVertex &refitVertex,
                            RefCountedKinematicParticle &refitKPi,
                            RefCountedKinematicParticle &refitKaon,
                            RefCountedKinematicParticle &refitPi_D0);
    
    bool BToKPiPiVertexRefitting(const pat::PackedCandidate &kaon,
                                 const pat::PackedCandidate &pi_D0,
                                 const pat::PackedCandidate &pi_Bu,
                                 edm::ESHandle<MagneticField> bFieldHandle,
                                 RefCountedKinematicVertex &refitVertex,
                                 RefCountedKinematicParticle &refitBToKPiPi,
                                 RefCountedKinematicParticle &refitKaon,
                                 RefCountedKinematicParticle &refitPi_D0,
                                 RefCountedKinematicParticle &refitPi_Bu);

    bool BToD0PiVertexRefitting(const RefCountedKinematicParticle refitKPi,
				const pat::PackedCandidate &pi_Bu,
				edm::ESHandle<MagneticField> bFieldHandle,
				RefCountedKinematicVertex &refitVertex,
				RefCountedKinematicParticle &refitBToKPiPi,
				RefCountedKinematicParticle &refitD0,
				RefCountedKinematicParticle &refitPi_Bu);
    
    pair<double,double> computeLS(RefCountedKinematicVertex refitVertex,
                           reco::BeamSpot beamSpot);
    
    double computeCosAlpha(RefCountedKinematicParticle refitBToKPiPi,
                           RefCountedKinematicVertex vertexFitTree,
                           reco::BeamSpot beamSpot);
    
    pair<double,double> computeDCA(const pat::PackedCandidate &kaon,
                                   edm::ESHandle<MagneticField> bFieldHandle,
                                   reco::BeamSpot beamSpot);
    
    
    // ----------member data ---------------------------
    
    edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
    edm::EDGetTokenT<edm::View<pat::PackedCandidate>> PFCandSrc_;
    
    double ptMin_;
    double ptMinPiB_;
    double etaMax_;
    double DCASigMin_;
    double massMinKPi_;
    double massMaxKPi_;
    double CLVtxMinKPi_;
    bool D0Charge_;
    bool KCharge_;
    double massMinB_;
    double massMaxB_;
    double CLVtxMinB_;
    double D0MassConstraint_;
    bool save3TrkRefit_;

    float KaonMass_ = 0.493677;
    float KaonMassErr_ = 1.6e-5;
    float PionMass_ = 0.139570;
    float PionMassErr_ = 3.5e-7;
    //float D0Mass_ = 1.86484;  //Configurable parameter
    float D0MassErr_ = 0.05;
    
};



BToKpipiProducer::BToKpipiProducer(const edm::ParameterSet &iConfig):
beamSpotSrc_( consumes<reco::BeamSpot> ( iConfig.getParameter<edm::InputTag>( "beamSpot" ) ) ),
PFCandSrc_( consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "PFCandCollection" ) ) ),
ptMin_( iConfig.getParameter<double>( "MinPt" ) ),
ptMinPiB_( iConfig.getParameter<double>( "MinPtPiB" ) ),
etaMax_( iConfig.getParameter<double>( "MaxEta" ) ),
DCASigMin_( iConfig.getParameter<double>( "MinDCASig") ),
massMinKPi_( iConfig.getParameter<double>( "D0MinMass" ) ),
massMaxKPi_( iConfig.getParameter<double>( "D0MaxMass" ) ),
CLVtxMinKPi_(iConfig.getParameter<double>( "D0MinCLVtx" ) ),
D0Charge_( iConfig.getParameter<bool>( "D0ChargeCheck" ) ),
KCharge_( iConfig.getParameter<bool>( "KaonChargeCheck" ) ),
massMinB_( iConfig.getParameter<double>( "BMinMass" ) ),
massMaxB_( iConfig.getParameter<double>( "BMaxMass" ) ),
CLVtxMinB_( iConfig.getParameter<double>( "BMinCLVtx" ) ),
D0MassConstraint_( iConfig.getParameter<double>( "D0MassConstraint" ) ),
save3TrkRefit_( iConfig.getParameter<bool>( "save3TrackRefit" ) )
{
    produces<pat::CompositeCandidateCollection>();
}


void BToKpipiProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
    edm::ESHandle<MagneticField> bFieldHandle;
    edm::Handle<reco::BeamSpot> beamSpotHandle;
    
    iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
    iEvent.getByToken(beamSpotSrc_, beamSpotHandle);
    
    if ( ! beamSpotHandle.isValid() ) {
        edm::LogError("BToKpipiProducer") << "No beam spot available from EventSetup" ;
    }
    
    reco::BeamSpot beamSpot = *beamSpotHandle;
    
    edm::Handle<edm::View<pat::PackedCandidate>> pfCandHandle;
    
    iEvent.getByToken(PFCandSrc_, pfCandHandle);
    
    unsigned int pfCandNumber = pfCandHandle->size();
    
    // Output collection
    std::unique_ptr<pat::CompositeCandidateCollection> result( new pat::CompositeCandidateCollection );
    
    
    if(pfCandNumber>1){        

      // loop on all the Kpipi triplets

      for (unsigned int i = 0; i < pfCandNumber; ++i) {

	const pat::PackedCandidate & kaon = (*pfCandHandle)[i];

	if(abs(kaon.pdgId())!=211) continue; //Charged hadrons
	if(!kaon.hasTrackDetails()) continue;
	if(kaon.pt()<ptMin_ || abs(kaon.eta())>etaMax_) continue;

	pair<double,double> kaon_DCA = computeDCA(kaon,
						  bFieldHandle,
						  beamSpot);
	double kaon_DCABS = kaon_DCA.first;
	double kaon_DCABSErr = kaon_DCA.second;

	if(fabs(kaon_DCABS/kaon_DCABSErr)<DCASigMin_) continue;

	for (unsigned int j = 0; j < pfCandNumber; ++j) {

	  if(j==i) continue;

	  const pat::PackedCandidate & piD0 = (*pfCandHandle)[j];

	  if(abs(piD0.pdgId())!=211) continue; //Charged hadrons
	  if(!piD0.hasTrackDetails()) continue;
	  if(piD0.pt()<ptMin_ || abs(piD0.eta())>etaMax_) continue;

	  if(D0Charge_ && kaon.charge()*piD0.charge()>0) continue;

	  pair<double,double> piD0_DCA = computeDCA(piD0,
						    bFieldHandle,
						    beamSpot);
	  double piD0_DCABS = piD0_DCA.first;
	  double piD0_DCABSErr = piD0_DCA.second;
	  if(fabs(piD0_DCABS/piD0_DCABSErr)<DCASigMin_) continue;

	  RefCountedKinematicVertex refitVertexKPi;
	  RefCountedKinematicParticle refitKPi;
	  RefCountedKinematicParticle refitKaon_KPi;
	  RefCountedKinematicParticle refitPiD0_KPi;

	  bool passed = KPiVertexRefitting(kaon, piD0,
					   bFieldHandle,
					   refitVertexKPi,
					   refitKPi,
					   refitKaon_KPi,
					   refitPiD0_KPi);

	  if ( !passed) continue;

	  pair<double,double> KPiLS = computeLS(refitVertexKPi,beamSpot);
	  double KPiLSBS = KPiLS.first;
	  double KPiLSBSErr = KPiLS.second;

	  double KPiVtx_CL = TMath::Prob((double)refitVertexKPi->chiSquared(),
					 int(rint(refitVertexKPi->degreesOfFreedom())));

	  if(KPiVtx_CL<CLVtxMinKPi_) continue;

	  double KPi_mass = refitKPi->currentState().mass();
	  if(KPi_mass<massMinKPi_ || KPi_mass>massMaxKPi_) continue;
  
	  double KPi_mass_err = sqrt(refitKPi->currentState().kinematicParametersError().matrix()(6,6));

	  math::XYZVector refitKaonV3D_KPi = refitKaon_KPi->refittedTransientTrack().track().momentum();
	  math::XYZVector refitPiD0V3D_KPi = refitPiD0_KPi->refittedTransientTrack().track().momentum();
	  math::XYZVector refitKPiV3D = refitKaonV3D_KPi + refitPiD0V3D_KPi;

	  for (unsigned int k = 0; k < pfCandNumber; ++k) {

	    if(k==i || k==j) continue;

	    const pat::PackedCandidate & piBu = (*pfCandHandle)[k];
	    if(abs(piBu.pdgId())!=211) continue; //Charged hadrons
	    if(!piBu.hasTrackDetails()) continue;
	    if(piBu.pt()<ptMinPiB_ || abs(piBu.eta())>etaMax_) continue;
	    if(KCharge_ && kaon.charge()*piBu.charge()<0) continue; //Removed suppressed decays

	    pair<double,double> piBu_DCA = computeDCA(piBu,
						      bFieldHandle,
						      beamSpot);
	    double piBu_DCABS = piBu_DCA.first;
	    double piBu_DCABSErr = piBu_DCA.second;

	    if(fabs(piBu_DCABS/piBu_DCABSErr)<DCASigMin_) continue;

	    RefCountedKinematicVertex refitVertexBToD0Pi;
	    RefCountedKinematicParticle refitBToD0Pi;
	    RefCountedKinematicParticle refitD0_D0Pi;
	    RefCountedKinematicParticle refitPiBu_D0Pi;

	    passed = BToD0PiVertexRefitting(refitKPi, piBu,
					    bFieldHandle,
					    refitVertexBToD0Pi,
					    refitBToD0Pi,
					    refitD0_D0Pi,
					    refitPiBu_D0Pi);

	    if (!passed) continue;

	    math::XYZVector refitD0V3D_D0Pi = refitD0_D0Pi->refittedTransientTrack().track().momentum();
	    math::XYZVector refitPiBuV3D_D0Pi = refitPiBu_D0Pi->refittedTransientTrack().track().momentum();
	    math::XYZVector refitBToD0PiV3D = refitD0V3D_D0Pi + refitPiBuV3D_D0Pi;

	    pair<double,double> LS = computeLS(refitVertexBToD0Pi,beamSpot);
	    double LSBS = LS.first;
	    double LSBSErr  = LS.second;

	    double BToD0PiVtx_CL = TMath::Prob((double)refitVertexBToD0Pi->chiSquared(),
					       int(rint(refitVertexBToD0Pi->degreesOfFreedom())));

            if(BToD0PiVtx_CL<CLVtxMinB_) continue;

            double cosAlpha = computeCosAlpha(refitBToD0Pi,refitVertexBToD0Pi,beamSpot);

            double mass = refitBToD0Pi->currentState().mass();
            if(mass<massMinB_ || mass>massMaxB_) continue;

            double mass_err = sqrt(refitBToD0Pi->currentState().kinematicParametersError().matrix()(6,6));

	    pat::CompositeCandidate BToKPiPiCand;
	    BToKPiPiCand.addDaughter( kaon , "kaon");
	    BToKPiPiCand.addDaughter( piD0 , "piD0");
	    BToKPiPiCand.addDaughter( piBu,  "piBu");

	    BToKPiPiCand.addUserFloat("kaon_pt",     sqrt(refitKaonV3D_KPi.perp2()));
	    BToKPiPiCand.addUserFloat("kaon_eta",    refitKaonV3D_KPi.eta());
	    BToKPiPiCand.addUserFloat("kaon_phi",    refitKaonV3D_KPi.phi());
	    BToKPiPiCand.addUserFloat("kaon_charge", refitKaon_KPi->currentState().particleCharge());
	    BToKPiPiCand.addUserFloat("kaon_DCASig", kaon_DCABS/kaon_DCABSErr);

	    BToKPiPiCand.addUserFloat("piD0_pt",     sqrt(refitPiD0V3D_KPi.perp2()));
	    BToKPiPiCand.addUserFloat("piD0_eta",    refitPiD0V3D_KPi.eta());
	    BToKPiPiCand.addUserFloat("piD0_phi",    refitPiD0V3D_KPi.phi());
	    BToKPiPiCand.addUserFloat("piD0_charge", refitPiD0_KPi->currentState().particleCharge());
	    BToKPiPiCand.addUserFloat("piD0_DCASig", piD0_DCABS/piD0_DCABSErr);

	    BToKPiPiCand.addUserFloat("piBu_pt",    sqrt(refitPiBuV3D_D0Pi.perp2()));
	    BToKPiPiCand.addUserFloat("piBu_eta",   refitPiBuV3D_D0Pi.eta());
	    BToKPiPiCand.addUserFloat("piBu_phi",   refitPiBuV3D_D0Pi.phi());
	    BToKPiPiCand.addUserFloat("piBu_charge",refitPiBu_D0Pi->currentState().particleCharge());
	    BToKPiPiCand.addUserFloat("piBu_DCASig", piBu_DCABS/piBu_DCABSErr);

	    BToKPiPiCand.addUserFloat("Kpi_pt",     sqrt(refitKPiV3D.perp2()));
	    BToKPiPiCand.addUserFloat("Kpi_eta",    refitKPiV3D.eta());
	    BToKPiPiCand.addUserFloat("Kpi_phi",    refitKPiV3D.phi());
	    BToKPiPiCand.addUserFloat("Kpi_mass",   KPi_mass);
	    BToKPiPiCand.addUserFloat("Kpi_mass_err",   KPi_mass_err);
	    BToKPiPiCand.addUserFloat("Kpi_Lxy", (float) KPiLSBS/KPiLSBSErr);
	    BToKPiPiCand.addUserFloat("Kpi_CL_vtx", (float) KPiVtx_CL);

	    BToKPiPiCand.addUserFloat("pt",     sqrt(refitBToD0PiV3D.perp2()));
	    BToKPiPiCand.addUserFloat("eta",    refitBToD0PiV3D.eta());
	    BToKPiPiCand.addUserFloat("phi",    refitBToD0PiV3D.phi());
	    BToKPiPiCand.addUserFloat("mass",   mass);
	    BToKPiPiCand.addUserFloat("mass_err", mass_err);
	    BToKPiPiCand.addUserFloat("Lxy", (float) LSBS/LSBSErr);
	    BToKPiPiCand.addUserFloat("CL_vtx", (float) BToD0PiVtx_CL);
	    BToKPiPiCand.addUserFloat("cosAlpha", (float) cosAlpha);

	    //Optional 3-track refitting variables

	    float pt_3trk = -9999.;
	    float eta_3trk = -9999.;
	    float phi_3trk = -9999.;
	    float mass_3trk = -9999.;
	    float mass_err_3trk = -9999.;
	    float Lxy_3trk = -9999.;
	    float CL_vtx_3trk = -9999.;
	    float cosAlpha_3trk = -9999.;

	    if(save3TrkRefit_){

	      RefCountedKinematicVertex refitVertexBToKPiPi;
	      RefCountedKinematicParticle refitBToKPiPi;
	      RefCountedKinematicParticle refitKaon_KPiPi;
	      RefCountedKinematicParticle refitPiD0_KPiPi;
	      RefCountedKinematicParticle refitPiBu_KPiPi;

	      passed = BToKPiPiVertexRefitting(kaon, piD0, piBu,
					       bFieldHandle,
					       refitVertexBToKPiPi,
					       refitBToKPiPi,
					       refitKaon_KPiPi,
					       refitPiD0_KPiPi,
					       refitPiBu_KPiPi);

	      if (passed){

		pair<double,double> BToKPiPiLS = computeLS(refitVertexBToKPiPi,beamSpot);
		double LSBS_3trk = BToKPiPiLS.first;
		double LSBSErr_3trk = BToKPiPiLS.second;

		math::XYZVector refitKaonV3D_KPiPi = refitKaon_KPiPi->refittedTransientTrack().track().momentum();
		math::XYZVector refitPiD0V3D_KPiPi = refitPiD0_KPiPi->refittedTransientTrack().track().momentum();
		math::XYZVector refitPiBuV3D_KPiPi = refitPiBu_KPiPi->refittedTransientTrack().track().momentum();
		math::XYZVector refitBToKPiPiV3D = refitKaonV3D_KPiPi + refitPiD0V3D_KPiPi + refitPiBuV3D_KPiPi;

		pt_3trk = sqrt(refitBToKPiPiV3D.perp2());
		eta_3trk = refitBToKPiPiV3D.eta();
		phi_3trk = refitBToKPiPiV3D.phi();
		mass_3trk = refitBToKPiPi->currentState().mass();
		mass_err_3trk = sqrt(refitBToKPiPi->currentState().kinematicParametersError().matrix()(6,6));
		Lxy_3trk = (float) LSBS_3trk/LSBSErr_3trk;
		CL_vtx_3trk = TMath::Prob((double)refitVertexBToKPiPi->chiSquared(),
					  int(rint(refitVertexBToKPiPi->degreesOfFreedom())));
		cosAlpha_3trk = computeCosAlpha(refitBToKPiPi,refitVertexBToKPiPi,beamSpot);

	      }

	    }

	    BToKPiPiCand.addUserFloat("pt_3trk", pt_3trk);
	    BToKPiPiCand.addUserFloat("eta_3trk", eta_3trk);
	    BToKPiPiCand.addUserFloat("phi_3trk", phi_3trk);
	    BToKPiPiCand.addUserFloat("mass_3trk", mass_3trk);
	    BToKPiPiCand.addUserFloat("mass_err_3trk", mass_err_3trk);
	    BToKPiPiCand.addUserFloat("Lxy_3trk", Lxy_3trk);
	    BToKPiPiCand.addUserFloat("CL_vtx_3trk", CL_vtx_3trk);
	    BToKPiPiCand.addUserFloat("cosAlpha_3trk", cosAlpha_3trk);

	    result->push_back(BToKPiPiCand);

	  }

	}
        
      }

    }
    
    iEvent.put(std::move(result));
    
}



bool BToKpipiProducer::KPiVertexRefitting(const pat::PackedCandidate &kaon,
                                          const pat::PackedCandidate &pi,
                                          edm::ESHandle<MagneticField> bFieldHandle,
                                          RefCountedKinematicVertex &refitVertex,
                                          RefCountedKinematicParticle &refitKPi,
                                          RefCountedKinematicParticle &refitKaon,
                                          RefCountedKinematicParticle &refitPi_D0){
    
    const reco::TransientTrack kaonTT(*(kaon.bestTrack()), &(*bFieldHandle));
    const reco::TransientTrack piTT(*(pi.bestTrack()), &(*bFieldHandle));

    KinematicParticleFactoryFromTransientTrack partFactory;
    KinematicParticleVertexFitter PartVtxFitter;

    std::vector<RefCountedKinematicParticle> KPiParticles;
    double chi = 0.;
    double ndf = 0.;
    KPiParticles.push_back(partFactory.particle(kaonTT,KaonMass_,chi,ndf,KaonMassErr_));
    KPiParticles.push_back(partFactory.particle(piTT,PionMass_,chi,ndf,PionMassErr_));
    RefCountedKinematicTree KPiVertexFitTree = PartVtxFitter.fit(KPiParticles);

    if ( !KPiVertexFitTree->isValid()) return false;

    KPiVertexFitTree->movePointerToTheTop();
    refitVertex = KPiVertexFitTree->currentDecayVertex();
    refitKPi = KPiVertexFitTree->currentParticle();

    if ( !refitVertex->vertexIsValid()) return false;

    // extract the re-fitted tracks
    KPiVertexFitTree->movePointerToTheTop();

    KPiVertexFitTree->movePointerToTheFirstChild();
    refitKaon = KPiVertexFitTree->currentParticle();

    KPiVertexFitTree->movePointerToTheNextChild();
    refitPi_D0 = KPiVertexFitTree->currentParticle();
    
    return true;
    
}




bool BToKpipiProducer::BToKPiPiVertexRefitting(const pat::PackedCandidate &kaon,
                                               const pat::PackedCandidate &pi_D0,
                                               const pat::PackedCandidate &pi_Bu,
                                               edm::ESHandle<MagneticField> bFieldHandle,
                                               RefCountedKinematicVertex &refitVertex,
                                               RefCountedKinematicParticle &refitBToKPiPi,
                                               RefCountedKinematicParticle &refitKaon,
                                               RefCountedKinematicParticle &refitPi_D0,
                                               RefCountedKinematicParticle &refitPi_Bu){

    const reco::TransientTrack kaonTT(*(kaon.bestTrack()), &(*bFieldHandle));
    const reco::TransientTrack piD0TT(*(pi_D0.bestTrack()), &(*bFieldHandle));
    const reco::TransientTrack piBuTT(*(pi_Bu.bestTrack()), &(*bFieldHandle));

    KinematicParticleFactoryFromTransientTrack partFactory;
    KinematicParticleVertexFitter PartVtxFitter;

    std::vector<RefCountedKinematicParticle> BToKPiPiParticles;
    double chi = 0.;
    double ndf = 0.;
    BToKPiPiParticles.push_back(partFactory.particle(kaonTT,KaonMass_,chi,ndf,KaonMassErr_));
    BToKPiPiParticles.push_back(partFactory.particle(piD0TT,PionMass_,chi,ndf,PionMassErr_));
    BToKPiPiParticles.push_back(partFactory.particle(piBuTT,PionMass_,chi,ndf,PionMassErr_));

    RefCountedKinematicTree BToKPiPiVertexFitTree = PartVtxFitter.fit(BToKPiPiParticles);

    if ( !BToKPiPiVertexFitTree->isValid()) return false;

    BToKPiPiVertexFitTree->movePointerToTheTop();
    refitVertex = BToKPiPiVertexFitTree->currentDecayVertex();
    refitBToKPiPi = BToKPiPiVertexFitTree->currentParticle();

    if ( !refitVertex->vertexIsValid()) return false;

    // extract the re-fitted tracks
    BToKPiPiVertexFitTree->movePointerToTheTop();

    BToKPiPiVertexFitTree->movePointerToTheFirstChild();
    refitKaon = BToKPiPiVertexFitTree->currentParticle();

    BToKPiPiVertexFitTree->movePointerToTheNextChild();
    refitPi_D0 = BToKPiPiVertexFitTree->currentParticle();

    BToKPiPiVertexFitTree->movePointerToTheNextChild();
    refitPi_Bu = BToKPiPiVertexFitTree->currentParticle();

    return true;

}







bool BToKpipiProducer::BToD0PiVertexRefitting(const RefCountedKinematicParticle refitKPi,
					      const pat::PackedCandidate &pi_Bu,
					      edm::ESHandle<MagneticField> bFieldHandle,
					      RefCountedKinematicVertex &refitVertex,
					      RefCountedKinematicParticle &refitBToKPiPi,
					      RefCountedKinematicParticle &refitD0,
					      RefCountedKinematicParticle &refitPi_Bu){

  const reco::TransientTrack D0TT = refitKPi->refittedTransientTrack();
  const reco::TransientTrack piBuTT(*(pi_Bu.bestTrack()), &(*bFieldHandle));

  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter PartVtxFitter;

  std::vector<RefCountedKinematicParticle> BToKPiPiParticles;
  double chi = 0.;
  double ndf = 0.;

  float KPi_mass = refitKPi->currentState().mass();
  float KPi_mass_err = sqrt(refitKPi->currentState().kinematicParametersError().matrix()(6,6));
  if(D0MassConstraint_ > 0){
    KPi_mass = D0MassConstraint_;
    KPi_mass_err = D0MassErr_;
  }

  BToKPiPiParticles.push_back(partFactory.particle(D0TT,KPi_mass,chi,ndf,KPi_mass_err));
  BToKPiPiParticles.push_back(partFactory.particle(piBuTT,PionMass_,chi,ndf,PionMassErr_));

  RefCountedKinematicTree BToKPiPiVertexFitTree = PartVtxFitter.fit(BToKPiPiParticles);

  if ( !BToKPiPiVertexFitTree->isValid()) return false;

  BToKPiPiVertexFitTree->movePointerToTheTop();
  refitVertex = BToKPiPiVertexFitTree->currentDecayVertex();
  refitBToKPiPi = BToKPiPiVertexFitTree->currentParticle();

  if ( !refitVertex->vertexIsValid()) return false;

  // extract the re-fitted tracks
  BToKPiPiVertexFitTree->movePointerToTheTop();

  BToKPiPiVertexFitTree->movePointerToTheFirstChild();
  refitD0 = BToKPiPiVertexFitTree->currentParticle();

  BToKPiPiVertexFitTree->movePointerToTheNextChild();
  refitPi_Bu = BToKPiPiVertexFitTree->currentParticle();

  return true;

}






pair<double,double> BToKpipiProducer::computeLS(RefCountedKinematicVertex refitVertex,
                                         reco::BeamSpot beamSpot){
    
    TVector v(2);
    v[0] = refitVertex->position().x()-beamSpot.position().x();
    v[1] = refitVertex->position().y()-beamSpot.position().y();
    
    TMatrix errVtx(2,2);
    errVtx(0,0) = refitVertex->error().cxx();
    errVtx(0,1) = refitVertex->error().matrix()(0,1);
    errVtx(1,0) = errVtx(0,1);
    errVtx(1,1) = refitVertex->error().cyy();
    
    TMatrix errBS(2,2);
    errBS(0,0) = beamSpot.covariance()(0,0);
    errBS(0,1) = beamSpot.covariance()(0,1);
    errBS(1,0) = beamSpot.covariance()(1,0);
    errBS(1,1) = beamSpot.covariance()(1,1);
    
    double LSBS = sqrt(v.Norm2Sqr());
    double LSBSErr = sqrt( v*(errVtx*v) + v*(errBS*v) ) / LSBS;
    
    pair<double,double> LS = make_pair(LSBS,LSBSErr);
    
    return LS;
    
}



double BToKpipiProducer::computeCosAlpha(RefCountedKinematicParticle refitBToKPiPi,
                                         RefCountedKinematicVertex refitVertex,
                                         reco::BeamSpot beamSpot){
    
    TVector v(2);
    v[0] = refitVertex->position().x()-beamSpot.position().x();
    v[1] = refitVertex->position().y()-beamSpot.position().y();
    
    TVector w(2);
    w[0] = refitBToKPiPi->currentState().globalMomentum().x();
    w[1] = refitBToKPiPi->currentState().globalMomentum().y();
    
    double cosAlpha = v*w/sqrt(v.Norm2Sqr()*w.Norm2Sqr());
    return cosAlpha;

}



pair<double,double> BToKpipiProducer::computeDCA(const pat::PackedCandidate &pfCand,
                                                 edm::ESHandle<MagneticField> bFieldHandle,
                                                 reco::BeamSpot beamSpot){
    
    const reco::TransientTrack trackTT((*(pfCand.bestTrack())), &(*bFieldHandle));
    
    TrajectoryStateClosestToPoint theDCAXBS = trackTT.trajectoryStateClosestToPoint( GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()) );
    
    double DCABS = theDCAXBS.perigeeParameters().transverseImpactParameter();
    double DCABSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
    
    pair<double,double> DCA = make_pair(DCABS,DCABSErr);
    
    return DCA;
}






DEFINE_FWK_MODULE(BToKpipiProducer);
