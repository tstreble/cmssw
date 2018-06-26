#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <TLorentzVector.h>
#include <TVector.h>
#include <TMatrix.h>


//
// class declaration
//

using namespace std;

class BToKeTrkProducer : public edm::EDProducer {
    
public:
    
    explicit BToKeTrkProducer(const edm::ParameterSet &iConfig);
    
    ~BToKeTrkProducer() override {};
    
    
private:
    
    virtual void produce(edm::Event&, const edm::EventSetup&);
    
    bool ETrkVertexRefitting(const pat::Electron & ele,
			     const reco::Track & trk,
			     edm::ESHandle<TransientTrackBuilder> theTTBuilder,
			     RefCountedKinematicVertex &refitVertex,
			     RefCountedKinematicParticle &refitETrk,
			     RefCountedKinematicParticle &refitEle,
			     RefCountedKinematicParticle &refitTrk);
    
    bool BToKETrkVertexRefitting(const pat::Electron &ele,
				 const reco::Track &trk,
				 const pat::PackedCandidate &kaon,
				 edm::ESHandle<TransientTrackBuilder> theTTBuilder,
				 RefCountedKinematicVertex &refitVertex,
				 RefCountedKinematicParticle &refitBToKETrk,
				 RefCountedKinematicParticle &refitEle,
				 RefCountedKinematicParticle &refitTrk,
				 RefCountedKinematicParticle &refitKaon);
    
    pair<double,double> computeLS(RefCountedKinematicVertex refitVertex,
				  reco::BeamSpot beamSpot);
    
    double computeCosAlpha(RefCountedKinematicParticle refitBToKETrk,
                           RefCountedKinematicVertex vertexFitTree,
                           reco::BeamSpot beamSpot);

    pair<double,double> computeDCA(const pat::PackedCandidate &kaon,
				   edm::ESHandle<MagneticField> bFieldHandle,
				   reco::BeamSpot beamSpot);
    
    // ----------member data ---------------------------
    
    edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
    edm::EDGetTokenT<std::vector<pat::Electron>> electronSrc_;
    edm::EDGetTokenT< std::vector<reco::Track> > generalTrackSrc_;
    edm::EDGetTokenT<edm::View<pat::PackedCandidate>> PFCandSrc_;
    
    double ptMinEle_;
    double etaMaxEle_;
    double ptMinTrk_;
    double etaMaxTrk_;
    double ptMinKaon_;
    double etaMaxKaon_;
    double DCASigMinKaon_;
    bool diEleCharge_;
    
    float ElectronMass_ = 0.5109989e-3;
    float ElectronMassErr_ = 3.1*1e-12;
    float KaonMass_ = 0.493677;
    float KaonMassErr_ = 1.6e-5;
    
};



BToKeTrkProducer::BToKeTrkProducer(const edm::ParameterSet &iConfig):
beamSpotSrc_( consumes<reco::BeamSpot> ( iConfig.getParameter<edm::InputTag>( "beamSpot" ) ) ),
electronSrc_( consumes<std::vector<pat::Electron>> ( iConfig.getParameter<edm::InputTag>( "electronCollection" ) ) ),
generalTrackSrc_( consumes<std::vector<reco::Track>> ( iConfig.getParameter<edm::InputTag>( "generalTrackCollection" ) ) ),
PFCandSrc_( consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "PFCandCollection" ) ) ),
ptMinEle_( iConfig.getParameter<double>( "ElectronMinPt" ) ),
etaMaxEle_( iConfig.getParameter<double>( "ElectronMaxEta" ) ),
ptMinTrk_( iConfig.getParameter<double>( "TrackMinPt" ) ),
etaMaxTrk_( iConfig.getParameter<double>( "TrackMaxEta" ) ),
ptMinKaon_( iConfig.getParameter<double>( "KaonMinPt" ) ),
etaMaxKaon_( iConfig.getParameter<double>( "KaonMaxEta" ) ),
DCASigMinKaon_( iConfig.getParameter<double>( "KaonMinDCASig" ) ),
diEleCharge_( iConfig.getParameter<bool>( "DiElectronChargeCheck" ) )
{
    produces<pat::CompositeCandidateCollection>();
}


void BToKeTrkProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
  edm::ESHandle<MagneticField> bFieldHandle;
  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  edm::Handle<reco::VertexCollection> vertexHandle;   
  
  iEvent.getByToken(beamSpotSrc_, beamSpotHandle);
  
  if ( ! beamSpotHandle.isValid() ) {
    edm::LogError("BToKeTrkProducer") << "No beam spot available from EventSetup" ;
  }
  
  reco::BeamSpot beamSpot = *beamSpotHandle;
  
  edm::Handle<std::vector<pat::Electron>> electronHandle;
  edm::Handle<edm::View<pat::PackedCandidate>> pfCandHandle;
  edm::Handle<std::vector<reco::Track> > generalTrackHandle;

  iEvent.getByToken(electronSrc_, electronHandle);
  iEvent.getByToken(PFCandSrc_, pfCandHandle);
  try { iEvent.getByToken(generalTrackSrc_, generalTrackHandle); } 
  catch (...) { std::cout << "Problem with generalTrackHandle" << std::endl; }    
  
  unsigned int electronNumber = electronHandle->size();    
  unsigned int pfCandNumber = pfCandHandle->size();
  unsigned int trackNumber = generalTrackHandle->size();
  
  // Output collection
  std::unique_ptr<pat::CompositeCandidateCollection> result( new pat::CompositeCandidateCollection );

  
  // loop on all the eeK triplets
  for (unsigned int i = 0; i < electronNumber; ++i) {
    
    const pat::Electron & ele = (*electronHandle)[i];
    
    //could implement ele ID criteria here
    if(ele.pt()<ptMinEle_ || abs(ele.eta())>etaMaxEle_) continue;
    
    for (unsigned int j = 0; j < trackNumber; ++j) {
      
      const reco::Track & trk = (*generalTrackHandle)[j];
      
      if(trk.pt()<ptMinTrk_ || abs(trk.eta())>etaMaxTrk_) continue;     
      if(diEleCharge_ && ele.charge()*trk.charge()>0) continue;
      if(deltaR(ele,trk)<0.1) continue; //rough overlap criteria, to be refined?
      
      RefCountedKinematicVertex refitVertexETrk;
      RefCountedKinematicParticle refitETrk;
      RefCountedKinematicParticle refitEle_ETrk;
      RefCountedKinematicParticle refitTrk_ETrk;
      
      bool passed = ETrkVertexRefitting(ele, trk,
					theTTBuilder,
					refitVertexETrk,
					refitETrk,
					refitEle_ETrk,
					refitTrk_ETrk);
      
      if ( !passed) continue;
      
      pair<double,double> ETrkLS = computeLS(refitVertexETrk,beamSpot);
      double ETrkLSBS = ETrkLS.first;
      double ETrkLSBSErr = ETrkLS.second;
      
      
      double ETrkVtx_CL = TMath::Prob((double)refitVertexETrk->chiSquared(),
				      int(rint(refitVertexETrk->degreesOfFreedom())));
      
      double ETrk_mass_err = sqrt(refitETrk->currentState().kinematicParametersError().matrix()(6,6));
      
      math::XYZVector refitEleV3D_ETrk = refitEle_ETrk->refittedTransientTrack().track().momentum();
      math::XYZVector refitTrkV3D_ETrk = refitTrk_ETrk->refittedTransientTrack().track().momentum();		
      math::XYZVector refitETrkV3D = refitEleV3D_ETrk + refitTrkV3D_ETrk;
      
      //Kaon
      for (unsigned int k = 0; k < pfCandNumber; ++k) {
	
	const pat::PackedCandidate & pfCand = (*pfCandHandle)[k];

	if(abs(pfCand.pdgId())!=211) continue; //Charged hadrons
	if(!pfCand.hasTrackDetails()) continue;
	if(pfCand.pt()<ptMinKaon_ || abs(pfCand.eta())>etaMaxKaon_) continue;
	if(deltaR(pfCand,trk)<0.1) continue; //rough overlap criteria, to be refined?
	
	pair<double,double> DCA = computeDCA(pfCand,
					     bFieldHandle,
					     beamSpot);
	double DCABS = DCA.first;
	double DCABSErr = DCA.second;
	
	if(fabs(DCABS/DCABSErr)<DCASigMinKaon_) continue;
          
	RefCountedKinematicVertex refitVertexBToKETrk;
	RefCountedKinematicParticle refitBToKETrk;
	RefCountedKinematicParticle refitEle;
	RefCountedKinematicParticle refitTrk;
	RefCountedKinematicParticle refitKaon;
	
        
	passed = BToKETrkVertexRefitting(ele, trk, pfCand,
					 theTTBuilder,
					 refitVertexBToKETrk,
					 refitBToKETrk,
					 refitEle,
					 refitTrk,
					 refitKaon);
	
	if (!passed) continue;
	  
	pair<double,double> BToKETrkLS = computeLS(refitVertexBToKETrk,beamSpot);
	double LSBS = BToKETrkLS.first;
	double LSBSErr = BToKETrkLS.second;
        
	double BToKETrkVtx_CL = TMath::Prob((double)refitVertexBToKETrk->chiSquared(),
					    int(rint(refitVertexBToKETrk->degreesOfFreedom())));

	
	double cosAlpha = computeCosAlpha(refitBToKETrk,refitVertexBToKETrk,beamSpot);
	
	double mass_err = sqrt(refitBToKETrk->currentState().kinematicParametersError().matrix()(6,6));
	
	pat::CompositeCandidate BToKETrkCand;
	BToKETrkCand.addDaughter( ele , "ele");
	BToKETrkCand.addDaughter( pfCand, "kaon");
	//BToKETrkCand.addDaughter( trk , "track"); //No Track to Candidate conversion method available
	
	math::XYZVector refitEleV3D = refitEle->refittedTransientTrack().track().momentum();
	BToKETrkCand.addUserFloat("ele_pt",     sqrt(refitEleV3D.perp2()));
	BToKETrkCand.addUserFloat("ele_eta",    refitEleV3D.eta());
	BToKETrkCand.addUserFloat("ele_phi",    refitEleV3D.phi());
	BToKETrkCand.addUserFloat("ele_charge", refitEle->currentState().particleCharge());
	
	math::XYZVector refitTrkV3D = refitTrk->refittedTransientTrack().track().momentum();
	BToKETrkCand.addUserFloat("trk_pt",     sqrt(refitTrkV3D.perp2()));
	BToKETrkCand.addUserFloat("trk_eta",    refitTrkV3D.eta());
	BToKETrkCand.addUserFloat("trk_phi",    refitTrkV3D.phi());
	BToKETrkCand.addUserFloat("trk_charge", refitTrk->currentState().particleCharge());		    
	
	math::XYZVector refitKaonV3D = refitKaon->refittedTransientTrack().track().momentum();
	BToKETrkCand.addUserFloat("kaon_pt",    sqrt(refitKaonV3D.perp2()));
	BToKETrkCand.addUserFloat("kaon_eta",   refitKaonV3D.eta());
	BToKETrkCand.addUserFloat("kaon_phi",   refitKaonV3D.phi());
	BToKETrkCand.addUserFloat("kaon_charge",refitKaon->currentState().particleCharge());
	BToKETrkCand.addUserFloat("kaon_DCASig", DCABS/DCABSErr);
	
	BToKETrkCand.addUserFloat("ee_pt",     sqrt(refitETrkV3D.perp2()));
	BToKETrkCand.addUserFloat("ee_eta",    refitETrkV3D.eta());
	BToKETrkCand.addUserFloat("ee_phi",    refitETrkV3D.phi());
	BToKETrkCand.addUserFloat("ee_mass",   refitETrk->currentState().mass());
	BToKETrkCand.addUserFloat("ee_mass_err",   ETrk_mass_err);
	
	BToKETrkCand.addUserFloat("ee_Lxy", (float) ETrkLSBS/ETrkLSBSErr);
	BToKETrkCand.addUserFloat("ee_CL_vtx", (float) ETrkVtx_CL);
          
	math::XYZVector refitBToKETrkV3D = refitEleV3D + refitTrkV3D + refitKaonV3D;
	BToKETrkCand.addUserFloat("pt",     sqrt(refitBToKETrkV3D.perp2()));
	BToKETrkCand.addUserFloat("eta",    refitBToKETrkV3D.eta());
	BToKETrkCand.addUserFloat("phi",    refitBToKETrkV3D.phi());
	BToKETrkCand.addUserFloat("mass",   refitBToKETrk->currentState().mass());
	BToKETrkCand.addUserFloat("mass_err", mass_err);
	  
	BToKETrkCand.addUserFloat("Lxy", (float) LSBS/LSBSErr);
	BToKETrkCand.addUserFloat("CL_vtx", (float) BToKETrkVtx_CL);
	BToKETrkCand.addUserFloat("cosAlpha", (float) cosAlpha);
        
	result->push_back(BToKETrkCand);
        
      
      }
	
    }
      
  }
    
  
  iEvent.put(std::move(result));
  
}



bool BToKeTrkProducer::ETrkVertexRefitting(const pat::Electron & ele,
					   const reco::Track & trk,
					   edm::ESHandle<TransientTrackBuilder> theTTBuilder,
					   RefCountedKinematicVertex &refitVertex,
					   RefCountedKinematicParticle &refitETrk,
					   RefCountedKinematicParticle &refitEle,
					   RefCountedKinematicParticle &refitTrk){
    
    const reco::TransientTrack eleTT = theTTBuilder->build(ele.gsfTrack());
    const reco::TransientTrack trkTT = theTTBuilder->build(trk);
    
    KinematicParticleFactoryFromTransientTrack partFactory;
    KinematicParticleVertexFitter PartVtxFitter;
    
    std::vector<RefCountedKinematicParticle> eleParticles;
    double chi = 0.;
    double ndf = 0.;
    eleParticles.push_back(partFactory.particle(eleTT,ElectronMass_,chi,ndf,ElectronMassErr_));
    eleParticles.push_back(partFactory.particle(trkTT,ElectronMass_,chi,ndf,ElectronMassErr_));
    RefCountedKinematicTree eeVertexFitTree = PartVtxFitter.fit(eleParticles);
    
    if ( !eeVertexFitTree->isValid()) return false;
    
    eeVertexFitTree->movePointerToTheTop();
    refitVertex = eeVertexFitTree->currentDecayVertex();
    refitETrk = eeVertexFitTree->currentParticle();
    
    if ( !refitVertex->vertexIsValid()) return false;

    // extract the re-fitted tracks
    eeVertexFitTree->movePointerToTheTop();
    
    eeVertexFitTree->movePointerToTheFirstChild();
    refitEle = eeVertexFitTree->currentParticle();
    
    eeVertexFitTree->movePointerToTheNextChild();
    refitTrk = eeVertexFitTree->currentParticle();
    
    return true;
    
}




bool BToKeTrkProducer::BToKETrkVertexRefitting(const pat::Electron &ele,
					       const reco::Track &trk,
					       const pat::PackedCandidate &kaon,
					       edm::ESHandle<TransientTrackBuilder> theTTBuilder,					   
					       RefCountedKinematicVertex &refitVertex,
					       RefCountedKinematicParticle &refitBToKETrk,
					       RefCountedKinematicParticle &refitEle,
					       RefCountedKinematicParticle &refitTrk,
					       RefCountedKinematicParticle &refitKaon){

  const reco::TransientTrack eleTT = theTTBuilder->build(ele.gsfTrack());
  const reco::TransientTrack trkTT = theTTBuilder->build(trk);
  const reco::TransientTrack kaonTT = theTTBuilder->build(kaon.bestTrack());
  
  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter PartVtxFitter;
    
  std::vector<RefCountedKinematicParticle> BToKETrkParticles;
  double chi = 0.;
  double ndf = 0.;
  BToKETrkParticles.push_back(partFactory.particle(eleTT,ElectronMass_,chi,ndf,ElectronMassErr_));
  BToKETrkParticles.push_back(partFactory.particle(trkTT,ElectronMass_,chi,ndf,ElectronMassErr_));
  BToKETrkParticles.push_back(partFactory.particle(kaonTT,KaonMass_,chi,ndf,KaonMassErr_));

  RefCountedKinematicTree BToKETrkVertexFitTree = PartVtxFitter.fit(BToKETrkParticles);
    
  if ( !BToKETrkVertexFitTree->isValid()) return false;
    
  BToKETrkVertexFitTree->movePointerToTheTop();
  refitVertex = BToKETrkVertexFitTree->currentDecayVertex();
  refitBToKETrk = BToKETrkVertexFitTree->currentParticle();
    
  if ( !refitVertex->vertexIsValid()) return false;
    
  // extract the re-fitted tracks
  BToKETrkVertexFitTree->movePointerToTheTop();
  
  BToKETrkVertexFitTree->movePointerToTheFirstChild();
  refitEle = BToKETrkVertexFitTree->currentParticle();
  
  BToKETrkVertexFitTree->movePointerToTheNextChild();
  refitTrk = BToKETrkVertexFitTree->currentParticle();
  
  BToKETrkVertexFitTree->movePointerToTheNextChild();
  refitKaon = BToKETrkVertexFitTree->currentParticle();
  
  return true;
    
    
}






pair<double,double> BToKeTrkProducer::computeLS(RefCountedKinematicVertex refitVertex,
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



double BToKeTrkProducer::computeCosAlpha(RefCountedKinematicParticle refitBToKETrk,
				       RefCountedKinematicVertex refitVertex,
				       reco::BeamSpot beamSpot){
    
  TVector v(2);
  v[0] = refitVertex->position().x()-beamSpot.position().x();
  v[1] = refitVertex->position().y()-beamSpot.position().y();

  TVector w(2);
  w[0] = refitBToKETrk->currentState().globalMomentum().x();
  w[1] = refitBToKETrk->currentState().globalMomentum().y();
    
  double cosAlpha = v*w/sqrt(v.Norm2Sqr()*w.Norm2Sqr());
  return cosAlpha;
}






pair<double,double> BToKeTrkProducer::computeDCA(const pat::PackedCandidate &kaon,
						 edm::ESHandle<MagneticField> bFieldHandle,
						 reco::BeamSpot beamSpot){
  
  const reco::TransientTrack trackTT((*(kaon.bestTrack())), &(*bFieldHandle));
  
  TrajectoryStateClosestToPoint theDCAXBS = trackTT.trajectoryStateClosestToPoint( GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()) );  
  
  double DCABS = theDCAXBS.perigeeParameters().transverseImpactParameter();
  double DCABSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
    
  pair<double,double> DCA = make_pair(DCABS,DCABSErr);
    
  return DCA;
}




DEFINE_FWK_MODULE(BToKeTrkProducer);
