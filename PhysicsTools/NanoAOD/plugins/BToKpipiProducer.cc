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
    double etaMax_;
    double DCASigMin_;
    bool KPiCharge_;
    
    float KaonMass_ = 0.493677;
    float KaonMassErr_ = 1.6e-5;
    float PionMass_ = 0.139570; //TBC
    float PionMassErr_ = 3.5e-7; //TBC

    
};



BToKpipiProducer::BToKpipiProducer(const edm::ParameterSet &iConfig):
beamSpotSrc_( consumes<reco::BeamSpot>                  ( iConfig.getParameter<edm::InputTag>( "beamSpot" ) ) ),
PFCandSrc_(   consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "PFCandCollection" ) ) ),
ptMin_(     iConfig.getParameter<double>( "MinPt" ) ),
etaMax_(    iConfig.getParameter<double>( "MaxEta" ) ),
DCASigMin_( iConfig.getParameter<double>( "MinDCASig") ),
KPiCharge_( iConfig.getParameter<bool>( "KPiChargeCheck" ) )
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
            
            if(kaon_DCABS/kaon_DCABSErr<DCASigMin_) continue;            

            for (unsigned int j = 0; j < pfCandNumber; ++j) {
                
                if(j==i) continue;
                
                const pat::PackedCandidate & piD0 = (*pfCandHandle)[j];

                if(abs(piD0.pdgId())!=211) continue; //Charged hadrons
                if(!piD0.hasTrackDetails()) continue;
                if(piD0.pt()<ptMin_ || abs(piD0.eta())>etaMax_) continue;

                if(KPiCharge_ && kaon.charge()*piD0.charge()>0) continue;
                
                pair<double,double> piD0_DCA = computeDCA(piD0,
                                                          bFieldHandle,
                                                          beamSpot);
                double piD0_DCABS = piD0_DCA.first;
                double piD0_DCABSErr = piD0_DCA.second;
                if(piD0_DCABS/piD0_DCABSErr<DCASigMin_) continue;
                
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
                
                double KPi_mass_err = sqrt(refitKPi->currentState().kinematicParametersError().matrix()(6,6));
                
                math::XYZVector refitKaonV3D_KPi = refitKaon_KPi->refittedTransientTrack().track().momentum();
                math::XYZVector refitPiD0V3D_KPi = refitPiD0_KPi->refittedTransientTrack().track().momentum();
                math::XYZVector refitKPiV3D = refitKaonV3D_KPi + refitPiD0V3D_KPi;                

                for (unsigned int k = 0; k < pfCandNumber; ++k) {
                    
                    if(k==i || k==j) continue;
                    
                    const pat::PackedCandidate & piBu = (*pfCandHandle)[k];
                    if(abs(piBu.pdgId())!=211) continue; //Charged hadrons
                    if(!piBu.hasTrackDetails()) continue;
                    if(piBu.pt()<ptMin_ || abs(piBu.eta())>etaMax_) continue;
                    
                    pair<double,double> piBu_DCA = computeDCA(piBu,
                                                              bFieldHandle,
                                                              beamSpot);
                    double piBu_DCABS = piBu_DCA.first;
                    double piBu_DCABSErr = piBu_DCA.second;
                    
                    if(piBu_DCABS/piBu_DCABSErr<DCASigMin_) continue;
                    
                    RefCountedKinematicVertex refitVertexBToKPiPi;
                    RefCountedKinematicParticle refitBToKPiPi;
                    RefCountedKinematicParticle refitKaon;
                    RefCountedKinematicParticle refitPiD0;
                    RefCountedKinematicParticle refitPiBu;
                    
                    
                    passed = BToKPiPiVertexRefitting(kaon, piD0, piBu,
                                                     bFieldHandle,
                                                     refitVertexBToKPiPi,
                                                     refitBToKPiPi,
                                                     refitKaon,
                                                     refitPiD0,
                                                     refitPiBu);
                    
                    if (!passed) continue;                    

                    pair<double,double> BToKPiPiLS = computeLS(refitVertexBToKPiPi,beamSpot);
                    double LSBS = BToKPiPiLS.first;
                    double LSBSErr = BToKPiPiLS.second;
                    
                    double BToKPiPiVtx_CL = TMath::Prob((double)refitVertexBToKPiPi->chiSquared(),
                                                    int(rint(refitVertexBToKPiPi->degreesOfFreedom())));

                    
                    double cosAlpha = computeCosAlpha(refitBToKPiPi,refitVertexBToKPiPi,beamSpot);
                    
                    double mass_err = sqrt(refitBToKPiPi->currentState().kinematicParametersError().matrix()(6,6));
                    
                    pat::CompositeCandidate BToKPiPiCand;
                    BToKPiPiCand.addDaughter( kaon , "kaon");
                    BToKPiPiCand.addDaughter( piD0 , "piD0");
                    BToKPiPiCand.addDaughter( piBu,  "piBu");
                    
                    math::XYZVector refitKaonV3D = refitKaon->refittedTransientTrack().track().momentum();
                    BToKPiPiCand.addUserFloat("kaon_pt",     sqrt(refitKaonV3D.perp2()));
                    BToKPiPiCand.addUserFloat("kaon_eta",    refitKaonV3D.eta());
                    BToKPiPiCand.addUserFloat("kaon_phi",    refitKaonV3D.phi());
                    BToKPiPiCand.addUserFloat("kaon_charge", refitKaon->currentState().particleCharge());
                    BToKPiPiCand.addUserFloat("kaon_DCASig", kaon_DCABS/kaon_DCABSErr);
                    
                    math::XYZVector refitPiD0V3D = refitPiD0->refittedTransientTrack().track().momentum();
                    BToKPiPiCand.addUserFloat("piD0_pt",     sqrt(refitPiD0V3D.perp2()));
                    BToKPiPiCand.addUserFloat("piD0_eta",    refitPiD0V3D.eta());
                    BToKPiPiCand.addUserFloat("piD0_phi",    refitPiD0V3D.phi());
                    BToKPiPiCand.addUserFloat("piD0_charge", refitPiD0->currentState().particleCharge());
                    BToKPiPiCand.addUserFloat("piD0_DCASig", piD0_DCABS/piD0_DCABSErr);
                
                    math::XYZVector refitPiBuV3D = refitPiBu->refittedTransientTrack().track().momentum();
                    BToKPiPiCand.addUserFloat("piBu_pt",    sqrt(refitPiBuV3D.perp2()));
                    BToKPiPiCand.addUserFloat("piBu_eta",   refitPiBuV3D.eta());
                    BToKPiPiCand.addUserFloat("piBu_phi",   refitPiBuV3D.phi());
                    BToKPiPiCand.addUserFloat("piBu_charge",refitPiBu->currentState().particleCharge());
                    BToKPiPiCand.addUserFloat("piBu_DCASig", piBu_DCABS/piBu_DCABSErr);
                    
                    BToKPiPiCand.addUserFloat("Kpi_pt",     sqrt(refitKPiV3D.perp2()));
                    BToKPiPiCand.addUserFloat("Kpi_eta",    refitKPiV3D.eta());
                    BToKPiPiCand.addUserFloat("Kpi_phi",    refitKPiV3D.phi());
                    BToKPiPiCand.addUserFloat("Kpi_mass",   refitKPi->currentState().mass());
                    BToKPiPiCand.addUserFloat("Kpi_mass_err",   KPi_mass_err);
                    BToKPiPiCand.addUserFloat("Kpi_Lxy", (float) KPiLSBS/KPiLSBSErr);
                    BToKPiPiCand.addUserFloat("Kpi_CL_vtx", (float) KPiVtx_CL);
                    
                    math::XYZVector refitBToKPiPiV3D = refitKaonV3D + refitPiD0V3D + refitPiBuV3D;
                    BToKPiPiCand.addUserFloat("pt",     sqrt(refitBToKPiPiV3D.perp2()));
                    BToKPiPiCand.addUserFloat("eta",    refitBToKPiPiV3D.eta());
                    BToKPiPiCand.addUserFloat("phi",    refitBToKPiPiV3D.phi());
                    BToKPiPiCand.addUserFloat("mass",   refitBToKPiPi->currentState().mass());
                    BToKPiPiCand.addUserFloat("mass_err", mass_err);
                    
                    BToKPiPiCand.addUserFloat("Lxy", (float) LSBS/LSBSErr);
                    BToKPiPiCand.addUserFloat("CL_vtx", (float) BToKPiPiVtx_CL);
                    BToKPiPiCand.addUserFloat("cosAlpha", (float) cosAlpha);
                    
                    result->push_back(BToKPiPiCand);
                    
                    
                }
                
            }
            
        }
        
    }
    
    //iEvent.put(result);
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
