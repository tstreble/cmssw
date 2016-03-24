///
/// \class l1t::GtExternalFakeProducer
///
/// Description: Fill uGT external condition to allow testing stage 2 algos, e.g. Bptx
///
/// 
/// \author: D. Puigh OSU
///


// system include files
#include <boost/shared_ptr.hpp>

// user include files

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//#include <vector>
#include "DataFormats/L1Trigger/interface/BXVector.h"

#include "DataFormats/L1TGlobal/interface/GlobalExtBlk.h"

using namespace std;
using namespace edm;


namespace l1t {

  //
  // class declaration
  //

  class GtExternalFakeProducer : public EDProducer {
  public:
    explicit GtExternalFakeProducer(const ParameterSet&);
    ~GtExternalFakeProducer();

    static void fillDescriptions(ConfigurationDescriptions& descriptions);

  private:
    virtual void produce(Event&, EventSetup const&);
    virtual void beginJob();
    virtual void endJob();
    virtual void beginRun(Run const&iR, EventSetup const&iE);
    virtual void endRun(Run const& iR, EventSetup const& iE);

    // ----------member data ---------------------------
    unsigned long long m_paramsCacheId; // Cache-ID from current parameters, to check if needs to be updated.
    //boost::shared_ptr<const CaloParams> m_dbpars; // Database parameters for the trigger, to be updated as needed.
    //boost::shared_ptr<const FirmwareVersion> m_fwv;
    //boost::shared_ptr<FirmwareVersion> m_fwv; //not const during testing.

    // BX parameters
    int bxFirst_;
    int bxLast_;

    bool setBptxAND_;
    bool setBptxPlus_;
    bool setBptxMinus_;
    bool setBptxOR_;

  };

  //
  // constructors and destructor
  //
  GtExternalFakeProducer::GtExternalFakeProducer(const ParameterSet& iConfig)
  {
    // register what you produce
    produces<GlobalExtBlkBxCollection>();

    // Setup parameters
    bxFirst_ = iConfig.getParameter<int>("bxFirst");
    bxLast_  = iConfig.getParameter<int>("bxLast");

    setBptxAND_ = iConfig.getParameter<bool>("setBptxAND");
    setBptxPlus_ = iConfig.getParameter<bool>("setBptxPlus");
    setBptxMinus_ = iConfig.getParameter<bool>("setBptxMinus");
    setBptxOR_ = iConfig.getParameter<bool>("setBptxOR");

    // set cache id to zero, will be set at first beginRun:
    m_paramsCacheId = 0;
  }


  GtExternalFakeProducer::~GtExternalFakeProducer()
  {
  }



  //
  // member functions
  //

  // ------------ method called to produce the data ------------
  void
  GtExternalFakeProducer::produce(Event& iEvent, const EventSetup& iSetup)
  {

    LogDebug("GtExternalFakeProducer") << "GtExternalFakeProducer::produce function called...\n";

    // Setup vectors
    GlobalExtBlk extCond_bx;

    // Set the range of BX....TO DO...move to Params or determine from param set.
    int bxFirst = bxFirst_;
    int bxLast  = bxLast_;


    //outputs
    std::auto_ptr<GlobalExtBlkBxCollection> extCond( new GlobalExtBlkBxCollection(0,bxFirst,bxLast));


    extCond_bx.reset();
    // Fill in some external conditions for testing
    if( setBptxAND_ ) extCond_bx.setExternalDecision(8,true);  //EXT_BPTX_plus_AND_minus.v0
    if( setBptxPlus_ ) extCond_bx.setExternalDecision(9,true);  //EXT_BPTX_plus.v0
    if( setBptxMinus_ ) extCond_bx.setExternalDecision(10,true); //EXT_BPTX_minus.v0
    if( setBptxOR_ ) extCond_bx.setExternalDecision(11,true); //EXT_BPTX_plus_OR_minus.v0

    // Fill Externals
    extCond->push_back(-2, extCond_bx);
    extCond->push_back(-1, extCond_bx);
    extCond->push_back(0,  extCond_bx);
    extCond->push_back(1,  extCond_bx);
    extCond->push_back(2,  extCond_bx); 
   

    iEvent.put(extCond);

  }

  // ------------ method called once each job just before starting event loop ------------
  void
  GtExternalFakeProducer::beginJob()
  {
  }

  // ------------ method called once each job just after ending the event loop ------------
  void
  GtExternalFakeProducer::endJob() {
  }

  // ------------ method called when starting to processes a run ------------

  void GtExternalFakeProducer::beginRun(Run const&iR, EventSetup const&iE){

    LogDebug("GtExternalFakeProducer") << "GtExternalFakeProducer::beginRun function called...\n";

  }

  // ------------ method called when ending the processing of a run ------------
  void GtExternalFakeProducer::endRun(Run const& iR, EventSetup const& iE){

  }

  // ------------ method fills 'descriptions' with the allowed parameters for the module ------------
  void
  GtExternalFakeProducer::fillDescriptions(ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
  }

} // namespace

//define this as a plug-in
DEFINE_FWK_MODULE(l1t::GtExternalFakeProducer);
