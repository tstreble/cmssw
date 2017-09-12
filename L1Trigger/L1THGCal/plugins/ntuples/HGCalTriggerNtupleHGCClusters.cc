
#include <algorithm>
#include "DataFormats/L1THGCal/interface/HGCalCluster.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerNtupleBase.h"



class HGCalTriggerNtupleHGCClusters : public HGCalTriggerNtupleBase
{

  public:
    HGCalTriggerNtupleHGCClusters(const edm::ParameterSet& conf);
    ~HGCalTriggerNtupleHGCClusters(){};
    virtual void initialize(TTree&, const edm::ParameterSet&, edm::ConsumesCollector&&) override final;
    virtual void fill(const edm::Event& e, const edm::EventSetup& es) override final;

  private:
    virtual void clear() override final;


    edm::EDGetToken clusters_token_;

    int cl_n_ ;
    std::vector<float> cl_pt_;
    std::vector<float> cl_energy_;
    std::vector<float> cl_eta_;
    std::vector<float> cl_phi_;
    std::vector<float> cl_mipPt_;
    std::vector<float> cl_x_;
    std::vector<float> cl_y_;
    std::vector<float> cl_z_;
    std::vector<float> cl_xN_;
    std::vector<float> cl_yN_;
    std::vector<int> cl_layer_;
    std::vector<int> cl_ncells_;   
    std::vector<int> cls_l_;
    std::vector<int> cls_NxLayer_;   
    std::vector<int> cls_zside_;
    std::vector<int> cls_Nseeds_;
    std::vector<int> clsTc_NseedsXlayer_;
    std::vector<int> clsTc_l_;
    std::vector<std::vector<unsigned>> cl_cells_;   

};

DEFINE_EDM_PLUGIN(HGCalTriggerNtupleFactory,
    HGCalTriggerNtupleHGCClusters,
    "HGCalTriggerNtupleHGCClusters" );


HGCalTriggerNtupleHGCClusters::
HGCalTriggerNtupleHGCClusters(const edm::ParameterSet& conf):HGCalTriggerNtupleBase(conf)
{
}

void
HGCalTriggerNtupleHGCClusters::
initialize(TTree& tree, const edm::ParameterSet& conf, edm::ConsumesCollector&& collector)
{
  clusters_token_ = collector.consumes<l1t::HGCalClusterBxCollection>(conf.getParameter<edm::InputTag>("Clusters"));

  tree.Branch("cl_n", &cl_n_, "cl_n/I");
  tree.Branch("cl_pt", &cl_pt_);
  tree.Branch("cl_energy", &cl_energy_);
  tree.Branch("cl_eta", &cl_eta_);
  tree.Branch("cl_phi", &cl_phi_);  
  tree.Branch("cl_layer", &cl_layer_);
  tree.Branch("cl_ncells", &cl_ncells_);
  tree.Branch("cl_cells", &cl_cells_);
  tree.Branch("cl_mipPt", &cl_mipPt_);  
  tree.Branch("cl_x", &cl_x_);  
  tree.Branch("cl_y", &cl_y_);  
  tree.Branch("cl_z", &cl_z_);  
  tree.Branch("cl_xN", &cl_xN_);  
  tree.Branch("cl_yN", &cl_yN_);  
  tree.Branch("cls_l", &cls_l_);
  tree.Branch("cls_NxLayer", &cls_NxLayer_);
  tree.Branch("cls_zside", &cls_zside_);
  tree.Branch("cls_Nseeds", &cls_Nseeds_);
  tree.Branch("clsTc_NseedsXlayer", &clsTc_NseedsXlayer_);
  tree.Branch("clsTc_l", &clsTc_l_);


}

void
HGCalTriggerNtupleHGCClusters::
fill(const edm::Event& e, const edm::EventSetup& es)
{

  // retrieve clusters
  edm::Handle<l1t::HGCalClusterBxCollection> clusters_h;
  e.getByToken(clusters_token_, clusters_h);
  const l1t::HGCalClusterBxCollection& clusters = *clusters_h;

  // retrieve geometry
  edm::ESHandle<HGCalTriggerGeometryBase> geometry;
  es.get<CaloGeometryRecord>().get(geometry);

  clear();
  for(auto cl_itr=clusters.begin(0); cl_itr!=clusters.end(0); cl_itr++)
  {
    cl_n_++;
    // physical values 
    cl_pt_.emplace_back(cl_itr->pt());
    cl_energy_.emplace_back(cl_itr->energy());
    cl_eta_.emplace_back(cl_itr->eta());
    cl_phi_.emplace_back(cl_itr->phi());
    cl_ncells_.emplace_back(cl_itr->constituents().size());
    cl_mipPt_.emplace_back(cl_itr->mipPt());
    cl_x_.emplace_back(cl_itr->centre().x());
    cl_y_.emplace_back(cl_itr->centre().y());
    cl_z_.emplace_back(cl_itr->centre().z());
    cl_xN_.emplace_back(cl_itr->centreProj().x());
    cl_yN_.emplace_back(cl_itr->centreProj().y());



    // Retrieve indices of trigger cells inside cluster
    cl_cells_.emplace_back(cl_itr->constituents().size());
    std::transform(cl_itr->constituents_begin(), cl_itr->constituents_end(),
        cl_cells_.back().begin(), [](const edm::Ptr<l1t::HGCalTriggerCell>& tc){return tc.key();}
        );
    if(cl_itr->subdetId()==3){
        cl_layer_.emplace_back(cl_itr->layer());
    }else if(cl_itr->subdetId()==4){
        cl_layer_.emplace_back(cl_itr->layer()+28);
    }else if(cl_itr->subdetId()==5){
        cl_layer_.emplace_back(cl_itr->layer()+12+28);
    }
    /* loop over the trigger cells pertinent to a cluster */
    const edm::PtrVector<l1t::HGCalTriggerCell> pertinentTC = cl_itr->constituents();
    int Nseeds = 0;
    for( edm::PtrVector<l1t::HGCalTriggerCell>::const_iterator it_tc=pertinentTC.begin(); it_tc<pertinentTC.end(); it_tc++){
        double tc_mipPt = (*it_tc)->mipPt();
        if(tc_mipPt>5){
            Nseeds++;
        }
    }
    cls_Nseeds_.emplace_back(Nseeds);
  }
  
  int zsides[2] = {-1, 1};
  for(int i_z(0); i_z<2; ++i_z){
      for( int l(0); l<=40; l++ ){
          int clNxLayer = 0;
          int layerN = -1;    
          for(auto cl_itr=clusters.begin(0); cl_itr!=clusters.end(0); cl_itr++){             
              if(cl_itr->subdetId()==3 && cl_itr->zside() == zsides[i_z]){
                  layerN=cl_itr->layer();
              }
              else if(cl_itr->subdetId()==4 && cl_itr->zside() == zsides[i_z]){
                  layerN=cl_itr->layer()+28;
              }                    
              if(layerN==l && cl_itr->zside() == zsides[i_z] ){
                  /* loop over the trigger cells pertinent to a cluster */
                  const edm::PtrVector<l1t::HGCalTriggerCell> pertinentTC = cl_itr->constituents();                  
                  int NseedsXcl = 0;                  
                  for( edm::PtrVector<l1t::HGCalTriggerCell>::const_iterator it_tc=pertinentTC.begin(); it_tc<pertinentTC.end(); it_tc++){
                      if( (*it_tc)->mipPt()>5 ){
                          NseedsXcl++;
                      }
                  }
                  clsTc_NseedsXlayer_.emplace_back(NseedsXcl);
                  clsTc_l_.emplace_back(l);
                  clNxLayer++;
              }              
          }
          /* in the event there are clNxlayer */
          cls_l_.emplace_back(l);
          cls_NxLayer_.emplace_back(clNxLayer);
          cls_zside_.emplace_back(zsides[i_z]);
      }      
  }
}


void
HGCalTriggerNtupleHGCClusters::
clear()
{
  cl_n_ = 0;
  cl_pt_.clear();
  cl_energy_.clear();
  cl_eta_.clear();
  cl_phi_.clear();
  cl_layer_.clear();
  cl_ncells_.clear();
  cl_cells_.clear();
  cl_mipPt_.clear();
  cl_x_.clear();
  cl_y_.clear();
  cl_z_.clear();
  cl_xN_.clear();
  cl_yN_.clear();
  cls_l_.clear();
  cls_NxLayer_.clear();
  cls_zside_.clear();
  cls_Nseeds_.clear();
  clsTc_NseedsXlayer_.clear();
  clsTc_l_.clear();

}




