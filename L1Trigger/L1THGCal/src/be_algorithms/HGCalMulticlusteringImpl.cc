

#include "L1Trigger/L1THGCal/interface/be_algorithms/HGCalMulticlusteringImpl.h"
#include "L1Trigger/L1THGCal/interface/be_algorithms/HGCalShowerShape.h"
#include "DataFormats/Math/interface/deltaR.h"


HGCalMulticlusteringImpl::HGCalMulticlusteringImpl( const edm::ParameterSet& conf ) :
    dr_(conf.getParameter<double>("dR_multicluster")),
    ptC3dThreshold_(conf.getParameter<double>("minPt_multicluster")),
    multiclusterAlgoType_(conf.getParameter<string>("type_multicluster")),
    distDbscan_(conf.getParameter<double>("dist_dbscan_multicluster")),
    minNDbscan_(conf.getParameter<unsigned>("minN_dbscan_multicluster")),
    memoryMultiCone_(conf.getParameter<unsigned>("memory_multicluster")),
    dRMultiCone_(conf.getParameter<vector<double>>("dR_MultiCone_multicluster"))
{    
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster dR for Near Neighbour search: " << dr_;  
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster minimum transverse-momentum: " << ptC3dThreshold_;
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster DBSCAN Clustering distance: " << distDbscan_;
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster clustering min number of subclusters: " << minNDbscan_;
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster type of multiclustering algortihm: " << multiclusterAlgoType_;
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster multi cone memory: " << memoryMultiCone_;
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster multi cone dR parameters: {";
    for(auto dR : dRMultiCone_) edm::LogInfo("HGCalMulticlusterParameters") << dR << ",";
    edm::LogInfo("HGCalMulticlusterParameters") << "}";
    id_.reset( HGCalTriggerClusterIdentificationFactory::get()->create("HGCalTriggerClusterIdentificationBDT") );
    id_->initialize(conf.getParameter<edm::ParameterSet>("EGIdentification")); 
}


bool HGCalMulticlusteringImpl::isPertinent( const l1t::HGCalCluster & clu, 
                                            const l1t::HGCalMulticluster & mclu, 
                                            double dR ) const
{
    HGCalDetId cluDetId( clu.detId() );
    HGCalDetId firstClusterDetId( mclu.detId() );
    
    if( cluDetId.zside() != firstClusterDetId.zside() ){
        return false;
    }
    if( ( mclu.centreProj() - clu.centreProj() ).mag() < dR ){
        return true;
    }
    return false;

}


double HGCalMulticlusteringImpl::distanceMultiCone( const l1t::HGCalCluster & clu,
						    const l1t::HGCalMulticluster & mclu,
						    int memory) const
{

  double dRmin = -1.;

  HGCalDetId cluDetId( clu.detId() );
  HGCalDetId firstClusterDetId( mclu.detId() );

  if( cluDetId.zside() != firstClusterDetId.zside() ){
    return dRmin;
  }

  int layer1 = triggerTools_.layerWithOffset(cluDetId);

  const std::unordered_map<uint32_t, edm::Ptr<l1t::HGCalCluster>>& clusters = mclu.constituents();
  for( const auto& id_cluster : clusters ){
    int layer2 = triggerTools_.layerWithOffset(id_cluster.first);
    if( abs(layer2-layer1) > memory ) continue;
    double dR = ( id_cluster.second->centreProj() - clu.centreProj() ).mag();
    if( dRmin < 0 || dR < dRmin ) dRmin = dR;
  }

  return dRmin;

}




void HGCalMulticlusteringImpl::findNeighbor( const std::vector<std::pair<unsigned int,double>>&  rankedList,
                                             unsigned int searchInd,
                                             const std::vector<edm::Ptr<l1t::HGCalCluster>>& clustersPtrs, 
                                             std::vector<unsigned int>& neighbors
                                            ){

    if(clustersPtrs.size() <= searchInd || clustersPtrs.size() < rankedList.size()){
        throw cms::Exception("IndexOutOfBound: clustersPtrs in 'findNeighbor'");
    }

    for(unsigned int ind = searchInd+1; ind < rankedList.size() && fabs(rankedList.at(ind).second - rankedList.at(searchInd).second) < distDbscan_ ; ind++){

        if(clustersPtrs.size() <= rankedList.at(ind).first){
            throw cms::Exception("IndexOutOfBound: clustersPtrs in 'findNeighbor'");

        } else if(((*(clustersPtrs[rankedList.at(ind).first])).centreProj() - (*(clustersPtrs[rankedList.at(searchInd).first])).centreProj()).mag() < distDbscan_){
            neighbors.push_back(ind);
        }
    }

    for(unsigned int ind = 0; ind < searchInd && fabs(rankedList.at(searchInd).second - rankedList.at(ind).second) < distDbscan_ ; ind++){

        if(clustersPtrs.size() <= rankedList.at(ind).first){
            throw cms::Exception("IndexOutOfBound: clustersPtrs in 'findNeighbor'");

        } else if(((*(clustersPtrs[rankedList.at(ind).first])).centreProj() - (*(clustersPtrs[rankedList.at(searchInd).first])).centreProj()).mag() < distDbscan_){
            neighbors.push_back(ind);
        }
    }
}


void HGCalMulticlusteringImpl::clusterizeDR( const std::vector<edm::Ptr<l1t::HGCalCluster>> & clustersPtrs, 
                                           l1t::HGCalMulticlusterBxCollection & multiclusters,
                                           const HGCalTriggerGeometryBase & triggerGeometry)
{

    std::vector<l1t::HGCalMulticluster> multiclustersTmp;

    int iclu = 0;
    for(std::vector<edm::Ptr<l1t::HGCalCluster>>::const_iterator clu = clustersPtrs.begin(); clu != clustersPtrs.end(); ++clu, ++iclu){
        
        int imclu=0;
        vector<int> tcPertinentMulticlusters;
        for( const auto& mclu : multiclustersTmp ){
            if( this->isPertinent(**clu, mclu, dr_) ){
                tcPertinentMulticlusters.push_back(imclu);
            }
            ++imclu;
        }
        if( tcPertinentMulticlusters.empty() ){
            multiclustersTmp.emplace_back( *clu );
        }
        else{
            unsigned minDist = 1;
            unsigned targetMulticlu = 0; 
            for( int imclu : tcPertinentMulticlusters ){
                double d = ( multiclustersTmp.at(imclu).centreProj() - (*clu)->centreProj() ).mag() ;
                if( d < minDist ){
                    minDist = d;
                    targetMulticlu = imclu;
                }
            } 

            multiclustersTmp.at( targetMulticlu ).addConstituent( *clu );
            
        }        
    }

    /* making the collection of multiclusters */
    finalizeClusters(multiclustersTmp, multiclusters, triggerGeometry);
    
}





void HGCalMulticlusteringImpl::clusterizeMultiCone( const std::vector<edm::Ptr<l1t::HGCalCluster>> & clustersPtrs,
						    l1t::HGCalMulticlusterBxCollection & multiclusters,
						    const HGCalTriggerGeometryBase & triggerGeometry)
{

  if(multiclusterAlgoType_=="MultiConeC3d" && dRMultiCone_.size()!=triggerTools_.lastLayerBH()+1){
    throw edm::Exception(edm::errors::Configuration, "Configuration")
      << "HGCalMulticlusteringImpl wrong number of dR_MultiCone_multicluster elements: "<<triggerTools_.lastLayerBH()+1<<" needed, "<<dRMultiCone_.size()<<" provided"<<endl;
  }

  std::unordered_map<int,std::vector<edm::Ptr<l1t::HGCalCluster> > > C2Ds;

  for(std::vector<edm::Ptr<l1t::HGCalCluster>>::const_iterator clu = clustersPtrs.begin(); clu != clustersPtrs.end(); ++clu){

    unsigned layer = triggerTools_.layerWithOffset((**clu).detId());
    C2Ds[layer].emplace_back(*clu);

  }


  std::vector<l1t::HGCalMulticluster> multiclustersTmp;

  for(unsigned layer=1;layer<=triggerTools_.lastLayerBH();layer++){

    double dR = dRMultiCone_[layer];

    for(auto clu : C2Ds[layer]){

      if(dR==0.){
	throw cms::Exception("BadConfiguration")
	  <<"3D cluster multicone dR forced to 0 by coefficients.\n"
	  <<"The configuration should be changed. "
	  <<"Discarded layers should be defined in hgcalTriggerGeometryESProducer.TriggerGeometry.DisconnectedLayers and not with multicone dR coefficients = 0\n";
      }      

      int imclu=0;
      vector<pair<double,int> > tcPertinentMulticlusters;
      for( const auto& mclu : multiclustersTmp ){
	double dist = this->distanceMultiCone(*clu, mclu, memoryMultiCone_);
	if( dist>0 && dist<dR ){
	  tcPertinentMulticlusters.push_back(make_pair(dist,imclu));
	}
	++imclu;
      }

      if( tcPertinentMulticlusters.empty() ){
	multiclustersTmp.emplace_back( clu );
      }
      else{
	unsigned minDist = 1;
	unsigned targetMulticlu = 0;
	for( auto dist_imclu : tcPertinentMulticlusters ){
	  double d = dist_imclu.first ;
	  if( d < minDist ){
	    minDist = d;
	    targetMulticlu = dist_imclu.second;
	  }
	}

	multiclustersTmp.at( targetMulticlu ).addConstituent( clu );

      }

    }

  }

  /* making the collection of multiclusters */
  finalizeClusters(multiclustersTmp, multiclusters, triggerGeometry);

}




void HGCalMulticlusteringImpl::clusterizeDBSCAN( const std::vector<edm::Ptr<l1t::HGCalCluster>> & clustersPtrs, 
                                                 l1t::HGCalMulticlusterBxCollection & multiclusters,
                                                 const HGCalTriggerGeometryBase & triggerGeometry)
{

    std::vector<l1t::HGCalMulticluster> multiclustersTmp;
    l1t::HGCalMulticluster mcluTmp;
    std::vector<bool> visited(clustersPtrs.size(),false);
    std::vector<bool> merged (clustersPtrs.size(),false);
    std::vector<std::pair<unsigned int,double>>  rankedList;
    rankedList.reserve(clustersPtrs.size());
    std::vector<std::vector<unsigned int>> neighborList;
    neighborList.reserve(clustersPtrs.size());

    int iclu = 0, imclu = 0, neighNo;
    double dist = 0.;

    for(std::vector<edm::Ptr<l1t::HGCalCluster>>::const_iterator clu = clustersPtrs.begin(); clu != clustersPtrs.end(); ++clu, ++iclu){
        dist = (*clu)->centreProj().mag()*HGCalDetId((*clu)->detId()).zside();
        rankedList.push_back(std::make_pair(iclu,dist));
    }  
    iclu = 0;
    std::sort(rankedList.begin(), rankedList.end(), [](auto &left, auto &right) {
            return left.second < right.second;
            });

    for(const auto& cluRanked: rankedList){
        std::vector<unsigned int> neighbors;      

        if(!visited.at(iclu)){
            visited.at(iclu) = true;
            findNeighbor(rankedList, iclu, clustersPtrs, neighbors);
            neighborList.push_back(std::move(neighbors));

            if(neighborList.at(iclu).size() >= minNDbscan_) {
                multiclustersTmp.emplace_back( clustersPtrs[cluRanked.first] );
                merged.at(iclu) = true;
                /* dynamic range loop: range-based loop syntax cannot be employed */
                for(unsigned int neighInd = 0; neighInd < neighborList.at(iclu).size(); neighInd++){
                    neighNo = neighborList.at(iclu).at(neighInd);
                    /* This condition also ensures merging of clusters visited by other clusters but not merged. */
                    if(!merged.at(neighNo) ){
                        merged.at(neighNo) = true;          
                        multiclustersTmp.at(imclu).addConstituent( clustersPtrs[rankedList.at(neighNo).first] );

                        if(!visited.at(neighNo)){
                            visited.at(neighNo) = true;
                            std::vector<unsigned int> secNeighbors;
                            findNeighbor(rankedList, neighNo,clustersPtrs, secNeighbors);

                            if(secNeighbors.size() >= minNDbscan_){
                                neighborList.at(iclu).insert(neighborList.at(iclu).end(), secNeighbors.begin(), secNeighbors.end());
                            }
                        }
                    }
                }
                imclu++;
            }
        }
        else neighborList.push_back(std::move(neighbors));
        iclu++;    
    }
    /* making the collection of multiclusters */
    finalizeClusters(multiclustersTmp, multiclusters, triggerGeometry);
}


void
HGCalMulticlusteringImpl::
finalizeClusters(std::vector<l1t::HGCalMulticluster>& multiclusters_in,
            l1t::HGCalMulticlusterBxCollection& multiclusters_out, 
            const HGCalTriggerGeometryBase& triggerGeometry) {
    for(auto& multicluster : multiclusters_in) {
        // compute the eta, phi from its barycenter
        // + pT as scalar sum of pT of constituents
        double sumPt=0.;
        const std::unordered_map<uint32_t, edm::Ptr<l1t::HGCalCluster>>& clusters = multicluster.constituents();
        for(const auto& id_cluster : clusters) sumPt += id_cluster.second->pt();

        math::PtEtaPhiMLorentzVector multiclusterP4(  sumPt,
                multicluster.centre().eta(),
                multicluster.centre().phi(),
                0. );
        multicluster.setP4( multiclusterP4 );

        if( multicluster.pt() > ptC3dThreshold_ ){
            //compute shower shapes
            multicluster.showerLength(shape_.showerLength(multicluster));
            multicluster.coreShowerLength(shape_.coreShowerLength(multicluster, triggerGeometry));
            multicluster.firstLayer(shape_.firstLayer(multicluster));
            multicluster.maxLayer(shape_.maxLayer(multicluster));
            multicluster.sigmaEtaEtaTot(shape_.sigmaEtaEtaTot(multicluster));
            multicluster.sigmaEtaEtaMax(shape_.sigmaEtaEtaMax(multicluster));
            multicluster.sigmaPhiPhiTot(shape_.sigmaPhiPhiTot(multicluster));
            multicluster.sigmaPhiPhiMax(shape_.sigmaPhiPhiMax(multicluster));
            multicluster.sigmaZZ(shape_.sigmaZZ(multicluster));
            multicluster.sigmaRRTot(shape_.sigmaRRTot(multicluster));
            multicluster.sigmaRRMax(shape_.sigmaRRMax(multicluster));
            multicluster.sigmaRRMean(shape_.sigmaRRMean(multicluster));
            multicluster.eMax(shape_.eMax(multicluster));
            // fill quality flag
            multicluster.setHwQual(id_->decision(multicluster));

            multiclusters_out.push_back( 0, multicluster);
        }
    }
}
