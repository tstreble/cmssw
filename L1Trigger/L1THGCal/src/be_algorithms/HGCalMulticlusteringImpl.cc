

#include "L1Trigger/L1THGCal/interface/be_algorithms/HGCalMulticlusteringImpl.h"
#include "L1Trigger/L1THGCal/interface/be_algorithms/HGCalShowerShape.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <unordered_map>

HGCalMulticlusteringImpl::HGCalMulticlusteringImpl( const edm::ParameterSet& conf ) :
    dr_(conf.getParameter<double>("dR_multicluster")),
    ptC3dThreshold_(conf.getParameter<double>("minPt_multicluster")),
    multiclusterAlgoType_(conf.getParameter<string>("type_multicluster")),
    distDbscan_(conf.getParameter<double>("dist_dbscan_multicluster")),
    minNDbscan_(conf.getParameter<unsigned>("minN_dbscan_multicluster")),
    nBinsHisto_(conf.getParameter<unsigned>("nBins_histo_multicluster")),
    superClustering_(conf.getParameter<bool>("superClustering"))
{    
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster dR for Near Neighbour search: " << dr_;  
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster minimum transverse-momentum: " << ptC3dThreshold_;
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster DBSCAN Clustering distance: " << distDbscan_;
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster clustering min number of subclusters: " << minNDbscan_;
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster number of bins for the histo algorithm: " << nBinsHisto_<<endl;
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster number of bins for the histo algorithm: " << superClustering_<<endl;
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster type of multiclustering algortihm: " << multiclusterAlgoType_;
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




bool HGCalMulticlusteringImpl::isPertinent( const l1t::HGCalMulticluster & clu,
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




bool HGCalMulticlusteringImpl::isPertinent( const l1t::HGCalCluster & clu,
                                            const GlobalPoint & seed,
                                            double dR ) const
{
    HGCalDetId cluDetId( clu.detId() );

    if( cluDetId.zside()*seed.z()<0) return false;
    else return (seed - clu.centreProj() ).mag() < dR;

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








void HGCalMulticlusteringImpl::clusterizeHistoDR( const std::vector<edm::Ptr<l1t::HGCalCluster>> & clustersPtrs,
						  l1t::HGCalMulticlusterBxCollection & multiclusters,
						  const HGCalTriggerGeometryBase & triggerGeometry)
{


  std::map<std::vector<int>,float> histoMultiCluster; //key[0] = z.side(), key[1] = bin_x, key[2] = bin_y

  for(std::vector<edm::Ptr<l1t::HGCalCluster>>::const_iterator clu = clustersPtrs.begin(); clu != clustersPtrs.end(); ++clu){

    int bin_x = int( (**clu).centreProj().x()*nBinsHisto_/kXYOverZMax_ );
    int bin_y = int( (**clu).centreProj().y()*nBinsHisto_/kXYOverZMax_ );
    std::vector<int> key = { (**clu).zside(), bin_x, bin_y };
    histoMultiCluster[key]+=(**clu).mipPt();

  }


  std::vector<GlobalPoint> seedPositions;
  std::map<int,l1t::HGCalMulticluster> mapSeedMulticluster;
  std::vector<l1t::HGCalMulticluster> multiclustersTmp;

  for(int z_side = -1; z_side<2; z_side++){
    if(z_side==0) continue;

    //Define seeds
    for(int bin_x = -nBinsHisto_+1; bin_x<int(nBinsHisto_); bin_x++){
      for(int bin_y = -nBinsHisto_+1; bin_y<int(nBinsHisto_); bin_y++){

	//Build seed with four adjacent bins
	float MIPT = histoMultiCluster[{ z_side, bin_x, bin_y }];
	MIPT +=  histoMultiCluster[{ z_side, bin_x+1, bin_y }];
	MIPT += histoMultiCluster[{ z_side, bin_x, bin_y+1 }];
	MIPT += histoMultiCluster[{ z_side, bin_x+1, bin_y+1 }];

	//Compare with corners
	float MIPT_NW = histoMultiCluster[{ z_side, bin_x, bin_y }];
	MIPT_NW += histoMultiCluster[{ z_side, bin_x-1, bin_y }];
	MIPT_NW += histoMultiCluster[{ z_side, bin_x, bin_y-1 }];
	MIPT_NW += histoMultiCluster[{ z_side, bin_x-1, bin_y-1 }];

	float MIPT_SW = histoMultiCluster[{ z_side, bin_x, bin_y+1 }];
	MIPT_SW += histoMultiCluster[{ z_side, bin_x, bin_y+2 }];
	MIPT_SW += histoMultiCluster[{ z_side, bin_x-1, bin_y+1 }];
	MIPT_SW += histoMultiCluster[{ z_side, bin_x-1, bin_y+2 }];

	float MIPT_SE = histoMultiCluster[{ z_side, bin_x+1, bin_y+1 }];
	MIPT_SE += histoMultiCluster[{ z_side, bin_x+1, bin_y+2 }];
	MIPT_SE += histoMultiCluster[{ z_side, bin_x+2, bin_y+1 }];
	MIPT_SE += histoMultiCluster[{ z_side, bin_x+2, bin_y+2 }];

	float MIPT_NE = histoMultiCluster[{ z_side, bin_x+1, bin_y }];
	MIPT_NE += histoMultiCluster[{ z_side, bin_x+1, bin_y-1 }];
	MIPT_NE += histoMultiCluster[{ z_side, bin_x+2, bin_y }];
	MIPT_NE += histoMultiCluster[{ z_side, bin_x+2, bin_y-1 }];

	bool isSeed = MIPT>=MIPT_NW && MIPT>=MIPT_SW && MIPT>MIPT_NE && MIPT>MIPT_SE;

	if(isSeed){
	  float x_seed = (bin_x+0.5)/float(nBinsHisto_)*kXYOverZMax_;
	  float y_seed = (bin_y+0.5)/float(nBinsHisto_)*kXYOverZMax_;
	  GlobalPoint seed(x_seed,y_seed,z_side);
	  seedPositions.emplace_back(seed);
	}

      }
    }
  }


  //Clustering
  for(std::vector<edm::Ptr<l1t::HGCalCluster>>::const_iterator clu = clustersPtrs.begin(); clu != clustersPtrs.end(); ++clu){

    int iseed=0;
    vector<int> tcPertinentSeeds;
    for( const auto& seed : seedPositions ){
      if( this->isPertinent(**clu, seed, dr_) ){
	tcPertinentSeeds.push_back(iseed);
      }
      ++iseed;
    }

    unsigned minDist = 1;
    unsigned targetSeed = 0;
    for( int iseed : tcPertinentSeeds ){
      GlobalPoint seed = seedPositions[iseed];
      double d = (seed - (**clu).centreProj() ).mag() ;
      if( d < minDist ){
	minDist = d;
	targetSeed = iseed;
      }
    }

    mapSeedMulticluster[targetSeed].addConstituent( *clu );

  }


  for(auto mclu : mapSeedMulticluster) multiclustersTmp.emplace_back(mclu.second);


  if(superClustering_){

    for(unsigned int i_mclu1=0; i_mclu1<multiclustersTmp.size(); i_mclu1++){

      l1t::HGCalMulticluster mclu1 = multiclustersTmp[i_mclu1];

      for(unsigned int i_mclu2=0; i_mclu2<multiclustersTmp.size(); i_mclu2++){
	if(i_mclu1==i_mclu2) continue;

	l1t::HGCalMulticluster mclu2 = multiclustersTmp[i_mclu2];

	if( this->isPertinent(mclu1, mclu2, dr_) ){

	  const std::unordered_map<uint32_t, edm::Ptr<l1t::HGCalCluster>>& clusters = mclu2.constituents();
	  for(const auto& id_cluster : clusters)  multiclustersTmp[i_mclu1].addConstituent(id_cluster.second);
	  multiclustersTmp.erase(multiclustersTmp.begin()+i_mclu2);
	  i_mclu2--;

	}

      }

    }

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
