

#include "L1Trigger/L1THGCal/interface/be_algorithms/HGCalMulticlusteringImpl.h"
#include "L1Trigger/L1THGCal/interface/be_algorithms/HGCalShowerShape.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

HGCalMulticlusteringImpl::HGCalMulticlusteringImpl( const edm::ParameterSet& conf ) :
    dr_(conf.getParameter<double>("dR_multicluster")),
    ptC3dThreshold_(conf.getParameter<double>("minPt_multicluster")),
    multiclusterAlgoType_(conf.getParameter<string>("type_multicluster")),
    distDbscan_(conf.getParameter<double>("dist_dbscan_multicluster")),
    minNDbscan_(conf.getParameter<unsigned>("minN_dbscan_multicluster")),
    nBinsRHisto_(conf.getParameter<unsigned>("nBins_R_histo_multicluster")),
    nBinsPhiHisto_(conf.getParameter<unsigned>("nBins_Phi_histo_multicluster")),
    superClustering_(conf.getParameter<bool>("superClustering"))
{    
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster dR for Near Neighbour search: " << dr_;  
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster minimum transverse-momentum: " << ptC3dThreshold_;
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster DBSCAN Clustering distance: " << distDbscan_;
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster clustering min number of subclusters: " << minNDbscan_;
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster number of R-bins for the histo algorithm: " << nBinsRHisto_<<endl;
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster number of Phi-bins for the histo algorithm: " << nBinsPhiHisto_<<endl;
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster superclustering flag: " << superClustering_<<endl;
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




float HGCalMulticlusteringImpl::dR( const l1t::HGCalCluster & clu,
				   const GlobalPoint & seed) const
{

    Basic3DVector<float> seed_3dv( seed );
    GlobalPoint seed_proj = GlobalPoint( seed_3dv / seed.z() );
    return (seed_proj - clu.centreProj() ).mag();

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

        double minDist = dr_;
        int targetMulticlu = -1;

        for(unsigned imclu=0; imclu<multiclustersTmp.size(); imclu++){

            if(!this->isPertinent(**clu, multiclustersTmp.at(imclu), dr_)) continue;

            double d = ( multiclustersTmp.at(imclu).centreProj() - (*clu)->centreProj() ).mag() ;
            if(d<minDist){
                minDist = d;
                targetMulticlu = int(imclu);
            }
        }

        if(targetMulticlu<0) multiclustersTmp.emplace_back( *clu );
        else multiclustersTmp.at( targetMulticlu ).addConstituent( *clu );

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




std::map<std::vector<int>,float> HGCalMulticlusteringImpl::fillHistoClusters( const std::vector<edm::Ptr<l1t::HGCalCluster>> & clustersPtrs ){


  std::map<std::vector<int>,float> histoClusters; //key[0] = z.side(), key[1] = bin_R, key[2] = bin_phi

  for(std::vector<edm::Ptr<l1t::HGCalCluster>>::const_iterator clu = clustersPtrs.begin(); clu != clustersPtrs.end(); ++clu){

    float ROverZ = sqrt( pow((**clu).centreProj().x(),2) + pow((**clu).centreProj().y(),2) );
    int bin_R = int( (ROverZ-kROverZMin_) * nBinsRHisto_ / (kROverZMax_-kROverZMin_) );
    int bin_phi = int( (reco::reduceRange((**clu).phi())+M_PI) * nBinsPhiHisto_ / (2*M_PI) );

    std::vector<int> key = { (**clu).zside(), bin_R, bin_phi };
    histoClusters[key]+=(**clu).mipPt();

  }

  return histoClusters;

}




vector<unsigned>  HGCalMulticlusteringImpl::binSumsPolarHisto( const unsigned nBins ){

  vector<unsigned> bins;

  //Hardcoded for 36 for now, to be fixed
  if(nBins==36){
    bins = {13,                                           // 0
	    11, 11, 11,                                   // 1 - 3
	    9, 9, 9,                                      // 4 - 6
	    7, 7, 7, 7, 7, 7,                             // 7 - 12
	    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,  //13 - 27
	    3, 3, 3, 3, 3, 3, 3, 3                        //28 - 35
    };
  }

  return bins;

}



std::map<std::vector<int>,float> HGCalMulticlusteringImpl::fillSmoothPhiHistoClusters( std::map<std::vector<int>,float> histoClusters,
										       const vector<unsigned> binSums){

  std::map<std::vector<int>,float> histoSumPhiClusters; //key[0] = z.side(), key[1] = bin_R, key[2] = bin_phi

  for(int z_side = -1; z_side<2; z_side++){
    if(z_side==0) continue;

    for(int bin_R = 0; bin_R<int(nBinsRHisto_); bin_R++){

      int nBinsSide = (binSums[bin_R]-1)/2;

      for(int bin_phi = 0; bin_phi<int(nBinsPhiHisto_); bin_phi++){

	float content = histoClusters[{z_side,bin_R,bin_phi}];

	for(int bin_phi2=1; bin_phi2<=nBinsSide; bin_phi2++ ){

	  int binToSumLeft = bin_phi - bin_phi2;
	  if( binToSumLeft<0 ) binToSumLeft += nBinsPhiHisto_;
	  int binToSumRight = bin_phi + bin_phi2;
	  if( binToSumRight>=int(nBinsPhiHisto_) ) binToSumRight -= nBinsPhiHisto_;

	  content += histoClusters[{z_side,bin_R,binToSumLeft}] / pow(2,bin_phi2); // quadratic kernel
	  content += histoClusters[{z_side,bin_R,binToSumRight}] / pow(2,bin_phi2); // quadratic kernel

	}

	float R1 = kROverZMin_ + bin_R*(kROverZMax_-kROverZMin_);
	float R2 = R1 + (kROverZMax_-kROverZMin_);
	double area = 0.5 * (pow(R2,2)-pow(R1,2)) * (1+0.5*(1-pow(0.5,nBinsSide))); // Takes into account different area of bins in different R-rings + sum of quadratic weights used
	histoSumPhiClusters[{z_side,bin_R,bin_phi}] = content/area;

      }

    }

  }

  return histoSumPhiClusters;

}



std::map<std::vector<int>,float> HGCalMulticlusteringImpl::fillSmoothRPhiHistoClusters( std::map<std::vector<int>,float> histoClusters){

 std::map<std::vector<int>,float> histoSumRPhiClusters; //key[0] = z.side(), key[1] = bin_R, key[2] = bin_phi

  for(int z_side = -1; z_side<2; z_side++){
    if(z_side==0) continue;

    for(int bin_R = 0; bin_R<int(nBinsRHisto_); bin_R++){

      float weight = (bin_R==0 || bin_R==int(nBinsRHisto_)-1) ? 1.5 : 2.;

      for(int bin_phi = 0; bin_phi<int(nBinsPhiHisto_); bin_phi++){

	float content = histoClusters[{z_side,bin_R,bin_phi}];
	float contentDown = histoClusters[{z_side,bin_R-1,bin_phi}]; //Non-allocated elements in maps return default 0 value
	float contentUp = histoClusters[{z_side,bin_R+1,bin_phi}];

	histoSumRPhiClusters[{z_side,bin_R,bin_phi}] = (content + 0.5*contentDown + 0.5*contentUp)/weight;

      }

    }

  }

  return histoSumRPhiClusters;

}



std::vector<GlobalPoint> HGCalMulticlusteringImpl::computeSeeds( std::map<std::vector<int>,float> histoClusters,
								 const vector<unsigned> binSums ){

  std::vector<GlobalPoint> seedPositions;

  for(int z_side = -1; z_side<2; z_side++){
    if(z_side==0) continue;

    for(int bin_R = 0; bin_R<int(nBinsRHisto_); bin_R++){

      int nBinsSide = (binSums[bin_R]-1)/2;
      int nBinsSideDown = bin_R==0 ? 0 : (binSums[bin_R-1]-1)/2;
      int nBinsSideUp = bin_R==int(nBinsRHisto_)-1 ? 0 : (binSums[bin_R+1]-1)/2;

      for(int bin_phi = 0; bin_phi<int(nBinsPhiHisto_); bin_phi++){

	float MIPT_seed = histoClusters[{z_side,bin_R,bin_phi}];
	bool isMax = MIPT_seed>0;

	float MIPT_S = histoClusters[{z_side,bin_R+1,bin_phi}];
	float MIPT_N = histoClusters[{z_side,bin_R-1,bin_phi}];

	isMax &= MIPT_seed>=MIPT_S;
	isMax &= MIPT_seed>MIPT_N;

	for(int bin_phi2=1; bin_phi2<=nBinsSideDown; bin_phi2++){ //Down side has always larger number of bins
	  int binLeft = bin_phi - bin_phi2;
	  if( binLeft<0 ) binLeft += nBinsPhiHisto_;
	  int binRight = bin_phi + bin_phi2;
	  if( binRight>=int(nBinsPhiHisto_) ) binRight -= nBinsPhiHisto_;

	  float MIPT_W = bin_phi2<=nBinsSide ? histoClusters[{z_side,bin_R,binLeft}] : 0;
	  float MIPT_E = bin_phi2<=nBinsSide ? histoClusters[{z_side,bin_R,binRight}] : 0;
	  float MIPT_SW = bin_phi2<=nBinsSideDown ? histoClusters[{z_side,bin_R-1,binLeft}] : 0;
	  float MIPT_SE = bin_phi2<=nBinsSideDown ? histoClusters[{z_side,bin_R-1,binRight}] : 0;
	  float MIPT_NW = bin_phi2<=nBinsSideUp ? histoClusters[{z_side,bin_R+1,binLeft}] : 0;
	  float MIPT_NE = bin_phi2<=nBinsSideUp ? histoClusters[{z_side,bin_R+1,binRight}] : 0;

	  isMax &= MIPT_seed>=MIPT_E;
	  isMax &= MIPT_seed>=MIPT_SE;
	  isMax &= MIPT_seed>=MIPT_NE;
	  isMax &= MIPT_seed>MIPT_W;
	  isMax &= MIPT_seed>MIPT_SW;
	  isMax &= MIPT_seed>MIPT_NW;

	}

	if(isMax){
	  float ROverZ_seed = kROverZMin_ + (bin_R+0.5) * (kROverZMax_-kROverZMin_)/nBinsRHisto_;
	  float phi_seed = -M_PI + (bin_phi+0.5) * 2*M_PI/nBinsPhiHisto_;
	  float x_seed = ROverZ_seed*cos(phi_seed);
	  float y_seed = ROverZ_seed*sin(phi_seed);
	  seedPositions.emplace_back(x_seed,y_seed,z_side);
	}

      }

    }

  }

  return seedPositions;

}


void HGCalMulticlusteringImpl::superClusteringDR(std::vector<l1t::HGCalMulticluster>& multiclusters){

  for(unsigned int i_mclu1=0; i_mclu1<multiclusters.size(); i_mclu1++){

    l1t::HGCalMulticluster mclu1 = multiclusters[i_mclu1];

    for(unsigned int i_mclu2=0; i_mclu2<multiclusters.size(); i_mclu2++){
      if(i_mclu1==i_mclu2) continue;

      l1t::HGCalMulticluster mclu2 = multiclusters[i_mclu2];

      if( this->isPertinent(mclu1, mclu2, dr_) ){

	const std::unordered_map<uint32_t, edm::Ptr<l1t::HGCalCluster>>& clusters = mclu2.constituents();
	for(const auto& id_cluster : clusters)  multiclusters[i_mclu1].addConstituent(id_cluster.second);
	multiclusters.erase(multiclusters.begin()+i_mclu2);
	i_mclu2--;

      }

    }

  }

  return;

}




void HGCalMulticlusteringImpl::clusterizePolarHistoDR( const std::vector<edm::Ptr<l1t::HGCalCluster>> & clustersPtrs,
						       l1t::HGCalMulticlusterBxCollection & multiclusters,
						       const HGCalTriggerGeometryBase & triggerGeometry)
{


  vector<unsigned> nBinsToSum = binSumsPolarHisto( nBinsRHisto_ );

  std::map<std::vector<int>,float> histoCluster = fillHistoClusters(clustersPtrs); //key[0] = z.side(), key[1] = bin_R, key[2] = bin_phi, content = MIPT
  std::map<std::vector<int>,float> smoothPhiHistoCluster = fillSmoothPhiHistoClusters(histoCluster,nBinsToSum); //content = MIPT-density, smoothened along phi
  std::map<std::vector<int>,float> smoothRPhiHistoCluster = fillSmoothRPhiHistoClusters(histoCluster); //content = MIPT-density, smoothened along R-phi

  // Seed positions determined with maximum finding
  std::vector<GlobalPoint> seedPositions = computeSeeds(smoothRPhiHistoCluster,nBinsToSum);


  //Clustering
  std::map<int,l1t::HGCalMulticluster> mapSeedMulticluster;
  std::vector<l1t::HGCalMulticluster> multiclustersTmp;

  for(std::vector<edm::Ptr<l1t::HGCalCluster>>::const_iterator clu = clustersPtrs.begin(); clu != clustersPtrs.end(); ++clu){

    HGCalDetId cluDetId( (**clu).detId() );
    int z_side = cluDetId.zside();

    double minDist = dr_;
    int targetSeed = -1;

    for( unsigned int iseed=0; iseed<seedPositions.size(); iseed++ ){

      if( z_side*seedPositions[iseed].z()<0) continue;

      double d = this->dR(**clu, seedPositions[iseed]);

      if(d<minDist){
	minDist = d;
	targetSeed = iseed;
      }
    }

    if(targetSeed<0) continue;

    if(mapSeedMulticluster[targetSeed].size()==0){
      l1t::HGCalMulticluster newMclu(*clu);
      mapSeedMulticluster[targetSeed] = newMclu;
    }
    else mapSeedMulticluster[targetSeed].addConstituent(*clu);

  }

  for(auto mclu : mapSeedMulticluster) multiclustersTmp.emplace_back(mclu.second);

  if(superClustering_) superClusteringDR(multiclustersTmp);

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
