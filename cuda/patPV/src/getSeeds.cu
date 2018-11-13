#include "getSeeds.cuh"
//simplficiations: no tracks2disable




__device__ double zCloseBeam( Velo::State track, const PatPV::XYZPoint& beamspot) {

  PatPV::XYZPoint tpoint(track.x, track.y, track.z);
  PatPV::XYZPoint tdir(track.tx, track.ty, 1.);

  double wx = ( 1. + tdir.x * tdir.x ) / track.c00;
  double wy = ( 1. + tdir.y * tdir.y ) / track.c11;

  double x0 = tpoint.x - tpoint.z * tdir.x - beamspot.x;
  double y0 = tpoint.y - tpoint.z * tdir.y - beamspot.y;
  double den = wx * tdir.x * tdir.x + wy * tdir.y * tdir.y;
  double zAtBeam = - ( wx * x0 * tdir.x + wy * y0 * tdir.y ) / den ;

  double xb = tpoint.x + tdir.x * ( zAtBeam - tpoint.z ) - beamspot.x;
  double yb = tpoint.y + tdir.y * ( zAtBeam - tpoint.z ) - beamspot.y;
  double r2AtBeam = xb*xb + yb*yb ;

  return r2AtBeam < 0.5*0.5 ? zAtBeam : 10e8;
}



__device__ void errorForPVSeedFinding(double tx, double ty, double &sigz2)  {

    // the seeding results depend weakly on this eror parametrization

    double pMean = 3000.; // unit: MeV

    double tanTheta2 =  tx * tx + ty * ty;
    double sinTheta2 =  tanTheta2 / ( 1. + tanTheta2 );

    // assume that first hit in VD at 8 mm
    double distr        = 8.; // unit: mm
    double dist2        = distr*distr/sinTheta2;
    double sigma_ms2    = PatPV::mcu_scatCons * PatPV::mcu_scatCons * dist2 / (pMean*pMean);
    double fslope2      = 0.0005*0.0005;
    double sigma_slope2 = fslope2*dist2;

    sigz2 = (sigma_ms2 + sigma_slope2) / sinTheta2;
    if(sigz2 == 0) sigz2 = 100.;

}




 __global__ void getSeeds(
  uint* dev_kalmanvelo_states,
  int * dev_atomics_storage,
  uint* dev_velo_track_hit_number,
  PatPV::XYZPoint * dev_seeds,
  uint * dev_number_seed) {

  PatPV::XYZPoint beamspot;
  beamspot.x = 0;
  beamspot.y = 0;
  beamspot.z = 0;

  const uint number_of_events = gridDim.x;
  const uint event_number = blockIdx.x;

  const Velo::Consolidated::Tracks velo_tracks {(uint*) dev_atomics_storage, dev_velo_track_hit_number, event_number, number_of_events};
  const Velo::Consolidated::States velo_states {dev_kalmanvelo_states, velo_tracks.total_number_of_tracks};
  const uint number_of_tracks_event = velo_tracks.number_of_tracks(event_number);
  const uint event_tracks_offset = velo_tracks.tracks_offset(event_number);

  



  

    PatPV::vtxCluster  vclusters[VeloTracking::max_tracks];



  int counter_number_of_clusters = 0;
  for (int i = 0; i < number_of_tracks_event; i++) {

    
    double sigsq;
    double zclu;
    Velo::State trk = velo_states.get(event_tracks_offset + i);


    zclu = zCloseBeam(trk,beamspot);
    errorForPVSeedFinding(trk.tx, trk.ty,sigsq);

    if ( fabs(zclu)>2000.) continue;
    PatPV::vtxCluster clu;
    clu.z = zclu;
    clu.sigsq = sigsq;
    clu.sigsqmin = clu.sigsq;
    clu.ntracks = 1;
    vclusters[counter_number_of_clusters] = clu;
    
    counter_number_of_clusters++;

  }

  double  zseeds[VeloTracking::max_tracks];

  int number_final_clusters = findClusters(vclusters, zseeds, counter_number_of_clusters);

  for(int i = 0; i < number_final_clusters; i++) dev_seeds[event_number * PatPV::max_number_vertices + i] = PatPV::XYZPoint{ beamspot.x, beamspot.y, zseeds[i]};
  

  dev_number_seed[event_number] = number_final_clusters;


 };

 __device__ int findClusters(PatPV::vtxCluster * vclus, double * zclusters, int number_of_clusters)  {


  
  

  
  for(int i = 0; i < number_of_clusters; i++) {
    vclus[i].sigsq *= PatPV::mcu_factorToIncreaseErrors*PatPV::mcu_factorToIncreaseErrors; // blow up errors
    vclus[i].sigsqmin = vclus[i].sigsq;
  }


  //maybe sort in z before merging? -> does not seem to help

  bool no_merges = false;
  while(!no_merges) {
    //reset merged flags
    for (int j = 0; j < number_of_clusters; j++) vclus[j].merged = false;

    no_merges = true;
    for(int index_cluster = 0; index_cluster < number_of_clusters - 1; index_cluster++) {


      //skip cluster which have already been merged
      if(vclus[index_cluster].ntracks == 0) continue;

   //sorting by chi2dist seems to increase efficiency in nominal code


      
      for(int index_second_cluster = 0; index_second_cluster < number_of_clusters ; index_second_cluster++) {
        if(vclus[index_second_cluster].merged || vclus[index_cluster].merged) continue;
        //skip cluster which have already been merged
        if(vclus[index_second_cluster].ntracks == 0) continue;
        if(index_cluster == index_second_cluster) continue;
        double z1 = vclus[index_cluster].z;
        double z2 = vclus[index_second_cluster].z;
        double s1 = vclus[index_cluster].sigsq;
        double s2 = vclus[index_second_cluster].sigsq;
        double s1min = vclus[index_cluster].sigsqmin;
        double s2min = vclus[index_second_cluster].sigsqmin;
        double sigsqmin = s1min;
        if(s2min<s1min) sigsqmin = s2min;


        double zdist = z1 - z2;
        double chi2dist = zdist*zdist/(s1+s2);
        //merge if chi2dist is smaller than max
        if (chi2dist<PatPV::mcu_maxChi2Merge ) {
          no_merges = false;
          double w_inv = (s1*s2/(s1+s2));
          double zmerge = w_inv*(z1/s1+z2/s2);

          vclus[index_cluster].z        = zmerge;
          vclus[index_cluster].sigsq    = w_inv;
          vclus[index_cluster].sigsqmin = sigsqmin;
          vclus[index_cluster].ntracks += vclus[index_second_cluster].ntracks;
          vclus[index_second_cluster].ntracks = 0;  // mark second cluster as used
          vclus[index_cluster].merged = true;
          vclus[index_second_cluster].merged = true;

          //break;
        } 
      }
    }
  }
 
  

  int return_number_of_clusters = 0;
  //count final number of clusters
  PatPV::vtxCluster pvclus[VeloTracking::max_tracks];
  for(int i = 0; i < number_of_clusters; i++) {
    if(vclus[i].ntracks != 0)     {pvclus[return_number_of_clusters] = vclus[i]; return_number_of_clusters++;}
  } 


  //clean up clusters, do we gain much from this?

  // Select good clusters.

  int number_good_clusters = 0;

  for(int index = 0; index < return_number_of_clusters; index++) {

    int n_tracks_close = 0;
    for(int i = 0; i < number_of_clusters; i++) if(fabs(vclus[i].z - pvclus[index].z ) < PatPV::mcu_dzCloseTracksInCluster ) n_tracks_close++;
  

    double dist_to_closest = 1000000.;
    if(return_number_of_clusters > 1) {
      for(int index2 = 0; index2 < return_number_of_clusters; index2++) {
        if( index!=index2 && ( fabs( pvclus[index2].z - pvclus[index].z) < dist_to_closest) )  dist_to_closest = fabs( pvclus[index2].z - pvclus[index].z);
      }
    }

    // ratio to remove clusters made of one low error track and many large error ones
    double rat = pvclus[index].sigsq/pvclus[index].sigsqmin;
    bool igood = false;
    int ntracks = pvclus[index].ntracks;
    if( ntracks >= PatPV::mcu_minClusterMult ) {
      if( dist_to_closest>10. && rat<0.95) igood=true;
      if( ntracks >= PatPV::mcu_highMult && rat < PatPV::mcu_ratioSig2HighMult)  igood=true;
      if( ntracks <  PatPV::mcu_highMult && rat < PatPV::mcu_ratioSig2LowMult )  igood=true;
    }
    // veto
    if( n_tracks_close < PatPV::mcu_minCloseTracksInCluster ) igood = false;
    if(igood) {zclusters[number_good_clusters] = pvclus[index].z; number_good_clusters++;}


  }



  //return return_number_of_clusters;
  return number_good_clusters;

}