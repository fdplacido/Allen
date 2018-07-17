#include "../include/run_PatPV_CPU.h"


/*
XYZPoint& seedPoint,
              std::vector<Track*>& rTracks,
             Vertex& vtx,
             std::vector<Track*>& tracks2remove

*/





bool reconstructMultiPVFromTracks( VeloState * tracks2use,
                                                       Vertex * outvtxvec, int host_number_of_tracks_pinned,
  uint * number_of_vertex) 
{
  

  VeloState * rtracks = tracks2use;

  //outvtxvec.clear();


  PVSeedTool seedtool;
  double m_beamSpotX = 0.02;
  double m_beamSpotY = -0.16;
  XYZPoint beamspot{m_beamSpotX, m_beamSpotY, 0.};
  

    


  int nvtx_before = -1;
  int nvtx_after  =  0;
  //for (int i = 0; i < 5 ; i++) {
  //do we really need this loop?
  //while ( nvtx_after > nvtx_before ) {
    nvtx_before = nvtx_after;
    // reconstruct vertices


  AdaptivePV3DFitter fitter;
  std::vector<XYZPoint> seeds = seedtool.getSeeds(rtracks, beamspot, host_number_of_tracks_pinned);
    for ( auto seed : seeds) {
      Vertex recvtx;


      //VeloState * tracks2remove;
      std::vector<VeloState> tracks2remove;
      // fitting
      bool scvfit = fitter.fitVertex( seed, rtracks, recvtx, tracks2remove, host_number_of_tracks_pinned);
      if (!scvfit) continue;
      
      

      
      outvtxvec[nvtx_after] = recvtx;
      nvtx_after++;
      //removeTracks(rtracks, tracks2remove);
    }//iterate on seeds
    *number_of_vertex = nvtx_after;
  //}//iterate on vtx

  return true;

}




int run_PatPV_on_CPU (
  VeloState * host_velo_states,
  int * host_accumulated_tracks,
  uint* host_velo_track_hit_number_pinned,
  VeloTracking::Hit<true>* host_velo_track_hits_pinned,
  int * host_number_of_tracks_pinned,
  const int &number_of_events,
  Vertex * outvtxvec,
  uint * number_of_vertex
) {

XYZPoint beamspot(0.,0.,0.);
PVSeedTool seedtool;
//std:std::vector<XYZPoint> seeds = seedtool.getSeeds(host_velo_states, beamspot, *host_number_of_tracks_pinned);

/*
AdaptivePV3DFitter fitter;
Vertex recvtx;
std::vector<VeloState> tracks2remove;
XYZPoint seed = seeds.at(0);
fitter.fitVertex(seed, host_velo_states, recvtx, tracks2remove, number_of_events);
*/
//Vertex  outvtxvec[100];
//std::vector<VeloState> velostate_vec;

//for(int i = 0; i < *host_number_of_tracks_pinned; i++)  velostate_vec.push_back(host_velo_states[i]); 


for(int i_event = 0; i_event < number_of_events; i_event++) {

  int number_of_tracks = host_number_of_tracks_pinned[i_event];
 VeloState * state_base_pointer = host_velo_states + 2 * host_accumulated_tracks[i_event];
VeloState  kalman_states[number_of_tracks];

//recovers previusoly found vertices in first event
for(int i = 0; i < number_of_tracks; i++) kalman_states[i] = state_base_pointer[2*i ];
std::cout << "least: " << kalman_states[1].x <<std::endl;
std::cout << "least: " << kalman_states[1].y <<std::endl;
std::cout << "least: " << kalman_states[1].z <<std::endl;
std::cout << "least: " << kalman_states[1].tx <<std::endl;
std::cout << "least: " << kalman_states[1].ty <<std::endl;
std::cout << "least: " << kalman_states[1].c00 <<std::endl;
std::cout << "least: " << kalman_states[1].c20 <<std::endl;
std::cout << "least: " << kalman_states[1].c22 <<std::endl;
std::cout << "least: " << kalman_states[1].c11 <<std::endl;
std::cout << "least: " << kalman_states[1].c31 <<std::endl;
std::cout << "least: " << kalman_states[1].c33 <<std::endl;

  //not workign yet
for(int i = 0; i < number_of_tracks; i++) kalman_states[i] = state_base_pointer[2*i +1];
std::cout << "kalman: " << kalman_states[1].x <<std::endl;
std::cout << "kalman: " << kalman_states[1].y <<std::endl;
std::cout << "kalman: " << kalman_states[1].z <<std::endl;
std::cout << "kalman: " << kalman_states[1].tx <<std::endl;
std::cout << "kalman: " << kalman_states[1].ty <<std::endl;
std::cout << "kalman: " << kalman_states[1].c00 <<std::endl;
std::cout << "kalman: " << kalman_states[1].c20 <<std::endl;
std::cout << "kalman: " << kalman_states[1].c22 <<std::endl;
std::cout << "kalman: " << kalman_states[1].c11 <<std::endl;
std::cout << "kalman: " << kalman_states[1].c31 <<std::endl;
std::cout << "least: " << kalman_states[1].c33 <<std::endl;

reconstructMultiPVFromTracks(kalman_states, outvtxvec, host_number_of_tracks_pinned[i_event], number_of_vertex);
}

  return 0;
}