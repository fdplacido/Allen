#include "RunMomentumForwardCPU.h"

#ifdef WITH_ROOT
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#endif

// straight line extrapolation of MiniState to other z position
MiniState state_at_z(const MiniState state, const float z) {
  MiniState extrap_state;
  extrap_state.tx = state.tx;
  extrap_state.ty = state.ty;
  extrap_state.x = state.x + (z-state.z) * state.tx;
  extrap_state.y = state.y + (z-state.z) * state.ty;
  extrap_state.z = z;
  return extrap_state;
}

// straight line extrapolation of y to other z position
float y_at_z(const MiniState state, const float z) {
  float yf = state.y + (z-state.z) * state.ty;
  return yf;
}

int run_momentum_forward_on_CPU(
  SciFi::TrackHits* host_scifi_tracks,
  int* host_scifi_n_tracks,
  const uint* host_scifi_hits,
  const uint* host_scifi_hit_count,
  const char* host_scifi_geometry,
  const std::array<float, 9>& host_inv_clus_res,
  const SciFi::Parameters& scifi_params,
  const uint* host_velo_tracks_atomics,
  const uint* host_velo_track_hit_number,
  const char* host_velo_states,
  const int * host_atomics_ut,
  const uint* host_ut_track_hit_number,
  const float* host_ut_qop,
  const float* host_ut_x,
  const float* host_ut_tx,
  const float* host_ut_z,
  const uint* host_ut_track_velo_indices,
  const std::vector< std::vector< std::vector< uint32_t > > > scifi_ids_ut_tracks,
  const std::vector< std::vector< float > > p_events,
  const uint number_of_events
) {

#ifdef WITH_ROOT
  // Histograms only for checking and debugging
  TFile *f = new TFile("../output/scifi.root", "RECREATE");
  TTree *t_Forward_tracks = new TTree("Forward_tracks", "Forward_tracks");
  TTree *t_statistics = new TTree("statistics", "statistics");
  TTree *t_scifi_hits = new TTree("scifi_hits","scifi_hits");
  TTree *t_extrap = new TTree("extrap","extrap");
  TTree *t_ut_tracks = new TTree("ut_tracks","ut_tracks");
  uint planeCode, LHCbID;
  float x0, z0, w, dxdy, dzdy, yMin, yMax;
  float qop;
  int n_tracks;
  float state_x, state_y, state_z, state_tx, state_ty;
  float xf, yf, txf, tyf, der_xf_qop, qop_update;
  float res_x_0, res_x_3, ut_qop, dx, x_extrap;
  float UT_x, UT_y, UT_z, UT_tx, UT_ty;
  float velo_x_extrap, velo_tx;
  int n_hits_in_window_0, n_hits_in_window_3, n_x_combinations;
  float p_diff_before_update, p_diff_after_update;

  t_Forward_tracks->Branch("qop", &qop);
  t_Forward_tracks->Branch("state_x", &state_x);
  t_Forward_tracks->Branch("state_y", &state_y);
  t_Forward_tracks->Branch("state_z", &state_z);
  t_Forward_tracks->Branch("state_tx", &state_tx);
  t_Forward_tracks->Branch("state_ty", &state_ty);

  t_statistics->Branch("n_tracks", &n_tracks);

  t_scifi_hits->Branch("planeCode", &planeCode);
  t_scifi_hits->Branch("LHCbID", &LHCbID);
  t_scifi_hits->Branch("x0", &x0);
  t_scifi_hits->Branch("z0", &z0);
  t_scifi_hits->Branch("w", &w);
  t_scifi_hits->Branch("dxdy", &dxdy);
  t_scifi_hits->Branch("dzdy", &dzdy);
  t_scifi_hits->Branch("yMin", &yMin);
  t_scifi_hits->Branch("yMax", &yMax);

  t_extrap->Branch("xf", &xf);
  t_extrap->Branch("yf", &yf);
  t_extrap->Branch("txf", &txf);
  t_extrap->Branch("tyf", &tyf);
  t_extrap->Branch("der_xf_qop", &der_xf_qop);
  t_extrap->Branch("res_x_0", &res_x_0);
  t_extrap->Branch("res_x_3", &res_x_3);
  t_extrap->Branch("x_extrap", &x_extrap);
  t_extrap->Branch("dx", &dx);
  t_extrap->Branch("ut_qop", &ut_qop);
  t_extrap->Branch("qop_update", &qop_update);
  t_extrap->Branch("p_diff_before_update", &p_diff_before_update);
  t_extrap->Branch("p_diff_after_update", &p_diff_after_update);
  t_extrap->Branch("n_hits_in_window_0", &n_hits_in_window_0);
  t_extrap->Branch("n_hits_in_window_3", &n_hits_in_window_3);
  t_extrap->Branch("n_x_combinations", &n_x_combinations);

  t_ut_tracks->Branch("ut_x", &UT_x);
  t_ut_tracks->Branch("ut_y", &UT_y);
  t_ut_tracks->Branch("ut_z", &UT_z);
  t_ut_tracks->Branch("ut_tx", &UT_tx);
  t_ut_tracks->Branch("ut_ty", &UT_ty);
  t_ut_tracks->Branch("velo_x_extrap", &velo_x_extrap);
  t_ut_tracks->Branch("velo_tx", &velo_tx);
  t_ut_tracks->Branch("ut_qop", &ut_qop);
#endif

  for ( uint i_event = 0; i_event < number_of_events; ++i_event ) {

    // Velo consolidated types
    const Velo::Consolidated::Tracks velo_tracks {(uint*) host_velo_tracks_atomics, (uint*)host_velo_track_hit_number, i_event, number_of_events};
    const uint velo_event_tracks_offset = velo_tracks.tracks_offset(i_event);
    const Velo::Consolidated::States velo_states {(char*)host_velo_states, velo_tracks.total_number_of_tracks};

    // UT consolidated types
    UT::Consolidated::Tracks ut_tracks {
     (uint*)host_atomics_ut,
     (uint*)host_ut_track_hit_number,
     (float*)host_ut_qop,
     (uint*)host_ut_track_velo_indices,
      i_event,
      number_of_events
      };
    const int n_veloUT_tracks_event = ut_tracks.number_of_tracks(i_event);
    const int ut_event_tracks_offset = ut_tracks.tracks_offset(i_event);
    
    // SciFi non-consolidated types
    int* n_forward_tracks = host_scifi_n_tracks + i_event;
    SciFi::TrackHits* scifi_tracks_event = host_scifi_tracks + i_event * SciFi::Constants::max_tracks;

    const uint total_number_of_hits = host_scifi_hit_count[number_of_events * SciFi::Constants::n_mat_groups_and_mats]; 
    SciFi::HitCount scifi_hit_count {(uint32_t*)host_scifi_hit_count, i_event};
    
    const SciFi::SciFiGeometry scifi_geometry(host_scifi_geometry);
    
    SciFi::Hits scifi_hits(
     (uint*)host_scifi_hits,
     total_number_of_hits,
     &scifi_geometry, 
     reinterpret_cast<const float*>(host_inv_clus_res.data()));

#ifdef WITH_ROOT
    // store hit variables in tree
    for(size_t zone = 0; zone < SciFi::Constants::n_zones; zone++) {
      const auto zone_offset = scifi_hit_count.zone_offset(zone);
      for(size_t hit = 0; hit < scifi_hit_count.zone_number_of_hits(zone); hit++) {
        const auto hit_offset = zone_offset + hit;

        planeCode = scifi_hits.planeCode(hit_offset);
        LHCbID = scifi_hits.LHCbID(hit_offset);
        x0 = scifi_hits.x0[hit_offset];
        z0 = scifi_hits.z0[hit_offset];
        w  = scifi_hits.w(hit_offset);
        dxdy = scifi_hits.dxdy(hit_offset);
        dzdy = scifi_hits.dzdy(hit_offset);
        yMin = scifi_hits.yMin(hit_offset);
        yMax = scifi_hits.yMax(hit_offset);
        t_scifi_hits->Fill();
      }
    }
#endif

    /* etrapolation to first SciFi station using parametrization*/
    
    // extrapolate veloUT tracks
    float tx,ty,qop;

    for(int i_veloUT_track = 0; i_veloUT_track < n_veloUT_tracks_event; ++i_veloUT_track ) {
      // veloUT track variables
      const float qop = ut_tracks.qop[i_veloUT_track];
      const int i_velo_track = ut_tracks.velo_track[i_veloUT_track];
      const MiniState velo_state {velo_states, velo_event_tracks_offset + i_velo_track};
      const int ut_track_index = ut_event_tracks_offset + i_veloUT_track;
      const float ut_x  = host_ut_x[ut_track_index];
      const float ut_tx = host_ut_tx[ut_track_index];
      const float ut_z  = host_ut_z[ut_track_index];

      // SciFi IDs for matched veloUT track
      const std::vector<uint32_t> true_scifi_ids = scifi_ids_ut_tracks[i_event][i_veloUT_track];
 
      // extrapolate velo y & ty to z of UT x and tx
      // use ty from Velo state
      MiniState state_UT;
      state_UT.x = ut_x;
      state_UT.tx = ut_tx;
      state_UT.z = ut_z;
      state_UT.ty = velo_state.ty;
      state_UT.y = y_at_z(velo_state, ut_z);
      
      // extrapolate state to last UT plane (needed as input for parametrization)
      const int z_last_UT_plane = 2642.f;
      MiniState UT_state_from_velo = state_at_z(velo_state, z_last_UT_plane);
      MiniState UT_state = state_at_z(state_UT, z_last_UT_plane);
      
#ifdef WITH_ROOT
      UT_x = UT_state.x;
      UT_y = UT_state.y;
      UT_z = UT_state.z;
      UT_tx = UT_state.tx;
      UT_ty = UT_state.ty;
      ut_qop = host_ut_qop[ut_track_index];
      t_ut_tracks->Fill();
#endif

      // propagation extrap() de ZINI a ZFIN
      float eloss = 2.f/3e5;
      
      int ret = extrap(
        UT_state_from_velo.x, UT_state_from_velo.y,
        UT_state_from_velo.tx, UT_state_from_velo.ty,
        qop, scifi_params,
        xf, yf, txf, tyf, der_xf_qop);
      
      if ( true_scifi_ids.size() == 0 )
        continue;

      if ( ret ) {
      // find x hit(s) in layer 0 that 
      // were truth matched to the veloUT track
      // first layer = zone 0 (y < 0) + zone 1 (y > 0)
        int x_zone_offset_begin_0 = scifi_hit_count.zone_offset(0);
        int n_hits_0 = scifi_hit_count.zone_number_of_hits(0);
        int x_zone_offset_begin_1 = scifi_hit_count.zone_offset(1);
        int n_hits_1 = scifi_hit_count.zone_number_of_hits(1);
        int n_hits, x_zone_offset_begin;
        if ( yf < 0 ) {
          n_hits = n_hits_0;
          x_zone_offset_begin = x_zone_offset_begin_0;
        } else {
          n_hits = n_hits_1;
          x_zone_offset_begin = x_zone_offset_begin_1;
        }
      
        bool match = false;
        float true_x_0;
        for ( const auto true_id : true_scifi_ids ) {
          for ( int i_hit = 0; i_hit < n_hits; ++i_hit ) { 
            const int hit_index = x_zone_offset_begin + i_hit;
            const uint32_t lhcbid = scifi_hits.LHCbID(hit_index);
            if ( true_id == lhcbid ) {
              res_x_0 = xf - scifi_hits.x0[hit_index];
              true_x_0 = scifi_hits.x0[hit_index];
              x_extrap = xf + der_xf_qop * (xf - true_x_0);
              match = true;
              break;
            }
          }
          if ( match ) break;
        }
        
        int ret_qop = update_qop_estimate(
          UT_state_from_velo, qop,
          true_x_0, scifi_params, 
          xf, yf, txf, tyf, der_xf_qop, qop_update);
        
        if ( ret_qop ) {
          // check momentum resolution
          float p_true = p_events[i_event][i_veloUT_track];
          p_diff_after_update = p_true - std::abs( 1.f/qop_update );
          p_diff_before_update = p_true - std::abs( 1.f/qop );
        }
        
        if (!match) {
          //debug_cout << "Size of matched ids = " << true_scifi_ids.size() << std::endl;
          res_x_0 = -10000;
        }
      
        // check combinatorics within search window in layer 0
        const float max_dx = 100;
        std::vector<int> candidates_x_0;
        for ( int i_hit = 0; i_hit < n_hits; ++i_hit ) { 
          const int hit_index = x_zone_offset_begin + i_hit;
          const float x = scifi_hits.x0[hit_index];
          if ( fabsf(x-xf) < max_dx )
            candidates_x_0.push_back(hit_index);
        }
        n_hits_in_window_0 = candidates_x_0.size();

        // find x hit(s) in layer 3  that 
        // were truth matched to the veloUT track
        // first layer = zone 6 + zone 7
        // int x_zone_offset_begin_6 = scifi_hit_count.zone_offset(6);
        // int n_hits_6 = scifi_hit_count.zone_number_of_hits(6);
        // int x_zone_offset_begin_7 = scifi_hit_count.zone_offset(7);
        // int n_hits_7 = scifi_hit_count.zone_number_of_hits(7);
        // int n_hits_3, x_zone_offset_begin_3;
        // if ( yf < 0 ) {
        //   n_hits_3 = n_hits_6;
        //   x_zone_offset_begin_3 = x_zone_offset_begin_6;
        // } else {
        //   n_hits_3 = n_hits_7;
        //   x_zone_offset_begin_3 = x_zone_offset_begin_7;
        // }
        
        //bool match_3 = false;
         // float true_x_3;
        // for ( const auto true_id : true_scifi_ids ) {
        //   for ( int i_hit = 0; i_hit < n_hits_3; ++i_hit ) { 
        //     const int hit_index = x_zone_offset_begin_3 + i_hit;
        //     const uint32_t lhcbid = scifi_hits.LHCbID(hit_index);
        //     if ( true_id == lhcbid ) {
        //       res_x_3 = xf - scifi_hits.x0[hit_index];
        //       true_x_3 = scifi_hits.x0[hit_index];
        //       match_3 = true;
        //       break;
        //     }
        //   }
        //   if ( match_3 ) break;
        // }
        // if (!match_3) {
        //   res_x_3 = -10000;
        // }
        // // check combinatorics within search window in layer 3
        // const float max_dx_3 = 200;
        // std::vector<int> candidates_x_3;
        // for ( int i_hit = 0; i_hit < n_hits_3; ++i_hit ) { 
        //   const int hit_index = x_zone_offset_begin_3 + i_hit;
        //   const float x = scifi_hits.x0[hit_index];
        //   if ( fabsf(x-xf) < max_dx_3 )
        //     candidates_x_3.push_back(hit_index);
        // }
        // n_hits_in_window_3 = candidates_x_3.size();
        // dx = true_x_3-true_x_0;

        // // check combinations of x hits from layer 0 and 3
        // n_x_combinations = 0;
        // for ( const auto index_0 : candidates_x_0 ) {
        //   const float x_0 =  scifi_hits.x0[index_0];
        //   for ( const auto index_3 : candidates_x_3 ) {
        //     const float x_3 =  scifi_hits.x0[index_3];
        //     if ( fabsf(x_0 - x_3 ) < 200 ) {
        //       n_x_combinations++;
        //     }
        //   }
        // }

        t_extrap->Fill();
        
      } // extrapolation worked
    } // loop over veloUT tracks
    
    
    *n_forward_tracks = 0;
    
#ifdef WITH_ROOT
    // store qop in tree
    for ( int i_track = 0; i_track < *n_forward_tracks; ++i_track ) {
      qop = scifi_tracks_event[i_track].qop;
      state_x  = scifi_tracks_event[i_track].state.x;
      state_y  = scifi_tracks_event[i_track].state.y;
      state_z  = scifi_tracks_event[i_track].state.z;
      state_tx = scifi_tracks_event[i_track].state.tx;
      state_ty = scifi_tracks_event[i_track].state.ty;
      t_Forward_tracks->Fill();
    }
    n_tracks = n_forward_tracks[i_event];
    t_statistics->Fill();
#endif 
  }
  
#ifdef WITH_ROOT
  f->Write();
  f->Close();
#endif

  return 0;
}
