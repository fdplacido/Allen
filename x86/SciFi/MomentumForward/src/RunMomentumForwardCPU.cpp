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

  // initialize parameters
  //char name_coef_T1[200] = "../input/test_UT_T1.tab";
  char name_coef_T1[200] = "../input/UT_T1_tilt_new.tab";
  debug_cout << "Reading coefs for extrapolation to T1: " << name_coef_T1 << std::endl;
  SciFi::Parameters scifi_params_T1 = SciFi::Parameters(name_coef_T1);
  
  char name_coef_T3[200] = "../input/test_UT_T3.tab";
  debug_cout << "Reading coefs for extrapolation to T3: " << name_coef_T3 << std::endl;
  SciFi::Parameters scifi_params_T3 = SciFi::Parameters(name_coef_T3);

#ifdef WITH_ROOT
  // Histograms only for checking and debugging
  TFile *f = new TFile("../output/scifi.root", "RECREATE");
  TTree *t_Forward_tracks = new TTree("Forward_tracks", "Forward_tracks");
  TTree *t_statistics = new TTree("statistics", "statistics");
  TTree *t_scifi_hits = new TTree("scifi_hits","scifi_hits");
  TTree *t_extrap_T1 = new TTree("extrap_T1","extrap_T1");
  TTree *t_extrap_T3 = new TTree("extrap_T3","extrap_T3");
  TTree *t_ut_tracks = new TTree("ut_tracks","ut_tracks");
  TTree* t_other_x_layer_t1 = new TTree("other_x_layer_T1","other_x_layer_T1");
  TTree* t_other_x_layer_t3 = new TTree("other_x_layer_T3","other_x_layer_T3");
  
  uint planeCode, LHCbID;
  float x0, z0, w, dxdy, dzdy, yMin, yMax;
  float qop;
  int n_tracks;
  float state_x, state_y, state_z, state_tx, state_ty;
  float xf_t1, yf_t1, txf_t1, tyf_t1, der_xf_qop_t1, qop_update_t1;
  float res_x_0_t1, res_x_3_t1, dx_t1, x_extrap_t1, true_x_t1, res_x_other_t1, true_z_t1;
  float UT_x, UT_y, UT_z, UT_tx, UT_ty, ut_qop;
  float velo_x_extrap, velo_tx;
  int n_hits_in_window_0_t1 = 0, n_hits_in_window_0_t1_true_p = 0, n_hits_in_window_3_t1 = 0;
  int n_hits_in_zone_t1 = 0, n_hits_in_window_other_t1 = 0;
  float p_diff_before_update_t1, p_diff_after_update_t1, p_diff_before_after_t1, p_resolution_after_update_t1;
  float qop_diff_before_update_t1, qop_diff_after_update_t1, qop_diff_before_after_t1, qop_resolution_after_update_t1;
  float tx_x_hits_t1;

  float xf_t3, yf_t3, txf_t3, tyf_t3, der_xf_qop_t3, qop_update_t3;
  float res_x_0_t3, res_x_3_t3, dx_t3, x_extrap_t3, true_x_t3, res_x_other_t3;;
  int n_hits_in_window_0_t3 = 0, n_hits_in_window_0_t3_true_p = 0, n_hits_in_window_3_t3 = 0;
  int n_hits_in_zone_t3 = 0, n_hits_in_window_other_t3 = 0;
  float p_diff_before_update_t3, p_diff_after_update_t3, p_resolution_after_update_t3;
  float qop_diff_before_update_t3, qop_diff_after_update_t3, qop_diff_before_after_t3, qop_resolution_after_update_t3;
  float tx_x_hits_t3;

  bool t1_extrap_worked, t3_extrap_worked, isLong;
  float p_true;
  bool match_t1, match_t3, match_t1_other, match_t3_other;

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

  t_extrap_T1->Branch("xf", &xf_t1);
  t_extrap_T1->Branch("yf", &yf_t1);
  t_extrap_T1->Branch("txf", &txf_t1);
  t_extrap_T1->Branch("tyf", &tyf_t1);
  t_extrap_T1->Branch("der_xf_qop", &der_xf_qop_t1);
  t_extrap_T1->Branch("res_x_0", &res_x_0_t1);
  t_extrap_T1->Branch("res_x_other", &res_x_other_t1);
  t_extrap_T1->Branch("true_x", &true_x_t1);
  t_extrap_T1->Branch("true_z", &true_z_t1);
  t_extrap_T1->Branch("x_extrap", &x_extrap_t1);
  t_extrap_T1->Branch("dx", &dx_t1);
  t_extrap_T1->Branch("qop_update", &qop_update_t1);
  t_extrap_T1->Branch("p_diff_before_update", &p_diff_before_update_t1);
  t_extrap_T1->Branch("p_diff_after_update", &p_diff_after_update_t1);
  t_extrap_T1->Branch("p_resolution_after_update", &p_resolution_after_update_t1);
  t_extrap_T1->Branch("qop_diff_before_update", &qop_diff_before_update_t1);
  t_extrap_T1->Branch("qop_diff_after_update", &qop_diff_after_update_t1);
  t_extrap_T1->Branch("qop_resolution_after_update", &qop_resolution_after_update_t1);
  t_extrap_T1->Branch("n_hits_in_window_0", &n_hits_in_window_0_t1);
  t_extrap_T1->Branch("n_hits_in_zone", &n_hits_in_zone_t1);
  t_extrap_T1->Branch("isLong", &isLong);
  t_extrap_T1->Branch("p_true", &p_true);
  t_extrap_T1->Branch("match", &match_t1);
  t_extrap_T1->Branch("match_other", &match_t1_other);
  t_extrap_T1->Branch("ut_qop", &ut_qop);
  t_extrap_T1->Branch("tx_x_hits", &tx_x_hits_t1);

  t_extrap_T3->Branch("xf", &xf_t3);
  t_extrap_T3->Branch("yf", &yf_t3);
  t_extrap_T3->Branch("txf", &txf_t3);
  t_extrap_T3->Branch("tyf", &tyf_t3);
  t_extrap_T3->Branch("der_xf_qop", &der_xf_qop_t3);
  t_extrap_T3->Branch("res_x_0", &res_x_0_t3);
  t_extrap_T3->Branch("res_x_other", &res_x_other_t3);
  t_extrap_T3->Branch("true_x", &true_x_t3);
  t_extrap_T3->Branch("x_extrap", &x_extrap_t3);
  t_extrap_T3->Branch("dx", &dx_t3);
  t_extrap_T3->Branch("qop_update", &qop_update_t3);
  t_extrap_T3->Branch("p_diff_before_update", &p_diff_before_update_t3);
  t_extrap_T3->Branch("p_diff_after_update", &p_diff_after_update_t3);
  t_extrap_T3->Branch("p_resolution_after_update", &p_resolution_after_update_t3);
  t_extrap_T3->Branch("qop_diff_before_update", &qop_diff_before_update_t3);
  t_extrap_T3->Branch("qop_diff_after_update", &qop_diff_after_update_t3);
  t_extrap_T3->Branch("qop_resolution_after_update", &qop_resolution_after_update_t3);
  t_extrap_T3->Branch("n_hits_in_window_0", &n_hits_in_window_0_t3);
  t_extrap_T3->Branch("n_hits_in_zone", &n_hits_in_zone_t3);
  t_extrap_T3->Branch("isLong", &isLong);
  t_extrap_T3->Branch("p_true", &p_true);
  t_extrap_T3->Branch("ut_qop", &ut_qop);
  t_extrap_T3->Branch("tx_x_hits", &tx_x_hits_t3);
  t_extrap_T3->Branch("match", &match_t3);
  t_extrap_T3->Branch("match_other", &match_t3_other);

  t_ut_tracks->Branch("ut_x", &UT_x);
  t_ut_tracks->Branch("ut_y", &UT_y);
  t_ut_tracks->Branch("ut_z", &UT_z);
  t_ut_tracks->Branch("ut_tx", &UT_tx);
  t_ut_tracks->Branch("ut_ty", &UT_ty);
  t_ut_tracks->Branch("velo_x_extrap", &velo_x_extrap);
  t_ut_tracks->Branch("velo_tx", &velo_tx);
  t_ut_tracks->Branch("ut_qop", &ut_qop);
  t_ut_tracks->Branch("t1_extrap_worked", &t1_extrap_worked);
  t_ut_tracks->Branch("t3_extrap_worked", &t3_extrap_worked);
  t_ut_tracks->Branch("isLong", &isLong);
  t_ut_tracks->Branch("p_true", &p_true);

  t_other_x_layer_t1->Branch("n_hits_in_window_other", &n_hits_in_window_other_t1);

  t_other_x_layer_t3->Branch("n_hits_in_window_other", &n_hits_in_window_other_t3);
#endif

  int n_veloUT_tracks = 0;
  int n_extrap_T1 = 0;
  int n_extrap_T3 = 0;
  
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
    n_veloUT_tracks += n_veloUT_tracks_event;
  
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
      const int z_last_UT_plane = 2642.;
      MiniState UT_state_from_velo = state_at_z(velo_state, z_last_UT_plane);
      MiniState UT_state = state_at_z(state_UT, z_last_UT_plane);
      
      //debug_cout << "ut_z = " << ut_z << std::endl;
      
      t1_extrap_worked = false;
      t3_extrap_worked = false;
      isLong = false;
      match_t1 = false;
      match_t3 = false;
      n_hits_in_window_0_t1 = 0;
      n_hits_in_window_0_t3 = 0;
      
      // sign of momentum -> charge of MCP
      // Caution: this is only set correctly if the veloUT track was matched to an MCP
      // set to 1e9 if not matched
      p_true = p_events[i_event][i_veloUT_track];

#ifdef WITH_ROOT
      UT_x = UT_state.x;
      UT_y = UT_state.y;
      UT_z = UT_state.z;
      UT_tx = UT_state.tx;
      UT_ty = UT_state.ty;
      ut_qop = qop;
#endif

      // propagation to first layer of T1
      int ret = extrap(
        UT_state.x, UT_state.y,
        UT_state.tx, UT_state.ty,
        qop, scifi_params_T1,
        xf_t1, yf_t1, txf_t1, tyf_t1, der_xf_qop_t1);
      
      if ( ret ) {
        n_extrap_T1++;
        t1_extrap_worked = true;
        
        // access hits in first layer of T1, layer 0 = zone 0 (y < 0) + zone 1 (y > 0)
        int x_zone_offset_begin_0 = scifi_hit_count.zone_offset(0);
        int n_hits_0 = scifi_hit_count.zone_number_of_hits(0);
        int x_zone_offset_begin_1 = scifi_hit_count.zone_offset(1);
        int n_hits_1 = scifi_hit_count.zone_number_of_hits(1);
        int n_hits, x_zone_offset_begin;
        if ( yf_t1 < 0 ) {
          n_hits = n_hits_0;
          x_zone_offset_begin = x_zone_offset_begin_0;
        } else {
          n_hits = n_hits_1;
          x_zone_offset_begin = x_zone_offset_begin_1;
        }
        // access hits in last layer of T1, layer 3 = zone 6 (y < 0) + zone 7 (y > 0)
        int x_zone_offset_begin_6 = scifi_hit_count.zone_offset(6);
        int n_hits_6 = scifi_hit_count.zone_number_of_hits(6);
        int x_zone_offset_begin_7 = scifi_hit_count.zone_offset(7);
        int n_hits_7 = scifi_hit_count.zone_number_of_hits(7);
        int n_hits_other, x_zone_offset_begin_other;
        if ( yf_t1 < 0 ) {
          n_hits_other = n_hits_6;
          x_zone_offset_begin_other = x_zone_offset_begin_6;
        } else {
          n_hits_other = n_hits_7;
          x_zone_offset_begin_other = x_zone_offset_begin_7;
        }

        // find x hit(s) in layer 0 that 
        // were truth matched to the veloUT track
        if ( true_scifi_ids.size() > 0 ) {
          isLong = true;
          match_t1 = false;
          for ( const auto true_id : true_scifi_ids ) {
            for ( int i_hit = 0; i_hit < n_hits; ++i_hit ) { 
              const int hit_index = x_zone_offset_begin + i_hit;
              const uint32_t lhcbid = scifi_hits.LHCbID(hit_index);
              if ( true_id == lhcbid ) {
                res_x_0_t1 = xf_t1 - scifi_hits.x0[hit_index];
                true_x_t1 = scifi_hits.x0[hit_index];
                true_z_t1 = scifi_hits.z0[hit_index];
                x_extrap_t1 = xf_t1 + der_xf_qop_t1 * (xf_t1 - true_x_t1);
                match_t1 = true;
                break;
              }
            }
            if ( match_t1 ) break;
          }
          
          // Update momentum estimate with true x hit position
          // caution: xf, yf etc. are changed in the qop update step
          // -> they are not the same as above when checking the resolution (res_x_0_t1)
          if ( match_t1 ) {

            int ret_qop = update_qop_estimate(
              UT_state, qop,
              true_x_t1, scifi_params_T1, 
              xf_t1, yf_t1, txf_t1, tyf_t1, der_xf_qop_t1, qop_update_t1);
            
            if ( ret_qop ) {
              // check momentum resolution
              p_diff_after_update_t1 = std::abs(p_true) - std::abs( 1.f/qop_update_t1 );
              p_diff_before_update_t1 = std::abs(p_true) - std::abs( 1.f/qop );
              p_resolution_after_update_t1 = (std::abs(p_true) - std::abs( 1.f/qop_update_t1 )) / std::abs(p_true);
              qop_diff_after_update_t1 = 1./p_true - qop_update_t1;
              qop_diff_before_update_t1 = 1./p_true - qop;
              qop_diff_before_after_t1 = qop - qop_update_t1;
              qop_resolution_after_update_t1 = (1./p_true - qop_update_t1) * p_true;
            }

            // Distance in x to correct hit in other x layer of station
            match_t1_other = false;
            float slope1, slope2;
            if ( qop < 0 ) {
              slope1 = 0.3e6;
              slope2 = 0.2e6;
            } else {
              slope1 = -0.2e6;
              slope2 = -0.3e6;
            }
            for ( const auto true_id : true_scifi_ids ) {
              for ( int i_hit = 0; i_hit < n_hits_other; ++i_hit ) { 
                const int hit_index = x_zone_offset_begin_other + i_hit;
                const uint32_t lhcbid = scifi_hits.LHCbID(hit_index);
                if ( true_id == lhcbid ) {
                  float true_x_t1_other = scifi_hits.x0[hit_index];
                  if ( fabsf(true_x_t1-true_x_t1_other) < 20.f - slope1 * qop 
                       && fabsf(true_x_t1-true_x_t1_other) > -20.f - slope2 * qop) {
                    float true_x_t1_other = scifi_hits.x0[hit_index];
                    res_x_other_t1 = true_x_t1 - true_x_t1_other;
                    match_t1_other = true;
                    tx_x_hits_t1 = (true_x_t1 - true_x_t1_other) / 210.; // dz of x-layers within one station = 210 mm
                    break;
                  }
                }
              }
              if ( match_t1_other ) break;
            }
            
          } // match in first layer of T1
        
        } // # of true SciFi IDs > 0
        
        // check combinatorics within search window in layer 0
        float slope1, slope2;
        if ( qop < 0 ) {
          slope1 = 0.3e6;
          slope2 = 0.2e6;
        } else {
          slope1 = -0.2e6;
          slope2 = -0.3e6;
        }
                
        for ( int i_hit = 0; i_hit < n_hits; ++i_hit ) { 
          const int hit_index = x_zone_offset_begin + i_hit;
          const float x = scifi_hits.x0[hit_index];
          if ( fabsf(x-xf_t1) < 20 + 1.e6 * fabsf(qop) ) {
            n_hits_in_window_0_t1++;
            // check combinatorics in other x-layer of last station
            n_hits_in_window_other_t1 = 0;
            for ( int i_hit_other = 0; i_hit_other < n_hits_other; ++i_hit_other ) { 
              const int hit_index_other = x_zone_offset_begin_other + i_hit_other;
              const float x_other = scifi_hits.x0[hit_index_other];
              if ( fabsf(x-x_other) < 20.f - slope1 * qop 
                   && fabsf(x-x_other) > -20.f - slope2 * qop) {
                n_hits_in_window_other_t1++;
              }
            }
            t_other_x_layer_t1->Fill();
          }
        }
        n_hits_in_zone_t1 = n_hits;
      
        t_extrap_T1->Fill();
        
      } // extrapolation to T1 worked
      
      
      // extrapolate to last layer of T3 (layer 11)
      ret = extrap(
        UT_state.x, UT_state.y,
        UT_state.tx, UT_state.ty,
        qop, scifi_params_T3,
        xf_t3, yf_t3, txf_t3, tyf_t3, der_xf_qop_t3);
      
      if ( ret ) {
        n_extrap_T3++;
        t3_extrap_worked = true;
        
        
        // Access hits in last layer of T3: layer 11 = zone 22 (y < 0) + zone 23 (y > 0)
        int x_zone_offset_begin_22 = scifi_hit_count.zone_offset(22);
        int n_hits_22 = scifi_hit_count.zone_number_of_hits(22);
        int x_zone_offset_begin_23 = scifi_hit_count.zone_offset(23);
        int n_hits_23 = scifi_hit_count.zone_number_of_hits(23);
        int n_hits, x_zone_offset_begin;
        if ( yf_t3 < 0 ) {
          n_hits = n_hits_22;
          x_zone_offset_begin = x_zone_offset_begin_22;
        } else {
          n_hits = n_hits_23;
          x_zone_offset_begin = x_zone_offset_begin_23;
        }
        // Access hits in first layer of T3, layer 8 = zone 16 (y < 0) + zone 17 (y > 0)
        int x_zone_offset_begin_16 = scifi_hit_count.zone_offset(16);
        int n_hits_16 = scifi_hit_count.zone_number_of_hits(16);
        int x_zone_offset_begin_17 = scifi_hit_count.zone_offset(17);
        int n_hits_17 = scifi_hit_count.zone_number_of_hits(17);
        int n_hits_other, x_zone_offset_begin_other;
        if ( yf_t3 < 0 ) {
          n_hits_other = n_hits_16;
          x_zone_offset_begin_other = x_zone_offset_begin_16;
        } else {
          n_hits_other = n_hits_17;
          x_zone_offset_begin_other = x_zone_offset_begin_17;
        }
        
        // find x hit(s) in layer 11 that 
        // were truth matched to the veloUT track
        if ( true_scifi_ids.size() > 0 ) {
          match_t3 = false;
          for ( const auto true_id : true_scifi_ids ) {
            for ( int i_hit = 0; i_hit < n_hits; ++i_hit ) { 
              const int hit_index = x_zone_offset_begin + i_hit;
              const uint32_t lhcbid = scifi_hits.LHCbID(hit_index);
              if ( true_id == lhcbid ) {
                res_x_0_t3 = xf_t3 - scifi_hits.x0[hit_index];
                true_x_t3 = scifi_hits.x0[hit_index];
                x_extrap_t3 = xf_t3 + der_xf_qop_t3 * (xf_t3 - true_x_t3);
                match_t3 = true;
                break;
              }
            }
            if ( match_t3 ) break;
          }
        
          // Update momentum estimate with true x hit position
          if ( match_t3 ) {
            int ret_qop = update_qop_estimate(
              UT_state, qop,
              true_x_t3, scifi_params_T3, 
              xf_t3, yf_t3, txf_t3, tyf_t3, der_xf_qop_t3, qop_update_t3);
            
            if ( ret_qop ) {
              // check momentum resolution
              p_diff_after_update_t3 = std::abs(p_true) - std::abs( 1.f/qop_update_t3 );
              p_diff_before_update_t3 = std::abs(p_true) - std::abs( 1.f/qop );
              p_resolution_after_update_t3 = (std::abs(p_true) - std::abs( 1.f/qop_update_t3 )) / std::abs(p_true);
              qop_diff_after_update_t3 = 1./p_true - qop_update_t3;
              qop_diff_before_update_t3 = 1./p_true - qop;
              qop_resolution_after_update_t3 = (1./p_true - qop_update_t3) * p_true;
            }
            
            // Distance in x to correct hit in other x layer of station
            match_t3_other = false;
            float slope1, slope2;
            if ( qop < 0 ) {
              slope1 = -0.2e6;
              slope2 = -0.3e6;
            } else {
              slope1 = 0.3e6;
              slope2 = 0.2e6;
            }
            for ( const auto true_id : true_scifi_ids ) {
              for ( int i_hit = 0; i_hit < n_hits_other; ++i_hit ) { 
                const int hit_index = x_zone_offset_begin_other + i_hit;
                const uint32_t lhcbid = scifi_hits.LHCbID(hit_index);
                if ( true_id == lhcbid ) {
                  float true_x_t3_other = scifi_hits.x0[hit_index];
                  if ( fabsf(true_x_t3-true_x_t3_other) < 20.f + slope1 * qop 
                       && fabsf(true_x_t3-true_x_t3_other) > -20.f + slope2 * qop) {
                    res_x_other_t3 = true_x_t3 - scifi_hits.x0[hit_index];
                    tx_x_hits_t3 = (true_x_t3 - true_x_t3_other) / 210.; // dz of x-layers within one station = 210 mm
                    match_t3_other = true;
                    break;
                  }
                }
              }
              if ( match_t3_other ) break;
            }
          }
          
        } // # of true SciFi IDs > 0
        
        // check combinatorics within search window in layer 11
        float slope1, slope2;
        if ( qop < 0 ) {
          slope1 = -0.2e6;
          slope2 = -0.3e6;
        } else {
          slope1 = 0.3e6;
          slope2 = 0.2e6;
        }
        for ( int i_hit = 0; i_hit < n_hits; ++i_hit ) { 
          const int hit_index = x_zone_offset_begin + i_hit;
          const float x = scifi_hits.x0[hit_index];
          if ( fabsf(x-xf_t3) < 40 + 1.5e6 * fabsf(qop) ) {
            n_hits_in_window_0_t3++;
            // check combinatorics in other x-layer of last station
            n_hits_in_window_other_t3 = 0;
            for ( int i_hit_other = 0; i_hit_other < n_hits_other; ++i_hit_other ) { 
              const int hit_index_other = x_zone_offset_begin_other + i_hit_other;
              const float x_other = scifi_hits.x0[hit_index_other];
              if ( fabsf(x-x_other) < 20.f + slope1 * qop 
                   && fabsf(x-x_other) > -20.f + slope2 * qop) {
                n_hits_in_window_other_t3++;
              }
            }
            t_other_x_layer_t3->Fill();
          }
        }
        n_hits_in_zone_t3 = n_hits;

        t_extrap_T3->Fill();
        
      } // extrapolation to T3 worked
#ifdef WITH_ROOT
      t_ut_tracks->Fill();
#endif
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

  info_cout << "Extrapolation to T1 worked: " << float(n_extrap_T1) / n_veloUT_tracks << std::endl;
  info_cout << "Extrapolation to T3 worked: " << float(n_extrap_T3) / n_veloUT_tracks << std::endl;

  return 0;
}
