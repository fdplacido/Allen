#include "ConsolidateSciFi.cuh"

__global__ void consolidate_scifi_tracks(
  uint* dev_scifi_hits,
  uint* dev_scifi_hit_count,
  char* dev_scifi_consolidated_hits,
  int* dev_atomics_scifi,
  uint* dev_scifi_track_hit_number,
  float* dev_scifi_qop,
  MiniState* dev_scifi_states,
  uint* dev_scifi_track_ut_indices,
  int* dev_atomics_ut,
  SciFi::TrackHits* dev_scifi_tracks,
  const uint* dev_scifi_selected_track_indices,
  const float* dev_scifi_lf_track_params,
  const char* dev_scifi_geometry,
  const float* dev_inv_clus_res)
{

  const uint number_of_events = gridDim.x;
  const uint event_number = blockIdx.x;

  const uint ut_event_tracks_offset = dev_atomics_ut[number_of_events + event_number];

  const SciFi::TrackHits* event_scifi_tracks = dev_scifi_tracks + ut_event_tracks_offset * SciFi::Constants::max_SciFi_tracks_per_UT_track;
  const uint* event_scifi_selected_track_indices = dev_scifi_selected_track_indices + ut_event_tracks_offset * SciFi::Constants::max_SciFi_tracks_per_UT_track;

  const uint total_number_of_scifi_hits =
    dev_scifi_hit_count[number_of_events * SciFi::Constants::n_mat_groups_and_mats];
  const SciFi::SciFiGeometry scifi_geometry {dev_scifi_geometry};
  SciFi::Hits scifi_hits(dev_scifi_hits, total_number_of_scifi_hits, &scifi_geometry, dev_inv_clus_res);
  const SciFi::HitCount scifi_hit_count {dev_scifi_hit_count, event_number};

  // Create consolidated SoAs.
  SciFi::Consolidated::Tracks scifi_tracks {(uint*) dev_atomics_scifi,
                                            dev_scifi_track_hit_number,
                                            dev_scifi_qop,
                                            dev_scifi_states,
                                            dev_scifi_track_ut_indices,
                                            event_number,
                                            number_of_events};
  const uint number_of_tracks_event = scifi_tracks.number_of_tracks(event_number);
  const uint event_offset = scifi_hit_count.event_offset();

  // Loop over tracks.
  for (uint i = threadIdx.x; i < number_of_tracks_event; i += blockDim.x) {
    scifi_tracks.ut_track[i] = event_scifi_tracks[i].ut_track_index;
    scifi_tracks.qop[i] = event_scifi_tracks[i].qop;
    const uint original_track = event_scifi_selected_track_indices[i];
    const float* trackParams = dev_scifi_lf_track_params + ut_event_tracks_offset * LookingForward::maximum_number_of_candidates_per_ut_track_after_x_filter * SciFi::Tracking::nTrackParams + original_track * SciFi::Tracking::nTrackParams;
    MiniState scifi_state(
        LookingForward::x_at_end_scifi(trackParams),
        LookingForward::y_at_end_scifi(trackParams),
        SciFi::Constants::ZEndT,
        LookingForward::tx_at_end_scifi(trackParams),
        LookingForward::ty_at_end_scifi(trackParams));
    scifi_tracks.states[i] = scifi_state;
    SciFi::Consolidated::Hits consolidated_hits =
      scifi_tracks.get_hits(dev_scifi_consolidated_hits, i, &scifi_geometry, dev_inv_clus_res);
    SciFi::TrackHits track = event_scifi_tracks[i];

    // Lambda for populating arrays.
    auto populate = [&track](uint32_t* __restrict__ a, uint32_t* __restrict__ b) {
      for (int i = 0; i < track.hitsNum; i++) {
        const auto hit_index = track.hits[i];
        a[i] = b[hit_index];
      }
    };

    // Populate the consolidated hits.
    populate((uint32_t*) consolidated_hits.x0, (uint32_t*) scifi_hits.x0 + event_offset);
    populate((uint32_t*) consolidated_hits.z0, (uint32_t*) scifi_hits.z0 + event_offset);
    populate((uint32_t*) consolidated_hits.m_endPointY, (uint32_t*) scifi_hits.m_endPointY + event_offset);
    populate((uint32_t*) consolidated_hits.channel, (uint32_t*) scifi_hits.channel + event_offset);
    populate(
      (uint32_t*) consolidated_hits.assembled_datatype, (uint32_t*) scifi_hits.assembled_datatype + event_offset);
  }
}
