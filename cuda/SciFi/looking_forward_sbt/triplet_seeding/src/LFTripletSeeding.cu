#include "LFTripletSeeding.cuh"
#include "LFTripletSeedingImpl.cuh"
#include "LookingForwardTools.cuh"

__global__ void lf_triplet_seeding(
  uint32_t* dev_scifi_hits,
  const uint32_t* dev_scifi_hit_count,
  const int* dev_atomics_ut,
  const float* dev_ut_qop,
  const MiniState* dev_ut_states,
  const char* dev_scifi_geometry,
  const float* dev_inv_clus_res,
  const uint* dev_scifi_lf_number_of_candidates,
  const short* dev_scifi_lf_candidates,
  const LookingForward::Constants* dev_looking_forward_constants,
  float* dev_scifi_lf_triplet_best_chi2,
  int8_t* dev_scifi_lf_triplet_best_h0h2)
{
  const uint number_of_events = gridDim.x;
  const uint event_number = blockIdx.x;

  // UT consolidated tracks
  const auto ut_event_tracks_offset = dev_atomics_ut[number_of_events + event_number];
  const auto ut_event_number_of_tracks = dev_atomics_ut[number_of_events + event_number + 1] - ut_event_tracks_offset;

  // SciFi hits
  const uint total_number_of_hits = dev_scifi_hit_count[number_of_events * SciFi::Constants::n_mat_groups_and_mats];
  const SciFi::HitCount scifi_hit_count {(uint32_t*) dev_scifi_hit_count, event_number};
  const SciFi::SciFiGeometry scifi_geometry {dev_scifi_geometry};
  const SciFi::Hits scifi_hits {dev_scifi_hits, total_number_of_hits, &scifi_geometry, dev_inv_clus_res};
  const auto event_offset = scifi_hit_count.event_offset();

  for (uint16_t i = blockIdx.y; i < ut_event_number_of_tracks; i += gridDim.y) {
    const auto current_ut_track_index = ut_event_tracks_offset + i;
    const auto scifi_lf_candidates = dev_scifi_lf_candidates + current_ut_track_index *
                                                                 LookingForward::number_of_x_layers *
                                                                 LookingForward::maximum_number_of_candidates;

    const auto qop = dev_ut_qop[current_ut_track_index];

    for (uint8_t relative_first_layer = 0; relative_first_layer < 4; ++relative_first_layer) {
      const uint8_t candidate_h0_size = dev_scifi_lf_number_of_candidates
        [current_ut_track_index * LookingForward::number_of_x_layers + relative_first_layer];
      const uint8_t candidate_h1_size = dev_scifi_lf_number_of_candidates
        [current_ut_track_index * LookingForward::number_of_x_layers + relative_first_layer + 1];
      const uint8_t candidate_h2_size = dev_scifi_lf_number_of_candidates
        [current_ut_track_index * LookingForward::number_of_x_layers + relative_first_layer + 2];

      const auto z0 = dev_looking_forward_constants->Zone_zPos_xlayers[relative_first_layer];
      const auto z1 = dev_looking_forward_constants->Zone_zPos_xlayers[relative_first_layer + 1];
      const auto z2 = dev_looking_forward_constants->Zone_zPos_xlayers[relative_first_layer + 2];

      lf_triplet_seeding_impl(
        scifi_hits.x0 + event_offset,
        candidate_h0_size,
        candidate_h1_size,
        candidate_h2_size,
        relative_first_layer,
        LookingForward::chi2_max_triplet_single,
        dev_scifi_lf_triplet_best_chi2 +
          (current_ut_track_index * 4 + relative_first_layer) * LookingForward::maximum_number_of_candidates,
        dev_scifi_lf_triplet_best_h0h2 +
          (current_ut_track_index * 4 + relative_first_layer) * 2 * LookingForward::maximum_number_of_candidates,
        scifi_lf_candidates,
        z1 - z0,
        z2 - z0,
        qop);
    }
  }
}
