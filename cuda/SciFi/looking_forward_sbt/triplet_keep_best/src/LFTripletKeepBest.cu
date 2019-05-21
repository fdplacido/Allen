#include "LFTripletKeepBest.cuh"

__global__ void lf_triplet_keep_best(
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
  SciFi::TrackHits* dev_scifi_tracks,
  int* dev_atomics_scifi,
  const float* dev_scifi_lf_triplet_best_chi2,
  const int8_t* dev_scifi_lf_triplet_best_h0h2)
{
  // Keep best for each h1 hit
  __shared__ float best_chi2[4 * LookingForward::maximum_number_of_candidates];
  __shared__ int8_t best_triplets[LookingForward::maximum_number_of_candidates_per_ut_track];

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
    const auto scifi_lf_candidates = dev_scifi_lf_candidates + current_ut_track_index * LookingForward::number_of_x_layers *
                                                           LookingForward::maximum_number_of_candidates;

    // Initialize shared memory buffers
    __syncthreads();

    // Initialize the best_ shared memory buffers
    for (uint16_t relative_first_layer = 0; relative_first_layer < 4; ++relative_first_layer) {
      for (uint16_t j = threadIdx.x; j < LookingForward::maximum_number_of_candidates; j += blockDim.x) {
        best_chi2[relative_first_layer * LookingForward::maximum_number_of_candidates + j] =
          dev_scifi_lf_triplet_best_chi2[(current_ut_track_index * 4 + relative_first_layer) * LookingForward::maximum_number_of_candidates + j];
      }
    }

    // Initialize best_triplets to -1
    for (uint16_t j = threadIdx.x; j < LookingForward::maximum_number_of_candidates_per_ut_track; j += blockDim.x) {
      best_triplets[j] = -1;
    }

    __syncthreads();

    // Now, we have the best candidates populated in best_chi2 and best_h0h2
    // Sort the candidates (insertion sort) into best_triplets
    for (uint16_t j = threadIdx.x; j < 4 * LookingForward::maximum_number_of_candidates; j += blockDim.x) {
      const float chi2 = best_chi2[j];
      if (chi2 < LookingForward::chi2_mean_triplet_single + 2.5f * LookingForward::chi2_stddev_triplet_single) {
        int16_t insert_position = 0;
        for (uint16_t k = 0; k < 4 * LookingForward::maximum_number_of_candidates; ++k) {
          const float other_chi2 = best_chi2[k];
          if (chi2 > other_chi2 || (chi2 == other_chi2 && j < k)) {
            ++insert_position;
          }
        }
        if (insert_position < LookingForward::maximum_number_of_candidates_per_ut_track) {
          best_triplets[insert_position] = j;
        }
      }
    }

    __syncthreads();

    // Save best triplet canidates as TrackHits candidates for further extrapolation
    for (uint16_t j = threadIdx.x; j < LookingForward::maximum_number_of_candidates_per_ut_track; j += blockDim.x) {
      const auto k = best_triplets[j];
      if (k != -1) {
        const auto relative_middle_layer = 1 + (k / LookingForward::maximum_number_of_candidates);
        const auto element = k % LookingForward::maximum_number_of_candidates;
        const auto h0_element = (current_ut_track_index * 4 + relative_middle_layer - 1) * 2 * LookingForward::maximum_number_of_candidates + element;
        const auto h2_element = (current_ut_track_index * 4 + relative_middle_layer - 1) * 2 * LookingForward::maximum_number_of_candidates + LookingForward::maximum_number_of_candidates + element;

        // Create triplet candidate with all information we have
        const int current_insert_index = atomicAdd(dev_atomics_scifi + current_ut_track_index, 1);
        assert(current_insert_index < LookingForward::maximum_number_of_candidates_per_ut_track);
        const uint16_t h0 = (uint16_t) scifi_lf_candidates[(relative_middle_layer - 1) * LookingForward::maximum_number_of_candidates + dev_scifi_lf_triplet_best_h0h2[h0_element]];
        const uint16_t h1 = (uint16_t) scifi_lf_candidates[relative_middle_layer * LookingForward::maximum_number_of_candidates + element];
        const uint16_t h2 = (uint16_t) scifi_lf_candidates[(relative_middle_layer + 1) * LookingForward::maximum_number_of_candidates +
          dev_scifi_lf_triplet_best_h0h2[h2_element]];

        const float x0 = scifi_hits.x0[event_offset + h0];
        const float x1 = scifi_hits.x0[event_offset + h1];
        const auto z0 = dev_looking_forward_constants->Zone_zPos_xlayers[relative_middle_layer - 1];
        const auto z1 = dev_looking_forward_constants->Zone_zPos_xlayers[relative_middle_layer];

        dev_scifi_tracks[current_ut_track_index * LookingForward::maximum_number_of_candidates_per_ut_track + current_insert_index] =
          SciFi::TrackHits {h0,
                            h1,
                            h2,
                            (uint16_t) (relative_middle_layer - 1),
                            (uint16_t) relative_middle_layer,
                            (uint16_t) (relative_middle_layer + 1),
                            best_chi2[k],
                            LookingForward::qop_update(
                              dev_ut_states[current_ut_track_index].tx,
                              x0,
                              z0,
                              x1,
                              z1,
                              dev_looking_forward_constants->ds_p_param_layer_inv[(relative_middle_layer - 1)]),
                            i};
      }
    }
  }

}