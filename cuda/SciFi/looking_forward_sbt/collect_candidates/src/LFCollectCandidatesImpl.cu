#include "LFCollectCandidatesImpl.cuh"
#include "BinarySearchTools.cuh"

__device__ void lf_collect_candidates_impl(
  const SciFi::Hits& scifi_hits,
  const int* initial_windows,
  const float qop,
  const int number_of_tracks,
  uint* number_of_candidates,
  short* candidates,
  const int event_offset)
{
  //float* initial_windows_f = (float*) &initial_windows[0];
  for (int i = threadIdx.y; i < LookingForward::number_of_x_layers; i += blockDim.y) {
    const auto window_start = initial_windows[i * 8 * number_of_tracks];
    const auto window_size = initial_windows[(i * 8 + 1) * number_of_tracks];
    const auto search_window_start = initial_windows[(i * 8 + 2) * number_of_tracks];
    const auto search_window_size = initial_windows[(i * 8 + 3) * number_of_tracks];
    // const auto param_uv_0 = initial_windows_f[(i * 8 + 4) * number_of_tracks];
    // const auto param_uv_1 = initial_windows_f[(i * 8 + 5) * number_of_tracks];
    // const auto param_uv_2 = initial_windows_f[(i * 8 + 6) * number_of_tracks];
    // const auto param_uv_3 = initial_windows_f[(i * 8 + 7) * number_of_tracks];
    uint8_t candidate_counter = 0;

    for (int j = 0; j < window_size && candidate_counter < LookingForward::maximum_number_of_candidates; ++j) {
      const auto hit_index = window_start + j;
      const auto xHit = scifi_hits.x0[hit_index];

      // const auto xPredUv = param_uv_0 + xHit * param_uv_1;
      // const auto maxDx = param_uv_2 + std::abs(xHit - param_uv_3) * SciFi::Tracking::tolYSlopeCollectX;
      // const auto xMinUV = xPredUv - maxDx;
      // const auto xMaxUV = xPredUv + maxDx;
      
      const float maxDx = 25 + 6e5 * std::abs(qop); 
      const float xMinUV = xHit - maxDx;
      const float xMaxUV = xHit + maxDx;

      if (binary_search_match_stereo_hit(
            scifi_hits,
            search_window_start,
            search_window_size,
            xMinUV,
            xMaxUV)) {
        candidates[i * LookingForward::maximum_number_of_candidates + candidate_counter++] =
          (short) (hit_index - event_offset);
      }
    }

    number_of_candidates[i] = candidate_counter;
  }
}
