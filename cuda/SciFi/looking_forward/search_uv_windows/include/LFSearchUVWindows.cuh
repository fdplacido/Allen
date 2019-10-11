#pragma once

#include "LookingForwardConstants.cuh"
#include "LookingForwardTools.cuh"
#include "SciFiEventModel.cuh"
#include "Handler.cuh"
#include "ArgumentsUT.cuh"
#include "ArgumentsSciFi.cuh"

__global__ void lf_search_uv_windows(
  const uint32_t* dev_scifi_hits,
  const uint32_t* dev_scifi_hit_count,
  const uint* dev_atomics_ut,
  const SciFi::TrackHits* dev_scifi_tracks,
  const uint* dev_atomics_scifi,
  const char* dev_scifi_geometry,
  const LookingForward::Constants* dev_looking_forward_constants,
  const float* dev_inv_clus_res,
  const MiniState* dev_ut_states,
  short* dev_scifi_lf_uv_windows,
  const int* dev_scifi_lf_initial_windows);

ALGORITHM(
  lf_search_uv_windows,
  lf_search_uv_windows_t,
  ARGUMENTS(
    dev_scifi_hits,
    dev_scifi_hit_count,
    dev_atomics_ut,
    dev_scifi_lf_x_filtered_tracks,
    dev_scifi_lf_x_filtered_atomics,
    dev_scifi_lf_number_of_candidates,
    dev_scifi_lf_candidates,
    dev_ut_states,
    dev_scifi_lf_uv_windows,
    dev_scifi_lf_initial_windows))