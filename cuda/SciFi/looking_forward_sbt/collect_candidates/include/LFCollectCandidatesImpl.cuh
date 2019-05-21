#pragma once

#include "SciFiDefinitions.cuh"
#include "SciFiEventModel.cuh"
#include "UTDefinitions.cuh"
#include "LookingForwardConstants.cuh"

__device__ void lf_collect_candidates_impl(
  const SciFi::Hits& scifi_hits,
  const int* initial_windows,
  const float qop,
  const int number_of_tracks,
  uint* number_of_candidates,
  short* candidates,
  const int event_offset,
  const int event_number_of_hits);

__device__ void lf_collect_candidates_p_impl(
  const SciFi::Hits& scifi_hits,
  const int* initial_windows,
  const float qop,
  const int number_of_tracks,
  uint* number_of_candidates,
  short* candidates,
  const int event_offset,
  const int event_number_of_hits);