#pragma once

#include "UTEventModel.cuh"
#include "UTDefinitions.cuh"
#include "Handler.cuh"
#include "ArgumentsCommon.cuh"
#include "ArgumentsUT.cuh"

__global__ void ut_find_permutation(
  uint32_t* dev_ut_hits,
  uint32_t* dev_ut_hit_offsets,
  uint* dev_hit_permutations,
  const uint* dev_unique_x_sector_layer_offsets);

ALGORITHM(
  ut_find_permutation,
  ut_find_permutation_t,
  ARGUMENTS(dev_ut_hits, dev_ut_hit_offsets, dev_ut_hit_permutations))
