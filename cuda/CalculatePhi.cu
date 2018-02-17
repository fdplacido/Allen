#include "SearchByTriplet.cuh"
#include "math_constants.h"

/**
 * @brief Calculate a single hit phi in odd sensor
 */
__device__ float hit_phi_odd(
  const float x,
  const float y
) {
  return atan2(y, x);
}

/**
 * @brief Calculate a single hit phi in even sensor
 */
__device__ float hit_phi_even(
  const float x,
  const float y
) {
  const auto phi = atan2(y, x);
  const auto less_than_zero = phi < 0.f;
  return phi + less_than_zero*2*CUDART_PI_F;
}

/**
 * @brief Calculates phi for each hit
 */
__device__ void calculatePhi(
  const unsigned int* module_hitStarts,
  const unsigned int* module_hitNums,
  const float* hit_Xs,
  const float* hit_Ys,
  float* hit_Phis,
  unsigned short* hit_permutations
) {
  __shared__ float shared_hit_phis [256];

  // Even modules
  for (unsigned int module=0; module<52; module += 2) {
    const auto hit_start = module_hitStarts[module];
    const auto hit_num = module_hitNums[module];

    // Calculate phis
    for (unsigned int i=0; i<(hit_num + blockDim.x - 1) / blockDim.x; ++i) {
      const auto hit_rel_id = i*blockDim.x + threadIdx.x;
      if (hit_rel_id < hit_num) {
        const auto hit_index = hit_start + hit_rel_id;
        shared_hit_phis[hit_rel_id] = hit_phi_even(hit_Xs[hit_index], hit_Ys[hit_index]);
      }
    }

    // shared_hit_phis
    __syncthreads();

    // Find the permutations given the phis in shared_hit_phis
    for (unsigned int i=0; i<(hit_num + blockDim.x - 1) / blockDim.x; ++i) {
      const auto hit_rel_id = i*blockDim.x + threadIdx.x;
      if (hit_rel_id < hit_num) {
        const auto hit_index = hit_start + hit_rel_id;
        const auto phi = shared_hit_phis[hit_rel_id];
        
        // Find out local position
        unsigned int position = 0;
        for (unsigned int j=0; j<hit_num; ++j) {
          const auto other_phi = shared_hit_phis[j];
          // Stable sorting
          position += phi>other_phi || (phi==other_phi && hit_rel_id>j);
        }
        ASSERT(position < MAX_NUMHITS_IN_MODULE)

        // Store it in hit permutations and in hit_Phis, already ordered
        const auto global_position = hit_start + position;
        hit_permutations[global_position] = hit_index;
        hit_Phis[global_position] = phi;
      }
    }

    // shared_hit_phis
    __syncthreads();
  }

  // Odd modules
  for (unsigned int module=1; module<52; module += 2) {
    const auto hit_start = module_hitStarts[module];
    const auto hit_num = module_hitNums[module];

    // Calculate phis
    for (unsigned int i=0; i<(hit_num + blockDim.x - 1) / blockDim.x; ++i) {
      const auto hit_rel_id = i*blockDim.x + threadIdx.x;
      if (hit_rel_id < hit_num) {
        const auto hit_index = hit_start + hit_rel_id;
        shared_hit_phis[hit_rel_id] = hit_phi_odd(hit_Xs[hit_index], hit_Ys[hit_index]);
      }
    }

    // shared_hit_phis
    __syncthreads();

    // Find the permutations given the phis in shared_hit_phis
    for (unsigned int i=0; i<(hit_num + blockDim.x - 1) / blockDim.x; ++i) {
      const auto hit_rel_id = i*blockDim.x + threadIdx.x;
      if (hit_rel_id < hit_num) {
        const auto hit_index = hit_start + hit_rel_id;
        const auto phi = shared_hit_phis[hit_rel_id];
        
        // Find out local position
        unsigned int position = 0;
        for (unsigned int j=0; j<hit_num; ++j) {
          const auto other_phi = shared_hit_phis[j];
          // Stable sorting
          position += phi>other_phi || (phi==other_phi && hit_rel_id>j);
        }
        ASSERT(position < MAX_NUMHITS_IN_MODULE)

        // Store it in hit permutations and in hit_Phis, already ordered
        const auto global_position = hit_start + position;
        hit_permutations[global_position] = hit_index;
        hit_Phis[global_position] = phi;
      }
    }

    // shared_hit_phis
    __syncthreads();
  }
  
}