#include "SearchByTriplet.cuh"
#include "ClusteringDefinitions.cuh"

__constant__ float Configuration::velo_search_by_triplet_t::forward_phi_tolerance;

__constant__ float Configuration::velo_search_by_triplet_t::max_chi2;

__constant__ float Configuration::velo_search_by_triplet_t::max_scatter_forwarding;
__constant__ float Configuration::velo_search_by_triplet_t::max_scatter_seeding;

__constant__ uint Configuration::velo_search_by_triplet_t::max_skipped_modules;

__constant__ uint Configuration::velo_search_by_triplet_t::max_weak_tracks;

__constant__ float Configuration::velo_search_by_triplet_t::phi_extrapolation_base;
__constant__ float Configuration::velo_search_by_triplet_t::phi_extrapolation_coef;
__constant__ uint Configuration::velo_search_by_triplet_t::ttf_modulo;

__constant__ int Configuration::velo_search_by_triplet_t::ttf_modulo_mask;

/**
 * @brief Track forwarding algorithm based on triplet finding
 * @detail For details, check out paper
 *         "A fast local algorithm for track reconstruction on parallel architectures"
 *
 *         Note: All hit arrays are contained in the dev_velo_cluster_container.
 *               By having a single array and offsetting it every time, less registers are
 *               required.
 *
 *               Hereby all the pointers and how to access them:
 *
 *         const float* hit_Xs = (float*) (dev_velo_cluster_container + 5 * number_of_hits + hit_offset);
 *         const float* hit_Ys = (float*) (dev_velo_cluster_container + hit_offset);
 *         const float* hit_Zs = (float*) (dev_velo_cluster_container + number_of_hits + hit_offset);
 *         const float* hit_Phis = (float*) (dev_velo_cluster_container + 4 * number_of_hits + hit_offset);
 *
 *         Note: Atomics is another case where we need several variables from an array.
 *               We just keep dev_atomics_velo, and access the required ones upon request.
 *
 *               Below are all atomics used by this algorithm:
 *
 *         const int ip_shift = gridDim.x + blockIdx.x * (Velo::num_atomics - 1);
 *         uint* tracks_insert_pointer = (uint*) dev_atomics_velo + event_number;
 *         uint* weaktracks_insert_pointer = (uint*) dev_atomics_velo + ip_shift;
 *         uint* tracklets_insert_pointer = (uint*) dev_atomics_velo + ip_shift + 1;
 *         uint* ttf_insert_pointer = (uint*) dev_atomics_velo + ip_shift + 2;
 *         uint* local_number_of_hits = (uint*) dev_atomics_velo + ip_shift + 3;
 */

__global__ void search_by_triplet(
  uint32_t* dev_velo_cluster_container,
  uint* dev_module_cluster_start,
  uint* dev_module_cluster_num,
  Velo::TrackHits* dev_tracks,
  Velo::TrackletHits* dev_tracklets,
  uint* dev_tracks_to_follow,
  Velo::TrackletHits* dev_weak_tracks,
  bool* dev_hit_used,
  uint* dev_atomics_velo,
  short* dev_h0_candidates,
  short* dev_h2_candidates,
  unsigned short* dev_rel_indices,
  const VeloGeometry* dev_velo_geometry)
{
  /* Data initialization */
  // Each event is treated with two blocks, one for each side.
  const uint event_number = blockIdx.x;
  const uint number_of_events = gridDim.x;
  const uint tracks_offset = event_number * Velo::Constants::max_tracks;

  // Pointers to data within the event
  const uint number_of_hits = dev_module_cluster_start[Velo::Constants::n_modules * number_of_events];
  const uint* module_hitStarts = dev_module_cluster_start + event_number * Velo::Constants::n_modules;
  const uint* module_hitNums = dev_module_cluster_num + event_number * Velo::Constants::n_modules;
  const uint hit_offset = module_hitStarts[0];
  assert((module_hitStarts[52] - module_hitStarts[0]) < Velo::Constants::max_number_of_hits_per_event);

  // Per event datatypes
  Velo::TrackHits* tracks = dev_tracks + tracks_offset;

  // Per side datatypes
  bool* hit_used = dev_hit_used + hit_offset;
  short* h0_candidates = dev_h0_candidates + 2 * hit_offset;
  short* h2_candidates = dev_h2_candidates + 2 * hit_offset;

  uint* tracks_to_follow = dev_tracks_to_follow + event_number * Configuration::velo_search_by_triplet_t::ttf_modulo;
  Velo::TrackletHits* weak_tracks =
    dev_weak_tracks + event_number * Configuration::velo_search_by_triplet_t::max_weak_tracks;
  Velo::TrackletHits* tracklets = dev_tracklets + event_number * Configuration::velo_search_by_triplet_t::ttf_modulo;
  unsigned short* h1_rel_indices = dev_rel_indices + event_number * Velo::Constants::max_numhits_in_module;

  // Shared memory size is defined externally
  __shared__ float module_data[12];

  process_modules(
    (Velo::Module*) &module_data[0],
    hit_used,
    h0_candidates,
    h2_candidates,
    module_hitStarts,
    module_hitNums,
    (float*) dev_velo_cluster_container + hit_offset,
    tracks_to_follow,
    weak_tracks,
    tracklets,
    tracks,
    number_of_hits,
    h1_rel_indices,
    hit_offset,
    dev_velo_geometry->module_zs,
    dev_atomics_velo);
}
