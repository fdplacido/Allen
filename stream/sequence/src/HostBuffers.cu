#include "HostBuffers.cuh"
#include "SciFiDefinitions.cuh"
#include "BeamlinePVConstants.cuh"

void HostBuffers::reserve(const uint max_number_of_events) {
  cudaCheck(cudaMallocHost((void**)&host_velo_tracks_atomics, (2 * max_number_of_events + 1) * sizeof(int)));
  cudaCheck(cudaMallocHost((void**)&host_velo_track_hit_number, max_number_of_events * VeloTracking::max_tracks * sizeof(uint)));
  cudaCheck(cudaMallocHost((void**)&host_velo_track_hits, max_number_of_events * VeloTracking::max_tracks * VeloTracking::max_track_size * sizeof(Velo::Hit)));
  cudaCheck(cudaMallocHost((void**)&host_total_number_of_velo_clusters, sizeof(uint)));
  cudaCheck(cudaMallocHost((void**)&host_number_of_reconstructed_velo_tracks, sizeof(uint)));
  cudaCheck(cudaMallocHost((void**)&host_accumulated_number_of_hits_in_velo_tracks, sizeof(uint)));
  cudaCheck(cudaMallocHost((void**)&host_velo_states, max_number_of_events * VeloTracking::max_tracks * sizeof(Velo::State)));
  cudaCheck(cudaMallocHost((void**)&host_kalmanvelo_states, max_number_of_events * VeloTracking::max_tracks * sizeof(Velo::State)));
  cudaCheck(cudaMallocHost((void**)&host_veloUT_tracks, max_number_of_events * VeloUTTracking::max_num_tracks * sizeof(VeloUTTracking::TrackUT)));
  cudaCheck(cudaMallocHost((void**)&host_atomics_veloUT, VeloUTTracking::num_atomics * max_number_of_events * sizeof(int)));
  cudaCheck(cudaMallocHost((void**)&host_accumulated_number_of_ut_hits, sizeof(uint)));
  cudaCheck(cudaMallocHost((void**)&host_accumulated_number_of_scifi_hits, sizeof(uint)));
  cudaCheck(cudaMallocHost((void**)&host_scifi_tracks, max_number_of_events * SciFi::max_tracks * sizeof(SciFi::Track)));
  cudaCheck(cudaMallocHost((void**)&host_n_scifi_tracks, max_number_of_events * sizeof(uint)));
  cudaCheck(cudaMallocHost((void**)&host_reconstructed_pvs,max_number_of_events * PatPV::max_number_vertices * sizeof(PV::Vertex)));
  cudaCheck(cudaMallocHost((void**)&host_number_of_vertex, max_number_of_events * sizeof(int)));

  cudaCheck(cudaMallocHost((void**)&host_number_of_seeds, max_number_of_events * sizeof(int)));
  cudaCheck(cudaMallocHost((void**)&host_zhisto, max_number_of_events * sizeof(float) * (m_zmax-m_zmin)/m_dz ) );

}


size_t HostBuffers::velo_track_hit_number_size() const {
  return host_number_of_reconstructed_velo_tracks[0] + 1;
}

uint32_t HostBuffers::scifi_hits_uints() const {
  return (sizeof(SciFi::SciFiHit)/sizeof(uint32_t) + 1) * host_accumulated_number_of_scifi_hits[0];
}
