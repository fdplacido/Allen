#include "ParKalmanVeloOnly.cuh"
#include "SequenceVisitor.cuh"

template<>
void SequenceVisitor::set_arguments_size<kalman_velo_only_t>(
  const kalman_velo_only_t& state,
  kalman_velo_only_t::arguments_t arguments,
  const RuntimeOptions& runtime_options,
  const Constants& constants,
  const HostBuffers& host_buffers)
{
  arguments.set_size<dev_kf_tracks>(host_buffers.host_number_of_reconstructed_scifi_tracks[0]);
}

template<>
void SequenceVisitor::visit<kalman_velo_only_t>(
  kalman_velo_only_t& state,
  const kalman_velo_only_t::arguments_t& arguments,
  const RuntimeOptions& runtime_options,
  const Constants& constants,
  HostBuffers& host_buffers,
  cudaStream_t& cuda_stream,
  cudaEvent_t& cuda_generic_event)
{
  state.set_opts(dim3(host_buffers.host_number_of_selected_events[0]), cuda_stream);
  state.set_arguments(
    arguments.offset<dev_atomics_velo>(),
    arguments.offset<dev_velo_track_hit_number>(),
    arguments.offset<dev_velo_track_hits>(),
    arguments.offset<dev_atomics_ut>(),
    arguments.offset<dev_ut_track_hit_number>(),
    arguments.offset<dev_ut_qop>(),
    arguments.offset<dev_ut_track_velo_indices>(),
    arguments.offset<dev_atomics_scifi>(),
    arguments.offset<dev_scifi_track_hit_number>(),
    arguments.offset<dev_scifi_qop>(),
    arguments.offset<dev_scifi_states>(),
    arguments.offset<dev_scifi_track_ut_indices>(),
    arguments.offset<dev_kf_tracks>(),
    constants.dev_scifi_geometry,
    constants.dev_kalman_params);
  state.invoke();
}
