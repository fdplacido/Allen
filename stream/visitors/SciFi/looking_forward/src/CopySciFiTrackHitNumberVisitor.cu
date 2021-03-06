#include "SequenceVisitor.cuh"
#include "PrefixSum.cuh"

template<>
void SequenceVisitor::set_arguments_size<copy_scifi_track_hit_number_t>(
  const copy_scifi_track_hit_number_t& state,
  copy_scifi_track_hit_number_t::arguments_t arguments,
  const RuntimeOptions& runtime_options,
  const Constants& constants,
  const HostBuffers& host_buffers)
{
  arguments.set_size<dev_scifi_track_hit_number>(host_buffers.scifi_track_hit_number_size());
}

template<>
void SequenceVisitor::visit<copy_scifi_track_hit_number_t>(
  copy_scifi_track_hit_number_t& state,
  const copy_scifi_track_hit_number_t::arguments_t& arguments,
  const RuntimeOptions& runtime_options,
  const Constants& constants,
  HostBuffers& host_buffers,
  cudaStream_t& cuda_stream,
  cudaEvent_t& cuda_generic_event)
{
  state.set_opts(dim3(host_buffers.host_number_of_selected_events[0]), cuda_stream);
  state.set_arguments(
    arguments.offset<dev_atomics_ut>(),
    arguments.offset<dev_scifi_tracks>(),
    arguments.offset<dev_atomics_scifi>(),
    arguments.offset<dev_scifi_track_hit_number>());

  state.invoke();
}
