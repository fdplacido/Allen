#include "LFQualityFilter.cuh"
#include "SequenceVisitor.cuh"

DEFINE_EMPTY_SET_ARGUMENTS_SIZE(lf_quality_filter_t)

template<>
void SequenceVisitor::visit<lf_quality_filter_t>(
  lf_quality_filter_t& state,
  const lf_quality_filter_t::arguments_t& arguments,
  const RuntimeOptions& runtime_options,
  const Constants& constants,
  HostBuffers& host_buffers,
  cudaStream_t& cuda_stream,
  cudaEvent_t& cuda_generic_event)
{
  // Scan:
  // 64 - 28.48%
  // 256 - 17.80%
  // 1024 - 

  state.set_opts(dim3(host_buffers.host_number_of_selected_events[0]), dim3(256), cuda_stream);
  state.set_arguments(
      arguments.offset<dev_scifi_hits>(),
      arguments.offset<dev_scifi_hit_count>(),
      arguments.offset<dev_atomics_velo>(),
      arguments.offset<dev_velo_track_hit_number>(),
      arguments.offset<dev_velo_states>(),
      arguments.offset<dev_atomics_ut>(),
      arguments.offset<dev_ut_track_hits>(),
      arguments.offset<dev_ut_track_hit_number>(),
      arguments.offset<dev_ut_qop>(),
      arguments.offset<dev_ut_track_velo_indices>(),
      arguments.offset<dev_scifi_lf_tracks>(),
      arguments.offset<dev_scifi_lf_atomics>(),
      constants.dev_scifi_geometry,
      constants.dev_inv_clus_res,
      arguments.offset<dev_ut_states>(),
      constants.dev_scifi_tmva1,
      constants.dev_scifi_tmva2,
      constants.dev_scifi_constArrays,
      arguments.offset<dev_atomics_scifi>(),
      arguments.offset<dev_scifi_tracks>());
  
  state.invoke();

  cudaCheck(cudaMemcpyAsync(
    host_buffers.host_atomics_scifi,
    arguments.offset<dev_atomics_scifi>(),
    arguments.size<dev_atomics_scifi>(),
    cudaMemcpyDeviceToHost,
    cuda_stream));

  cudaCheck(cudaMemcpyAsync(
    host_buffers.host_scifi_tracks,
    arguments.offset<dev_scifi_tracks>(),
    arguments.size<dev_scifi_tracks>(),
    cudaMemcpyDeviceToHost,
    cuda_stream));

  cudaEventRecord(cuda_generic_event, cuda_stream);
  cudaEventSynchronize(cuda_generic_event);

}
