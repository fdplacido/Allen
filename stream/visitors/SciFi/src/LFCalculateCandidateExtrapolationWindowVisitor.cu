#include "LFCalculateCandidateExtrapolationWindow.cuh"
#include "SequenceVisitor.cuh"

template<>
void SequenceVisitor::set_arguments_size<lf_calculate_candidate_extrapolation_window_t>(
  lf_calculate_candidate_extrapolation_window_t::arguments_t arguments,
  const RuntimeOptions& runtime_options,
  const Constants& constants,
  const HostBuffers& host_buffers)
{
  arguments.set_size<dev_extrapolation_layer_candidates>(
    2 * host_buffers.host_number_of_selected_events[0] * SciFi::Constants::max_track_candidates);
}

template<>
void SequenceVisitor::visit<lf_calculate_candidate_extrapolation_window_t>(
  lf_calculate_candidate_extrapolation_window_t& state,
  const lf_calculate_candidate_extrapolation_window_t::arguments_t& arguments,
  const RuntimeOptions& runtime_options,
  const Constants& constants,
  HostBuffers& host_buffers,
  cudaStream_t& cuda_stream,
  cudaEvent_t& cuda_generic_event)
{
  cudaCheck(
    cudaMemsetAsync(
      arguments.offset<dev_extrapolation_layer_candidates>(),
      0,
      arguments.size<dev_extrapolation_layer_candidates>(),
      cuda_stream));

  state.set_opts(dim3(host_buffers.host_number_of_selected_events[0]), dim3(256), cuda_stream);
  state.set_arguments(
    arguments.offset<dev_scifi_hits>(),
    arguments.offset<dev_scifi_hit_count>(),
    arguments.offset<dev_atomics_ut>(),
    arguments.offset<dev_scifi_track_candidates>(),
    arguments.offset<dev_atomics_scifi>(),
    constants.dev_scifi_geometry,
    constants.dev_looking_forward_constants,
    constants.dev_inv_clus_res,
    arguments.offset<dev_ut_states>(),
    7,
    arguments.offset<dev_extrapolation_layer_candidates>());

  state.invoke();

  std::vector<unsigned short> extrapolation_layer_candidates (arguments.size<dev_extrapolation_layer_candidates>() / sizeof(unsigned short));
  std::vector<int> scifi_atomics(arguments.size<dev_atomics_scifi>() / sizeof(int));

  cudaCheck(cudaMemcpyAsync(
    scifi_atomics.data(),
    arguments.offset<dev_atomics_scifi>(),
    arguments.size<dev_atomics_scifi>(),
    cudaMemcpyDeviceToHost,
    cuda_stream));

  cudaCheck(cudaMemcpyAsync(extrapolation_layer_candidates.data(),
    arguments.offset<dev_extrapolation_layer_candidates>(),
    arguments.size<dev_extrapolation_layer_candidates>(),
    cudaMemcpyDeviceToHost,
    cuda_stream));

  cudaEventRecord(cuda_generic_event, cuda_stream);
  cudaEventSynchronize(cuda_generic_event);

  for (uint i=0; i<host_buffers.host_number_of_selected_events[0]; ++i) {
    const auto number_of_candidates = scifi_atomics[i];
    info_cout << "Event " << i << ":" << std::endl;

    for (int j=0; j<number_of_candidates; ++j) {
      info_cout << " " << extrapolation_layer_candidates[i * SciFi::Constants::max_track_candidates + j]
        << ", " << extrapolation_layer_candidates[
            host_buffers.host_number_of_selected_events[0] * SciFi::Constants::max_track_candidates + 
            i * SciFi::Constants::max_track_candidates + j] << std::endl;
    }
  }
}
