#include "SequenceVisitor.cuh"
#include "PrefixSum.cuh"

template<>
void SequenceVisitor::set_arguments_size<prefix_sum_reduce_scifi_hits_t>(
  const RuntimeOptions& runtime_options,
  const Constants& constants,
  const HostBuffers& host_buffers,
  argument_manager_t& arguments)
{
  const uint total_number_of_zones = runtime_options.number_of_events * SciFi::Constants::n_zones;
  const size_t prefix_sum_auxiliary_array_size = (total_number_of_zones + 511) / 512;
  arguments.set_size<dev_prefix_sum_auxiliary_array_4>(prefix_sum_auxiliary_array_size);
}

DEFINE_EMPTY_SET_ARGUMENTS_SIZE(prefix_sum_single_block_scifi_hits_t)
DEFINE_EMPTY_SET_ARGUMENTS_SIZE(prefix_sum_scan_scifi_hits_t)

template<>
void SequenceVisitor::visit<prefix_sum_reduce_scifi_hits_t>(
  prefix_sum_reduce_scifi_hits_t& state,
  const RuntimeOptions& runtime_options,
  const Constants& constants,
  argument_manager_t& arguments,
  HostBuffers& host_buffers,
  cudaStream_t& cuda_stream,
  cudaEvent_t& cuda_generic_event)
{
  const uint total_number_of_zones = runtime_options.number_of_events * SciFi::Constants::n_zones;
  const size_t prefix_sum_auxiliary_array_size = (total_number_of_zones + 511) / 512;
  state.set_opts(dim3(prefix_sum_auxiliary_array_size), dim3(256), cuda_stream);
  state.set_arguments(
    arguments.offset<dev_scifi_hit_count>(),
    arguments.offset<dev_prefix_sum_auxiliary_array_4>(),
    total_number_of_zones
  );

  state.invoke();
}

template<>
void SequenceVisitor::visit<prefix_sum_single_block_scifi_hits_t>(
  prefix_sum_single_block_scifi_hits_t& state,
  const RuntimeOptions& runtime_options,
  const Constants& constants,
  argument_manager_t& arguments,
  HostBuffers& host_buffers,
  cudaStream_t& cuda_stream,
  cudaEvent_t& cuda_generic_event)
{
  const uint total_number_of_zones = runtime_options.number_of_events * SciFi::Constants::n_zones;
  const size_t prefix_sum_auxiliary_array_size = (total_number_of_zones + 511) / 512;
  state.set_opts(dim3(1), dim3(1024), cuda_stream);
  state.set_arguments(
    arguments.offset<dev_scifi_hit_count>() + total_number_of_zones,
    arguments.offset<dev_prefix_sum_auxiliary_array_4>(),
    prefix_sum_auxiliary_array_size
  );

  state.invoke();
}

template<>
void SequenceVisitor::visit<prefix_sum_scan_scifi_hits_t>(
  prefix_sum_scan_scifi_hits_t& state,
  const RuntimeOptions& runtime_options,
  const Constants& constants,
  argument_manager_t& arguments,
  HostBuffers& host_buffers,
  cudaStream_t& cuda_stream,
  cudaEvent_t& cuda_generic_event)
{
  const uint total_number_of_zones = runtime_options.number_of_events * SciFi::Constants::n_zones;
  const size_t prefix_sum_auxiliary_array_size = (total_number_of_zones + 511) / 512;
  const uint pss_scifi_hits_blocks = prefix_sum_auxiliary_array_size==1 ? 1 : (prefix_sum_auxiliary_array_size-1);
  state.set_opts(dim3(pss_scifi_hits_blocks), dim3(512), cuda_stream);
  state.set_arguments(
    arguments.offset<dev_scifi_hit_count>(),
    arguments.offset<dev_prefix_sum_auxiliary_array_4>(),
    total_number_of_zones
  );

  state.invoke();

  // Fetch total number of hits
  cudaCheck(cudaMemcpyAsync(host_buffers.host_accumulated_number_of_scifi_hits,
    arguments.offset<dev_scifi_hit_count>() + total_number_of_zones,
    sizeof(uint),
    cudaMemcpyDeviceToHost,
    cuda_stream));

  cudaEventRecord(cuda_generic_event, cuda_stream);
  cudaEventSynchronize(cuda_generic_event);

  // info_cout << "Total SciFi cluster estimate: " << *host_accumulated_number_of_scifi_hits << std::endl;
}
