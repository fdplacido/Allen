#include "SequenceVisitor.cuh"
#include "UTDecodeRawBanksInOrder.cuh"

DEFINE_EMPTY_SET_ARGUMENTS_SIZE(ut_decode_raw_banks_in_order_t)

template<>
void SequenceVisitor::visit<ut_decode_raw_banks_in_order_t>(
  ut_decode_raw_banks_in_order_t& state,
  const RuntimeOptions& runtime_options,
  const Constants& constants,
  argument_manager_t& arguments,
  HostBuffers& host_buffers,
  cudaStream_t& cuda_stream,
  cudaEvent_t& cuda_generic_event)
{
  state.set_opts(dim3(runtime_options.number_of_events, UT::Constants::n_layers), dim3(64), cuda_stream);
  state.set_arguments(
    arguments.offset<dev_ut_raw_input>(),
    arguments.offset<dev_ut_raw_input_offsets>(),
    constants.dev_ut_boards,
    constants.dev_ut_geometry,
    constants.dev_ut_region_offsets,
    constants.dev_unique_x_sector_layer_offsets,
    constants.dev_unique_x_sector_offsets,
    arguments.offset<dev_ut_hit_offsets>(),
    arguments.offset<dev_ut_hits>(),
    arguments.offset<dev_ut_hit_count>(),
    arguments.offset<dev_ut_hit_permutations>()
  );

  state.invoke();
}
