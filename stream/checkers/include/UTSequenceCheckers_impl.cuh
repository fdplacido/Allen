#include "VeloUT.cuh"

/**
 * @brief Specialization for any Velo reconstruction algorithm invoking
 *        consolidate_ut_tracks_t as last step.
 */
template<>
void SequenceVisitor::check<consolidate_ut_tracks_t>(
  const uint& start_event_offset,
  const uint& number_of_events_requested,
  const HostBuffers& host_buffers,
  const Constants& constants,
  const CheckerInvoker& checker_invoker) const
{
  info_cout << "Checking Velo+UT tracks" << std::endl;

  const auto tracks = prepareUTTracks(
    host_buffers.host_atomics_velo,
    host_buffers.host_velo_track_hit_number,
    host_buffers.host_velo_track_hits,
    host_buffers.host_atomics_ut,
    host_buffers.host_ut_track_hit_number,
    host_buffers.host_ut_track_hits,
    host_buffers.host_ut_track_velo_indices,
    host_buffers.host_ut_qop,
    number_of_events_requested);

  std::vector< std::vector< std::vector< uint32_t > > > scifi_ids_events = checker_invoker.check<TrackCheckerVeloUT>(start_event_offset, tracks);
  
  for ( int j = 0; j < number_of_events_requested; j++ ) {
    debug_cout << "UT: At event " << j << " found " << scifi_ids_events[j].size() << " tracks" <<  std::endl;
    for ( int i = 0; i < scifi_ids_events[j].size(); ++i ) {
      debug_cout << "at track " << i << " found " << scifi_ids_events[j][i].size() << " IDs" << std::endl;
      // for ( const auto id : scifi_ids_events[j][i] ) {
      //   debug_cout << "\t ID = " << std::hex << id << std::dec << std::endl;
      // }
    }
  }
  
}
