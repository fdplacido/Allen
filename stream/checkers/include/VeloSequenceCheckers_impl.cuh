#include "ConsolidateVelo.cuh"
#include "TrackChecker.h"
#include "PrepareTracks.h"

/**
 * @brief Specialization for any Velo reconstruction algorithm invoking
 *        consolidate_tracks_t as last step.
 */
template<>
void SequenceVisitor::check<consolidate_velo_tracks_t>(
  HostBuffers& host_buffers,
  const Constants&,
  const CheckerInvoker& checker_invoker,
  MCEvents const& mc_events) const
{
  const auto tracks = prepareVeloTracks(
    host_buffers.host_atomics_velo,
    host_buffers.host_velo_track_hit_number,
    host_buffers.host_velo_track_hits,
    host_buffers.host_number_of_selected_events[0]);

  auto& checker = checker_invoker.checker<TrackCheckerVelo>("\nVelo tracks:", "PrCheckerPlots.root");
  checker.accumulate<TrackCheckerVelo>(mc_events, tracks);
}
