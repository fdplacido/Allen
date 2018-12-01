#include "blpv_extrapolate.cuh"

__global__ void blpv_extrapolate(
  char* dev_kalmanvelo_states,
  int* dev_atomics_storage,
  uint* dev_velo_track_hit_number,
  PVTrack* dev_pvtracks)
{
  const uint number_of_events = gridDim.x;
  const uint event_number = blockIdx.x;

  const Velo::Consolidated::Tracks velo_tracks {
    (uint*) dev_atomics_storage, dev_velo_track_hit_number, event_number, number_of_events};
  const Velo::Consolidated::States velo_states =
    Velo::Consolidated::States(dev_kalmanvelo_states, velo_tracks.total_number_of_tracks);
  const uint number_of_tracks_event = velo_tracks.number_of_tracks(event_number);
  const uint event_tracks_offset = velo_tracks.tracks_offset(event_number);

  for (int i = 0; i < number_of_tracks_event / blockDim.x + 1; i++) {
    int index = blockDim.x * i + threadIdx.x;
    if (index < number_of_tracks_event) {
      VeloState s = velo_states.get(event_tracks_offset + index);
      PatPV::XYZPoint beamline(0., 0., 0.);
      const auto tx = s.tx;
      const auto ty = s.ty;
      const double dz = (tx * (beamline.x - s.x) + ty * (beamline.y - s.y)) / (tx * tx + ty * ty);
      const double newz = s.z + dz;
      PVTrack pvtrack = PVTrack {s, dz, index};
      dev_pvtracks[event_tracks_offset + index] = pvtrack;
    }
  }
}