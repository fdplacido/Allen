#include "RateChecker.h"

void checkHlt1Rate(
  const bool* one_track_decisions,
  const bool* two_track_decisions,
  const int* track_atomics,
  const uint* sv_atomics,
  const uint number_of_events)
{

  // Event counters.
  uint n_evts_one_track = 0;
  uint n_evts_two_track = 0;
  uint n_evts_inc = 0;
  
  // Event loop.
  for (uint i_event = 0; i_event < number_of_events; i_event++) {
    
    // Check one track decisions.
    bool one_track_pass = false;
    const int* event_tracks_offsets = track_atomics + number_of_events;
    const bool* event_one_track_decisions = one_track_decisions + event_tracks_offsets[i_event];
    const int n_tracks_event = track_atomics[i_event];
    for (int i_track = 0; i_track < n_tracks_event; i_track++) {
      if (event_one_track_decisions[i_track]) {
        one_track_pass = true;
        break;
      }
    }

    // Check two track decisions.
    bool two_track_pass = false;
    const bool* event_two_track_decisions = two_track_decisions + sv_atomics[i_event];
    const int n_svs_event = sv_atomics[i_event + 1] - sv_atomics[i_event];
    for (int i_sv = 0; i_sv < n_svs_event; i_sv++) {
      if (event_two_track_decisions[i_sv]) {
        two_track_pass = true;
        break;
      }
    }

    n_evts_one_track += one_track_pass;
    n_evts_two_track += two_track_pass;
    n_evts_inc += one_track_pass || two_track_pass;
    
  }

  // Get the total number of passing SVs.
  uint n_total_sv_passes = 0;
  for (int i_sv = 0; i_sv < sv_atomics[number_of_events]; i_sv++) {
    if (two_track_decisions[i_sv]) n_total_sv_passes++;
  }

  printf("Total events: %i \n", number_of_events);
  printf("One track: %i\n", n_evts_one_track);
  printf("Two track: %i\n", n_evts_two_track);
  printf("Inclusive: %i\n\n", n_evts_inc);
}
