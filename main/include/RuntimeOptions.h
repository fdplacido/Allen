#pragma once

#include <vector>
#include "MuonDefinitions.cuh"

/**
 * @brief Runtime options singleton.
 */
struct RuntimeOptions {
  char* host_velopix_events;
  uint* host_velopix_event_offsets;
  size_t host_velopix_events_size;
  size_t host_velopix_event_offsets_size;
  char* host_ut_events;
  uint* host_ut_event_offsets;
  size_t host_ut_events_size;
  size_t host_ut_event_offsets_size;
  char* host_scifi_events;
  uint* host_scifi_event_offsets;
  size_t host_scifi_events_size;
  size_t host_scifi_event_offsets_size;
  char* host_muon_events;
  uint* host_muon_event_offsets;
  size_t host_muon_events_size;
  size_t host_muon_event_offsets_size;
  uint number_of_events;
  uint number_of_selected_events;
  uint number_of_repetitions;
  bool do_check;
  bool cpu_offload;

  RuntimeOptions() = default;

  RuntimeOptions(
    char* param_host_velopix_events,
    uint* param_host_velopix_event_offsets,
    size_t param_host_velopix_events_size,
    size_t param_host_velopix_event_offsets_size,
    char* param_host_ut_events,
    uint* param_host_ut_event_offsets,
    size_t param_host_ut_events_size,
    size_t param_host_ut_event_offsets_size,
    char* param_host_scifi_events,
    uint* param_host_scifi_event_offsets,
    size_t param_host_scifi_events_size,
    size_t param_host_scifi_event_offsets_size,
    char* param_host_muon_events,
    uint* param_host_muon_event_offsets,
    size_t param_host_muon_events_size,
    size_t param_host_muon_event_offsets_size,
    uint param_number_of_events,
    uint param_number_of_repetitions,
    bool param_do_check,
    bool param_cpu_offload) :
    host_velopix_events(param_host_velopix_events),
    host_velopix_event_offsets(param_host_velopix_event_offsets),
    host_velopix_events_size(param_host_velopix_events_size),
    host_velopix_event_offsets_size(param_host_velopix_event_offsets_size), host_ut_events(param_host_ut_events),
    host_ut_event_offsets(param_host_ut_event_offsets), host_ut_events_size(param_host_ut_events_size),
    host_ut_event_offsets_size(param_host_ut_event_offsets_size), host_scifi_events(param_host_scifi_events),
    host_scifi_event_offsets(param_host_scifi_event_offsets), host_scifi_events_size(param_host_scifi_events_size),
    host_scifi_event_offsets_size(param_host_scifi_event_offsets_size), host_muon_events(param_host_muon_events),
    host_muon_event_offsets(param_host_muon_event_offsets), host_muon_events_size(param_host_muon_events_size),
    host_muon_event_offsets_size(param_host_muon_event_offsets_size), number_of_events(param_number_of_events),
    number_of_selected_events(param_number_of_events), number_of_repetitions(param_number_of_repetitions),
    do_check(param_do_check), cpu_offload(param_cpu_offload)
  {}
};
