#pragma once

// Associate Velo tracks to PVs using their impact parameter and store
// the calculated values.
#include "PV_Definitions.cuh"
#include "AssociateConsolidated.cuh"
#include "Common.h"
#include "Handler.cuh"
#include "ArgumentsCommon.cuh"
#include "ArgumentsVelo.cuh"
#include "ArgumentsPV.cuh"
#include "ArgumentsKalmanFilter.cuh"

__global__ void velo_pv_ip(
  char* dev_kalman_velo_states,
  uint* dev_atomics_velo,
  uint* dev_velo_track_hit_number,
  PV::Vertex* dev_multi_fit_vertices,
  uint* dev_number_of_multi_fit_vertices,
  char* dev_velo_pv_ip);

ALGORITHM(
  velo_pv_ip,
  velo_pv_ip_t,
  ARGUMENTS(
    dev_velo_kalman_beamline_states,
    dev_atomics_velo,
    dev_velo_track_hit_number,
    dev_multi_fit_vertices,
    dev_number_of_multi_fit_vertices,
    dev_velo_pv_ip))
