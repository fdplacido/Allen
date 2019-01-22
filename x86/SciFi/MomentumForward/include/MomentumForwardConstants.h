#pragma once

namespace SciFi {
  namespace MomentumForward {
    
    // cut on difference between extrapolated x position and x position in 
    // first layer of T1 / last layer of T3 as a function of qop
    // The cut is a straight line with offset and slope
    constexpr float dx_extrap_qop_offset_T1 = 20.f;
    constexpr float dx_extrap_qop_slope_T1 = 1.e6;
    constexpr float dx_extrap_qop_offset_T3 = 40.f;
    constexpr float dx_extrap_qop_slope_T3 = 1.5e6;
    
    // cut on difference between the x hit positions of the two x layers
    // in one station
    // The cut is based on two straight lines with different slopes containing 
    // the interesting region
    constexpr float x_diff_layer_qop_offset = 20.f;
    constexpr float x_diff_layer_qop_slope_a = 0.3e6;
    constexpr float x_diff_layer_qop_slope_b = 0.2e6;

    // cut on the difference between tx from the extrapolation and
    // tx from the hits in the two x layers
    constexpr float max_tx_diff = 0.05f;

    // z distance between the x layers of a station
    constexpr float dz_x_layers = 210.f;
    
  }
}
