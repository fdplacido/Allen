#pragma once

#include "SystemOfUnits.h"

namespace LookingForward {
  constexpr float dx_slope = 2e4f;
  constexpr float dx_min = 500.f;
  constexpr float tx_slope = 1250.f;
  constexpr float tx_min = 500.f;
  constexpr float tx_weight = 0.6f;
  constexpr float dx_weight = 0.4f;
  constexpr float max_window_layer0 = 1500.f;
  constexpr float max_window_layer1 = 20.f;
  constexpr float max_window_layer2 = 20.f;
  constexpr float max_window_layer3 = 20.f;
  constexpr float chi2_cut = 5.f;

  /**
   * Station where seeding starts from
   */
  constexpr uint seeding_station = 3;
  constexpr int seeding_first_layer = 8;
  constexpr int seeding_second_layer = 11;

  constexpr float x_diff_layer_qop_offset = 20.f * Gaudi::Units::mm;
  constexpr float x_diff_layer_qop_slope_a = 0.3e6;
  constexpr float x_diff_layer_qop_slope_b = 0.2e6;

  // cut on the difference between tx from the extrapolation and
  // tx from the hits in the two x layers
  constexpr float max_tx_diff = 0.05f * Gaudi::Units::mm;

  // z distance between various layers of a station
  constexpr float dz_layers_station = 70. * Gaudi::Units::mm;
  constexpr float dz_x_layers = 3.f * dz_layers_station;
  constexpr float dz_x_u_layers = 1.f * dz_layers_station;
  constexpr float dz_u_v_layers = 1.f * dz_layers_station;
  constexpr float dz_x_v_layers = 2.f * dz_layers_station;

  // z at the center of the magnet
  constexpr float z_magnet = 5212.38f;

  // z distance between various layers of different stations
  constexpr float dz_x_T1_0_T2_0 = 682 * Gaudi::Units::mm;
  constexpr float dz_x_T1_0_T2_3 = 892 * Gaudi::Units::mm;
  constexpr float dz_x_T1_0_T3_0 = 1367 * Gaudi::Units::mm;
  constexpr float dz_x_T1_0_T3_3 = 1577 * Gaudi::Units::mm;

  // cut on x difference between x- and u-/v-layers
  constexpr float dx_x_uv_layers = 200.f * Gaudi::Units::mm;
  constexpr float dx_x_uv_layers_slope = 2.f * Gaudi::Units::mm;

  // cut on x difference between T1 and T2/T3 x-layers
  constexpr float dx_x_T2_T3_offset = 500 * Gaudi::Units::mm;
  constexpr float dx_x_T2_T3_slope = 6.e6f;

  constexpr float z_last_UT_plane = 2642.f;

  // combinatorics cut-offs, to be tuned!!
  // max # of quadruplets per veloUT input track
  constexpr int max_quadruplets = 100;

  // detector limits
  constexpr float xMin = -4090.f;
  constexpr float xMax = 4090.f;
  constexpr float yUpMin = -50.f;
  constexpr float yUpMax = 3030.f;
  constexpr float yDownMin = -3030.f;
  constexpr float yDownMax = 50.f;

  struct Constants {
    float Zone_zPos[12] = {7826., 7896., 7966., 8036., 8508., 8578., 8648., 8718., 9193., 9263., 9333., 9403.};

    float zMagnetParams[4] = {5212.38, 406.609, -1102.35, -498.039};

    float ds_p_param[12] =
      // updated value for layer 8
      {1307.09, 1288.22, 899.152, 1304.91, 1304.63, 1293.6, 50.6114, 1303.6, 1207.72, 1297.08, 1299.11, 1308.75};
    // original values from the velo
    // {1307.09, 1288.22, 899.152, 1304.91, 1304.63, 1293.6, 50.6114, 1303.6, 1309.73, 1297.08, 1299.11, 1308.75};

    float dp_y_mag_plus[12][3] = {
      {-9.57809, 0.0665787, -0.000144284},
      {-16.1944, 0.191599, -0.000158976},
      {-146.023, -0.393702, -6.67197e-05},
      {-10.1971, 0.0698613, -0.000158352},
      {-13.3151, 0.0820391, -0.000193986},
      {-18.2638, 0.213891, -0.000199823},
      {-529.343, -1.32182, 0.000139385},
      {-13.6997, 0.0828196, -0.000196252},
      {
        -0.535394,
        0.0347712,
        -0.000186249,
      },
      {-23.6321, 0.229475, -0.000238147},
      {-23.6499, -0.068727, -0.000218933},
      {-16.7274, 0.0898987, -0.000222276},
    };
    float dp_y_mag_minus[12][3] = {
      {7.77962, -0.0671629, 0.000148067},
      {14.8172, 0.0731748, 0.000152873},
      {138.486, 0.151621, 4.48781e-05},
      {8.5616, -0.0718798, 0.000163705},
      {10.8283, -0.0812554, 0.000196347},
      {17.4195, 0.0633952, 0.000197882},
      {507.165, 1.13525, -0.000218559},
      {11.5582, -0.0861381, 0.00020702},
      {
        1.68943,
        -0.0542836,
        0.000202545,
      },
      {22.9984, 0.0674907, 0.000229801},
      {20.8404, -0.230712, 0.000236968},
      {13.1927, -0.0868799, 0.000227065},
    };
    float dp_x_mag_plus[12][5] = {{
                                          3.03633,
                                          -0.0129859,
                                          -5.12209e-06,
                                          7.35636e-09,
                                          3.00793e-12,
                                        },
                                        {
                                          1.58566,
                                          -0.00722185,
                                          -3.39866e-06,
                                          6.82563e-09,
                                          1.37203e-13,
                                        },
                                        {
                                          33.9066,
                                          0.0410481,
                                          -0.000210066,
                                          7.00692e-08,
                                          1.0448e-10,
                                        },
                                        {
                                          3.26,
                                          -0.0127267,
                                          -4.65083e-06,
                                          7.75563e-09,
                                          1.88199e-12,
                                        },
                                        {
                                          4.08483,
                                          -0.0130299,
                                          -6.08248e-06,
                                          8.9185e-09,
                                          1.92318e-12,
                                        },
                                        {
                                          2.95945,
                                          -0.00471798,
                                          -1.13934e-05,
                                          5.22784e-09,
                                          4.91834e-12,
                                        },
                                        {
                                          126.581,
                                          0.211745,
                                          -0.000884383,
                                          2.57066e-07,
                                          4.36623e-10,
                                        },
                                        {
                                          5.23607,
                                          -0.00723468,
                                          -1.02187e-05,
                                          7.26014e-09,
                                          3.22778e-12,
                                        },
                                        {
                                          9.24386,
                                          -0.0554029,
                                          -1.57925e-05,
                                          3.09113e-08,
                                          1.15065e-11,
                                        },
                                        {
                                          4.05552,
                                          -0.00266585,
                                          -1.52019e-05,
                                          7.85394e-09,
                                          5.00958e-12,
                                        },
                                        {
                                          9.89762,
                                          -0.0140104,
                                          -2.08032e-05,
                                          1.13114e-08,
                                          7.47106e-12,
                                        },
                                        {
                                          6.59536,
                                          -0.0149863,
                                          -1.75225e-05,
                                          1.2962e-08,
                                          6.08399e-12,
                                        }};

    float dp_x_mag_minus[12][5] = {{
                                           -2.43899,
                                           -0.0129165,
                                           2.52919e-06,
                                           6.16789e-09,
                                           -1.09957e-12,
                                         },
                                         {
                                           -6.08958,
                                           -0.0134028,
                                           7.64373e-06,
                                           8.50007e-09,
                                           -2.39752e-12,
                                         },
                                         {
                                           -27.1364,
                                           0.04193,
                                           0.000193579,
                                           7.36596e-08,
                                           -9.7807e-11,
                                         },
                                         {
                                           -2.68477,
                                           -0.0122401,
                                           3.12353e-06,
                                           5.43385e-09,
                                           -1.48007e-12,
                                         },
                                         {
                                           -3.69344,
                                           -0.0131073,
                                           4.58863e-06,
                                           6.65455e-09,
                                           -1.72323e-12,
                                         },
                                         {
                                           -6.98459,
                                           -0.0111156,
                                           6.42559e-06,
                                           7.09871e-09,
                                           -1.51384e-12,
                                         },
                                         {
                                           -116.517,
                                           0.174343,
                                           0.000862418,
                                           3.34216e-07,
                                           -4.77696e-10,
                                         },
                                         {
                                           -3.70516,
                                           -0.00888972,
                                           2.62343e-06,
                                           6.22499e-09,
                                           7.23911e-13,
                                         },
                                         {
                                           -13.5286,
                                           -0.0546777,
                                           1.94416e-05,
                                           3.40368e-08,
                                           -1.90795e-11,
                                         },
                                         {
                                           -9.13909,
                                           -0.0155821,
                                           1.37127e-05,
                                           1.1462e-08,
                                           -4.09602e-12,
                                         },
                                         {
                                           -2.99785,
                                           -0.00743007,
                                           5.08675e-06,
                                           1.04737e-08,
                                           1.4279e-12,
                                         },
                                         {
                                           -4.54212,
                                           -0.0182558,
                                           6.24697e-06,
                                           1.3144e-08,
                                           -9.98149e-13,
                                         }};
    float dp_plus_offset[12] =
      {19.3354, 22.3182, 46.002, 21.7518, 27.7576, 30.3437, 552.816, 28.3234, 37.037, 40.018, 37.6157, 35.6917};
    float dp_minus_offset[12] = {-18.1339,
                                       -21.0436,
                                       -45.0591,
                                       -20.4936,
                                       -25.8877,
                                       -28.3128,
                                       -444.456,
                                       -27.4285,
                                       -35.0231,
                                       -36.6656,
                                       -38.7877,
                                       -34.5437};

    // dxdy for the layers in one station
    float Zone_dxdy[4] = {0, 0.0874892, -0.0874892, 0};
  };
} // namespace LookingForward
