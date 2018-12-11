#pragma once

#include "SystemOfUnits.h"

constexpr uint N_LAYERS = VeloUTTracking::n_layers;
constexpr uint NUM_ELEMS = 5 * 2; // num windows * 2

namespace CompassUT {

constexpr uint max_considered_before_found = 4;

}

//=========================================================================
// Point to correct position for windows pointers
//=========================================================================
struct LayerCandidates {
  int from0;
  int size0;
  int from1;
  int size1;
  int from2;
  int size2;
  int from3;
  int size3;
  int from4;
  int size4;
};

struct TrackCandidates {
  const int* m_base_pointer;
  const int  m_num_tracks_event;
  const int  m_track;
  LayerCandidates layers[N_LAYERS];

  __host__ __device__ TrackCandidates(
    const int* base_pointer,
    const int num_tracks_event,
    const int track) : 
    m_base_pointer(base_pointer),
    m_num_tracks_event(num_tracks_event),
    m_track(track) {
      for (int i=0; i<N_LAYERS; ++i) {
        fill_layer(i);
        // printf("TC - t: %i - (%i, %i)(%i, %i)(%i, %i)(%i, %i)(%i, %i)\n",
        //   m_track,
        //   layers[i].from0,
        //   layers[i].size0,
        //   layers[i].from1,
        //   layers[i].size1,
        //   layers[i].from2,
        //   layers[i].size2,
        //   layers[i].from3,
        //   layers[i].size3,
        //   layers[i].from4,
        //   layers[i].size4);
      }
    }

  __host__ __device__ void fill_layer(const int layer) {
    layers[layer].from0 = m_base_pointer[m_track + (m_num_tracks_event * (NUM_ELEMS * layer + 0))];
    layers[layer].size0 = m_base_pointer[m_track + (m_num_tracks_event * (NUM_ELEMS * layer + 1))];
    layers[layer].from1 = m_base_pointer[m_track + (m_num_tracks_event * (NUM_ELEMS * layer + 2))];
    layers[layer].size1 = m_base_pointer[m_track + (m_num_tracks_event * (NUM_ELEMS * layer + 3))];
    layers[layer].from2 = m_base_pointer[m_track + (m_num_tracks_event * (NUM_ELEMS * layer + 4))];
    layers[layer].size2 = m_base_pointer[m_track + (m_num_tracks_event * (NUM_ELEMS * layer + 5))];
    layers[layer].from3 = m_base_pointer[m_track + (m_num_tracks_event * (NUM_ELEMS * layer + 6))];
    layers[layer].size3 = m_base_pointer[m_track + (m_num_tracks_event * (NUM_ELEMS * layer + 7))];
    layers[layer].from4 = m_base_pointer[m_track + (m_num_tracks_event * (NUM_ELEMS * layer + 8))];
    layers[layer].size4 = m_base_pointer[m_track + (m_num_tracks_event * (NUM_ELEMS * layer + 9))];    
  };
};

//=========================================================================
// Save the best q/p, chi2 and number of hits
//=========================================================================
struct BestParams {
  float qp;
  float chi2UT;
  int n_hits;

  __host__ __device__ BestParams () 
  {
    qp = 0.0f;
    chi2UT = PrVeloUTConst::maxPseudoChi2;
    n_hits = 0;
  }
};

/**
   *Constants mainly taken from PrVeloUT.h from Rec
   *  @author Mariusz Witek
   *  @date   2007-05-08
   *  @update for A-Team framework 2007-08-20 SHM
   *
   *  2017-03-01: Christoph Hasse (adapt to future framework)
   *  2018-05-05: Plácido Fernández (make standalone)
   *  2018-07:    Dorothea vom Bruch (convert to C, and then to CUDA code)

 */
namespace VeloUTConst {
  
  // zMidUT is a position of normalization plane which should
  // to be close to z middle of UT ( +- 5 cm ).
  // No need to update with small UT movement.
  static constexpr float zMidUT = 2484.6f;
  //  distToMomentum is properly recalculated in PrUTMagnetTool when B field changes
  static constexpr float distToMomentum = 4.0212e-05f;
  static constexpr float sigmaVeloSlope = 0.10f*Gaudi::Units::mrad;
  static constexpr float invSigmaVeloSlope = 1.0f/sigmaVeloSlope;
  static constexpr float zKink = 1780.0f;
 
  static constexpr float minMomentum =       1.5f * Gaudi::Units::GeV;
  static constexpr float minPT =             0.3f * Gaudi::Units::GeV;
  static constexpr float maxPseudoChi2 =     1280.0f;
  static constexpr float yTol =              0.5f * Gaudi::Units::mm;
  static constexpr float yTolSlope =         0.08f;
  static constexpr float hitTol1 =           6.0f * Gaudi::Units::mm;
  static constexpr float hitTol2 =           0.8f * Gaudi::Units::mm;
  static constexpr float deltaTx1 =          0.035f;
  static constexpr float deltaTx2 =          0.018f;
  static constexpr float maxXSlope =         0.350f;
  static constexpr float maxYSlope =         0.300f;
  static constexpr float centralHoleSize =   33.0f * Gaudi::Units::mm;
  static constexpr float intraLayerDist =    15.0f * Gaudi::Units::mm;
  static constexpr float overlapTol =        0.7f* Gaudi::Units::mm;
  static constexpr float passHoleSize =      40.0f * Gaudi::Units::mm;
  static constexpr int   minHighThres =      1;
  static constexpr bool  printVariables =    false;
  static constexpr bool  passTracks =        false;
  static constexpr bool  doTiming =          false;
  // Scale the z-component, to not run into numerical problems with floats
  // first add to sum values from hit at xMidField, zMidField hit
  static constexpr float zDiff =             0.001f * (zKink - zMidUT);

}