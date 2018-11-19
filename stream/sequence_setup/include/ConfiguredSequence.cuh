#pragma once

#include <tuple>

// All includes of all algorithms
#include "PrefixSum.cuh"
#include "EstimateInputSize.cuh"
#include "MaskedVeloClustering.cuh"
#include "CalculatePhiAndSort.cuh"
#include "SearchByTriplet.cuh"
#include "ConsolidateTracks.cuh"
#include "UTCalculateNumberOfHits.cuh"
#include "UTPreDecode.cuh"
#include "UTFindPermutation.cuh"
#include "UTDecodeRawBanksInOrder.cuh"
#include "VeloUT.cuh"
#include "SciFiCalculateClusterCount.cuh"
#include "SciFiPreDecode.cuh"
#include "SciFiRawBankDecoder.cuh"
#include "PrForward.cuh"
#include "VeloKalmanFilter.cuh"
#include "getSeeds.cuh"
#include "fitSeeds.cuh"
#include "RunForwardCPU.h"
#include "MuonFeaturesExtraction.cuh"


#define SEQUENCE_T(...) typedef std::tuple<__VA_ARGS__> configured_sequence_t;

// SEQUENCE must be defined at compile time.
// Values passed at compile time should match
// the name of the file in "sequences/<filename>.cuh":
//
// "cmake -DSEQUENCE=<sequence_name> .." matches "sequences/<sequence_name>.cuh"
//
// eg.
// "cmake -DSEQUENCE=DefaultSequence .." (or just "cmake ..") matches "sequences/DefaultSequence.cuh"
// "cmake -DSEQUENCE=Velo .." matches "sequences/Velo.cuh"
// "cmake -DSEQUENCE=VeloUT .." matches "sequences/VeloUT.cuh"


#include "sequences/ConfiguredSequence.cuh"

