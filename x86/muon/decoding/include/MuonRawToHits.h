#pragma once

#include <array>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cassert>
#include "MuonTable.h"
#include "MuonRaw.h"
#include "MuonGeometry.h"
#include "MuonDefinitions.cuh"

static const unsigned int MuonDAQHelper_maxTell1Number = 14;

struct Digit {
  Muon::MuonTileID tile;
  unsigned int     tdc;
};
using Digits      = std::vector<Digit>;
using DigitsRange = std::pair<Digits::iterator, Digits::iterator>;

/** @class MuonRawToHits MuonRawToHits.h
 *  This is the muon reconstruction algorithm
 *  This just crosses the logical strips back into pads
 */
class MuonRawToHits {
public:
  MuonRawToHits(MuonTable* pad_, MuonTable* stripX_, MuonTable* stripY_, Muon::MuonGeometry* muonGeometry_) {
    pad = pad_;
    stripX = stripX_;
    stripY = stripY_;
    muonGeometry = muonGeometry_;
  }

  MuonRawToHits() {}

  void operator()(Muon::MuonRawEvent& event, Muon::HitsSoA* hitsSoA) const;

private:
  void decodeTileAndTDC(Muon::MuonRawEvent& , std::array<std::vector<Digit>, 4>& ) const;

  std::array<MuonLayout, 2> makeStripLayouts( const unsigned int, const unsigned int ) const;

  void addCoordsCrossingMap( DigitsRange&, Muon::HitsSoA*, size_t& ) const;

  size_t m_nStations = 4;
  MuonTable* pad;
  MuonTable* stripX;
  MuonTable* stripY;
  Muon::MuonGeometry* muonGeometry;
};
