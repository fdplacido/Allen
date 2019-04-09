#pragma once

#include <array>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include "MuonTable.h"
#include "raw_bank.hpp"
#include "MuonDefinitions.cuh"
#include <cassert>

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
  MuonRawToHits(MuonTable* pad_, MuonTable* stripX_, MuonTable* stripY_) {
    pad = pad_;
    stripX = stripX_;
    stripY = stripY_;
  }

  MuonRawToHits() {}

  void operator()( const LHCb::RawEvent& event, Muon::HitsSoA* hitsSoA) const;

private:
  void decodeTileAndTDC( const LHCb::RawEvent&, std::array<std::vector<Digit>, 4>& ) const;

  std::array<MuonLayout, 2> makeStripLayouts( const unsigned int, const unsigned int ) const;

  void addCoordsCrossingMap( DigitsRange&, Muon::HitsSoA*, size_t& ) const;

  size_t m_nStations = 4;
  MuonTable* pad;
  MuonTable* stripX;
  MuonTable* stripY;

};
