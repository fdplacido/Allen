/*****************************************************************************\
* (c) Copyright 2000-2018 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#ifndef ALLEN_MUONRAWTOHITS_H
#define ALLEN_MUONRAWTOHITS_H

// Include files
// from STL
#include <array>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include "MuonTable.h"
#include "MuonTileID.h"
#include "RawEvent.h"
#include "CommonMuonHit.h"
#include "RawBank.h"
#include <cassert>

/** @class MuonRawToHits MuonRawToHits.h
 *  This is the muon reconstruction algorithm
 *  This just crosses the logical strips back into pads
 */

static const unsigned int MuonDAQHelper_maxTell1Number = 14;

struct Digit {
  LHCb::MuonTileID tile;
  unsigned int     tdc;
};
using Digits      = std::vector<Digit>;
using DigitsRange = std::pair<Digits::iterator, Digits::iterator>;

class MuonRawToHits {
public:
  MuonRawToHits(MuonTable* pad_, MuonTable* stripX_, MuonTable* stripY_) {
    pad = pad_;
    stripX = stripX_;
    stripY = stripY_;
  }

  MuonRawToHits() {}

  std::array<std::array<CommonMuonHits, 4>, 4> operator()( const LHCb::RawEvent& event ) const;

private:
  void decodeTileAndTDC( const LHCb::RawEvent&, std::array<std::vector<Digit>, 4>& ) const;

  std::array<MuonLayout, 2> makeStripLayouts( const unsigned int, const unsigned int ) const;

  void addCoordsCrossingMap( DigitsRange&, CommonMuonHits& ) const;

  size_t m_nStations = 4;
  MuonTable* pad;
  MuonTable* stripX;
  MuonTable* stripY;

};
#endif // ALLEN_MUONRAWTOHITS_H
