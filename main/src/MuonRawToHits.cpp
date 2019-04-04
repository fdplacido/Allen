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
#include "MuonRawToHits.h"

namespace MuonRawHits {

  enum class ErrorCode {
    BAD_MAGIC = 10,
    BANK_TOO_SHORT,
    PADDING_TOO_LONG,
    TOO_MANY_HITS,
    INVALID_TELL1
  };
  struct ErrorCategory {
    const char* name() const { return "MuonRawBankDecoding"; }
    std::string message( ErrorCode code ) const {
      switch ( code ) {
        case ErrorCode::BAD_MAGIC:
          return "Incorrect Magic pattern in raw bank";
        case ErrorCode::BANK_TOO_SHORT:
          return "Muon bank is too short";
        case ErrorCode::PADDING_TOO_LONG:
          return "Muon bank has too much padding for its size";
        case ErrorCode::TOO_MANY_HITS:
          return "Muon bank has too many hits for its size";
        case ErrorCode::INVALID_TELL1:
          return "Invalid TELL1 source ID";
        default:
          return "unknown exception";
      }
    }
  };
} // namespace MuonRawHits

namespace {
  std::array<std::array<double, 8>, 4> m_xRegions{
      {{-4900, -2400, -1200, -600,
           600, 1200, 2400, 4900},
          {-5252, -2576, -1288, -644,
              644, 1288, 2576, 5252},
          {-5668, -2784, -1392, -696,
              696, 1392, 2784, 5668},
          {-6052, -2976, -1488, -744,
              744, 1488, 2976, 6052}}};

  int nDigits(const LHCb::RawBank &rb) {
    auto range = rb.range<unsigned short>();
    if (range.empty()) return 0;
    auto preamble_size = 2 * ((range[0] + 3) / 2);
    auto overhead = preamble_size + 4;
    return range.size() > overhead ? range.size() - overhead : 0;
  }

  int nDigits(gsl::span<const LHCb::RawBank *> rbs) {
    return std::accumulate(rbs.begin(), rbs.end(), 0,
                           [](int s, const LHCb::RawBank *rb) { return rb ? s + nDigits(*rb) : s; });
  }

  [[gnu::noreturn]] void throw_exception(MuonRawHits::ErrorCode ec, const char *tag) {
    throw tag;
  }
}
#define OOPS( x ) throw_exception( x, __PRETTY_FUNCTION__ )

//=============================================================================
// Main execution
//=============================================================================
std::array<std::array<CommonMuonHits, 4>, 4> MuonRawToHits::operator()( const LHCb::RawEvent& raw ) const {

  const auto&                       mb = raw.banks( LHCb::RawBank::Muon );
  std::array<std::vector<Digit>, 4> decoding;
  for ( auto& decode : decoding ) { decode.reserve( nDigits( mb ) ); }

  std::array<std::array<CommonMuonHits, 4>, 4> commonMuonHitsByStationAndRegion;

  // decode tha data
  decodeTileAndTDC( raw, decoding );

  // sort the digits to ease the crossing
  // the hits come directly sorted by station due to tell1 reading
  // each element of the array represent one station
  constexpr auto regionAndQuarter = []( const Digit& i ) { return i.tile.region() * 4 + i.tile.quarter(); };
  for ( auto& decode : decoding ) {
    std::sort( decode.begin(), decode.end(),
               [&]( const Digit& a, const Digit& b ) { return regionAndQuarter( a ) < regionAndQuarter( b ); } );
  }

  for ( auto& decode : decoding ) {
    CommonMuonHits commonHits;
    commonHits.reserve( decode.size() );
    std::array<DigitsRange, 16> perRegQua;

    unsigned nReg = 0;
    auto     it   = decode.begin();
    for ( auto jt = it; jt != decode.end(); ++jt ) {
      if ( regionAndQuarter( *jt ) != regionAndQuarter( *it ) ) {
        perRegQua[nReg++] = make_pair(it, jt);
        it                = jt;
      }
    }
    perRegQua[nReg++] = make_pair(it, decode.end());

    // do the crossing
    for ( auto& coordsPerRegQua : perRegQua) {
      // return coords directly
      addCoordsCrossingMap( coordsPerRegQua, commonHits );
    }
    auto station = it->tile.station();
    auto region  = m_xRegions.size() - m_nStations + station;

    commonMuonHitsByStationAndRegion[station][region] = std::move(commonHits);
  }
  return commonMuonHitsByStationAndRegion;
}

//=============================================================================

std::array<MuonLayout, 2> MuonRawToHits::makeStripLayouts( const unsigned int station,
                                                           const unsigned int region ) const {
  unsigned int x1 = getLayoutX( 0, station, region );
  unsigned int y1 = getLayoutY( 0, station, region );
  unsigned int x2 = getLayoutX( 1, station, region );
  unsigned int y2 = getLayoutY( 1, station, region );
  if ( x1 > x2 ) {
    return {MuonLayout( x1, y1 ), MuonLayout( x2, y2 )};
  } else {
    return {MuonLayout( x2, y2 ), MuonLayout( x1, y1 )};
  }
}


void MuonRawToHits::addCoordsCrossingMap( DigitsRange& digits, CommonMuonHits& commonHits ) const {
  // need to calculate the shape of the horizontal and vertical logical strips

  // get local MuonLayouts for strips
  const auto& [layoutOne, layoutTwo] = makeStripLayouts( (*(digits.first)).tile.station(), (*(digits.first)).tile.region() );

  // used flags
  std::vector<bool> used( std::distance(digits.first, digits.second), false );

  // partition into the two directions of digits
  // vertical and horizontal stripes
  const auto mid       = std::partition( digits.first, digits.second,
                                         [&layoutOne]( const Digit& digit ) { return digit.tile.layout() == layoutOne; } );
  auto       digitsOne = make_pair( digits.first, mid);
  auto       digitsTwo = make_pair( mid, digits.second );

  // check how many cross
  unsigned i         = 0;
  int      thisGridX = layoutOne.xGrid();
  int      thisGridY = layoutOne.yGrid();

  int otherGridX = layoutTwo.xGrid();
  int otherGridY = layoutTwo.yGrid();
  for (Digits::iterator digits_it = digitsOne.first; digits_it != digitsOne.second; digits_it++) {
    const Digit& one = *digits_it;
    unsigned j = mid - digits.first;
    for (Digits::iterator digits_it2 = digitsTwo.first; digits_it2 != digitsTwo.second; digits_it2++) {
      const Digit& two = *digits_it2;
      if ( ( one.tile.nX() / thisGridX == two.tile.nX() / otherGridX ) &&
           ( one.tile.nY() / thisGridY == two.tile.nY() / otherGridY ) ) {
        unsigned int calcX = one.tile.nX() * otherGridX / thisGridX;
        if ( calcX != two.tile.nX() ) {
          ++j;
          continue;
        }

        unsigned int calcY = two.tile.nY() * thisGridY / otherGridY;
        if ( calcY != one.tile.nY() ) {
          ++j;
          continue;
        }

        LHCb::MuonTileID padTile( one.tile );
        padTile.setY( two.tile.nY() );
        padTile.setLayout( MuonLayout( thisGridX, otherGridY ) );

        double x = 0., dx = 0., y = 0., dy = 0., z = 0., dz = 0.;
        calcTilePos(pad, padTile, x, dx, y, dy, z);

        commonHits.emplace_back( std::move( padTile ), x, dx, y, dy, z, 0, 0, one.tdc, one.tdc - two.tdc );
        used[i] = used[j] = true;
      }
      ++j;
    }
    ++i;
  }

  // copy over "uncrossed" digits

  unsigned m = 0;
  for (Digits::iterator digits_it = digitsOne.first; digits_it != digitsOne.second; digits_it++) {
    Digit& digit = *digits_it;
    if ( !used[m] ) {
      double     x = 0., dx = 0., y = 0., dy = 0., z = 0., dz = 0.;
      if ( digit.tile.station() > ( m_nStations - 3 ) && digit.tile.region() == 0 ) {
        calcTilePos(pad, digit.tile, x, dx, y, dy, z);
      } else {
        calcStripXPos(stripX, digit.tile, x, dx, y, dy, z);
      }
      commonHits.emplace_back( digit.tile, x, dx, y, dy, z, 0, 1, digit.tdc, digit.tdc );
    }
    ++m;
  }
  for (Digits::iterator digits_it = digitsTwo.first; digits_it != digitsTwo.second; digits_it++) {
    Digit& digit = *digits_it;
    if ( !used[m] ) {
      double     x = 0., dx = 0., y = 0., dy = 0., z = 0., dz = 0.;
      if ( digit.tile.station() > ( m_nStations - 3 ) && digit.tile.region() == 0 ) {
        calcTilePos(pad, digit.tile, x, dx, y, dy, z);
      } else {
        calcStripYPos(stripX, digit.tile, x, dx, y, dy, z);
      }
      commonHits.emplace_back( digit.tile, x, dx, y, dy, z, 0, 1, digit.tdc, digit.tdc );
    }
    ++m;
  }

}


void MuonRawToHits::decodeTileAndTDC( const LHCb::RawEvent&             rawdata,
                                            std::array<std::vector<Digit>, 4>& storage ) const {

  // array of vectors of hits
  // each element of the array correspond to hits from a single station
  // this will ease the sorting after

  const auto& mb = rawdata.banks( LHCb::RawBank::Muon );
  if ( mb.empty() ) {
    //warning() << "Ther is no Muon raw Bank in this event" << endmsg;
    return;
  }

  for ( const auto& r : mb ) {
    if ( LHCb::RawBank::MagicPattern != r->magic() ) { OOPS( MuonRawHits::ErrorCode::BAD_MAGIC ); }
    unsigned int tell1Number = r->sourceID();
    if ( tell1Number >= MuonDAQHelper_maxTell1Number ) { OOPS( MuonRawHits::ErrorCode::INVALID_TELL1 ); }

    // decide in which array put the digits according to the Tell1 they come from
    const int inarray = ( tell1Number < 4 ? 0 : tell1Number < 6 ? 1 : tell1Number < 8 ? 2 : 3 );

    // minimum length is 3 words --> 12 bytes
    if ( r->size() < 12 ) { OOPS( MuonRawHits::ErrorCode::BANK_TOO_SHORT ); }
    auto range         = r->range<unsigned short>();
    auto preamble_size = 2 * ( ( range[0] + 3 ) / 2 );
    if ( range.size() < preamble_size ) { OOPS( MuonRawHits::ErrorCode::PADDING_TOO_LONG ); }
    range = range.subspan( preamble_size );

    for ( int i = 0; i < 4; i++ ) {
      //if ( UNLIKELY( range.empty() ) ) { OOPS( MuonRawHits::ErrorCode::BANK_TOO_SHORT ); }
      //if ( UNLIKELY( msgLevel( MSG::VERBOSE ) ) ) { verbose() << " hit in PP " << range[0] << endmsg; }
      //if ( UNLIKELY( range.size() < 1 + range[0] ) ) { OOPS( MuonRawHits::ErrorCode::TOO_MANY_HITS ); }
      for ( unsigned int pp : range.subspan( 1, range[0] ) ) {
        unsigned int     add       = ( pp & 0x0FFF );
        unsigned int     tdc_value = ( ( pp & 0xF000 ) >> 12 );

        //LHCb::MuonTileID tile      = m_muonDetector->getDAQInfo()->getADDInTell1( tell1Number, add );???
        LHCb::MuonTileID tile = LHCb::MuonTileID(0);

        //if ( UNLIKELY( msgLevel( MSG::VERBOSE ) ) ) verbose() << " add " << add << " " << tile << endmsg;
        unsigned int pippo = tile.id();
        if ( pippo != 0 ) {
          //if ( UNLIKELY( msgLevel( MSG::DEBUG ) ) ) debug() << " valid  add " << add << " " << tile << endmsg;
          storage[inarray].emplace_back( Digit{tile, tdc_value} );
        } else {
          //info() << "invalid add " << add << " " << tile << endmsg;
        }
      }
      range = range.subspan( 1 + range[0] );
    }
    assert( range.empty() );
  }
}