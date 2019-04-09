#pragma once

#include <iostream>
#include <vector>

/** @class MuonLayout MuonLayout.h Kernel/MuonLayout.h

   Defines a Muon station single plane logical layout. Layouts in
   all the regions of the station are the same with the scaling factor 2
   when passing from a region to a larger one. The class also implements
   various layout/tiles manipulation functions.

   @author A.Tsaregorodtsev
   @date 6 April 2001
*/

#include "MuonTileID.h"

namespace LHCb {
  class MuonTileID;
}

class MuonLayout {

public:
  /// Default constructor
  MuonLayout() {
    m_xgrid = 0;
    m_ygrid = 0;
  }

  /** Constructor taking X and Y granularities
    @param   xgrid  granularity in X
    @param   ygrid  granularity in Y
  */
 MuonLayout(unsigned int xgrid, unsigned int ygrid) {
    m_xgrid = xgrid;
    m_ygrid = ygrid;
  }

  /** Constructor taking X and Y granularities as std::pair
    @param   xygrid  granularity in X
  */
  MuonLayout(std::pair<unsigned int, unsigned int> xygrid) {
    m_xgrid = xygrid.first;
    m_ygrid = xygrid.second;
  }

  /// Accessor to X granularity
  unsigned int xGrid() const { return m_xgrid; }

  /// Accessor to Y granularity
  unsigned int yGrid() const { return m_ygrid; }

private:
  unsigned int m_xgrid;
  unsigned int m_ygrid;
};

inline bool operator==(const MuonLayout &ml1, const MuonLayout &ml2) {
  return ml1.xGrid() == ml2.xGrid() &&
         ml1.yGrid() == ml2.yGrid();
}

inline bool operator!=(const MuonLayout &ml1, const MuonLayout &ml2) {
  return !(ml1 == ml2);
}
