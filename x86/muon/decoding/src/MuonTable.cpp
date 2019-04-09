#include "MuonTable.h"

void calcTilePos(MuonTable* pad, Muon::MuonTileID& tile, double& x, double& deltax, double& y, double& deltay,
                 double& z) {
  int          station    = tile.station();
  int          region     = tile.region();
  int          quarter    = tile.quarter();
  int          perQuarter = 3 * padGridX[station] * padGridY[station];
  unsigned int xpad       = tile.nX();
  unsigned int ypad       = tile.nY();
  auto index              = static_cast<unsigned int>((region * 4 + quarter ) * perQuarter);

  if ( ypad < padGridY[station] ) {
    index = index + padGridX[station] * ypad + xpad - padGridX[station];
  } else {
    index = index + padGridX[station] * padGridY[station] +
            2 * padGridX[station] * ( ypad - padGridY[station] ) + xpad;
  }

  auto& p = (pad -> points)[station][index];
  x                        = p[0];
  y                        = p[1];
  z                        = p[2];
  deltax                   = (pad -> sizeX)[station * 4 + region];
  deltay                   = (pad -> sizeY)[station * 4 + region];
}

void calcStripXPos(MuonTable* stripX, Muon::MuonTileID& tile, double& x, double& deltax, double& y, double& deltay,
                   double& z) {
  int station        = tile.station();
  int region         = tile.region();
  int quarter        = tile.quarter();
  int perQuarter     = 3 * stripXGridX[station * 4 + region] * stripXGridY[station * 4 + region];
  unsigned int xpad  = tile.nX();
  unsigned int ypad  = tile.nY();
  unsigned int index = (stripX -> offset)[station * 4 + region] + quarter * perQuarter;

  if ( ypad < stripXGridY[station * 4 + region] ) {
    index = index + stripXGridX[station * 4 + region] * ypad + xpad -
            stripXGridX[station * 4 + region];
  } else {
    index = index + stripXGridX[station * 4 + region] * stripXGridY[station * 4 + region] +
            2 * stripXGridX[station * 4 + region] * ( ypad - stripXGridY[station * 4 + region] ) +
            xpad;
  }

  auto& p = (stripX -> points)[station][index];
  x                        = p[0];
  y                        = p[1];
  z                        = p[2];
  deltax                   = (stripX -> sizeX)[station * 4 + region];
  deltay                   = (stripX -> sizeY)[station * 4 + region];
}

void calcStripYPos(MuonTable* stripY, Muon::MuonTileID& tile, double& x, double& deltax, double& y, double& deltay,
                   double& z) {
  int station        = tile.station();
  int region         = tile.region();
  int quarter        = tile.quarter();
  int perQuarter     = 3 * stripYGridX[station * 4 + region] * stripYGridY[station * 4 + region];
  unsigned int xpad  = tile.nX();
  unsigned int ypad  = tile.nY();
  unsigned int index = (stripY -> offset)[station * 4 + region] + quarter * perQuarter;

  if ( ypad < stripYGridY[station * 4 + region] ) {
    index = index + stripYGridX[station * 4 + region] * ypad + xpad -
            stripYGridX[station * 4 + region];
  } else {
    index = index + stripYGridX[station * 4 + region] * stripYGridY[station * 4 + region] +
            2 * stripYGridX[station * 4 + region] * ( ypad - stripYGridY[station * 4 + region] ) +
            xpad;
  }

  auto& p = (stripY -> points)[station][index];
  x                        = p[0];
  y                        = p[1];
  z                        = p[2];
  deltax                   = (stripY -> sizeX)[station * 4 + region];
  deltay                   = (stripY -> sizeY)[station * 4 + region];
}

int transform_for_uncrossed_hits(Muon::MuonTileID& tile, MuonTable* pad, MuonTable* stripX, MuonTable* stripY,
                                 double& x, double& dx, double& y, double& dy, double& z) {
  unsigned int x1 = getLayoutX(0, tile.station(), tile.region());
  unsigned int y1 = getLayoutY(0, tile.station(), tile.region());
  MuonLayout layoutOne = MuonLayout(x1, y1);
  if (tile.station() > Muon::Constants::n_stations && tile.region() == 0) {
    calcTilePos(pad, tile, x, dx, y, dy, z);
    return 1;
  } else {
    if (tile.layout() == layoutOne) {
      calcStripXPos(stripX, tile, x, dx, y, dy, z);
      return 2;
    } else {
      calcStripYPos(stripY, tile, x, dx, y, dy, z);
      return 3;
    }
  }
}
