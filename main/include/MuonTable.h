//
// Created by Daniil on 03.04.2019.
//

#ifndef ALLEN_MUONTABLE_H
#define ALLEN_MUONTABLE_H
#include "Tools.h"
#include "MuonTileID.h"
#include "MuonLayout.h"
#include <array>

constexpr unsigned int padGridX[4]{48, 48, 12, 12};
constexpr unsigned int padGridY[4]{8, 8, 8, 8};
constexpr std::array<unsigned int, 16> stripXGridX{48, 48, 48, 48, 48, 48, 48, 48, 12, 12, 12, 12, 12, 12, 12, 12};
constexpr std::array<unsigned int, 16> stripXGridY{1, 2, 2, 2, 1, 2, 2, 2, 8, 2, 2, 2, 8, 2, 2, 2};
constexpr std::array<unsigned int, 16> stripYGridX{8, 4, 2, 2, 8, 4, 2, 2, 12, 4, 2, 2, 12, 4, 2, 2};
constexpr std::array<unsigned int, 16> stripYGridY{8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8};

unsigned int getLayoutX(unsigned int i, unsigned int station, unsigned int region) {
  return (i == 0) ? stripXGridX[station * 4 + region] : stripYGridX[station * 4 +region];
}

unsigned int getLayoutY(unsigned int i, unsigned int station, unsigned int region) {
  return (i == 0) ? stripXGridY[station * 4 + region] : stripYGridY[station * 4 + region];
}

struct MuonTable {
  int gridX[16]{}, gridY[16]{};
  unsigned int offset[16]{};
  float sizeX[16]{}, sizeY[16]{};
  std::vector<std::array<float, 3>> points[4];
};

class MuonTableReader {
public:
  void read(const char* raw_input, MuonTable* pad, MuonTable* stripX, MuonTable* stripY);
};

void calcTilePos(MuonTable* pad, LHCb::MuonTileID& tile, double& x, double& deltax, double& y, double& deltay,
                 double& z);

void calcStripXPos(MuonTable* stripX, LHCb::MuonTileID& tile, double& x, double& deltax, double& y, double& deltay,
                   double& z);

void calcStripYPos(MuonTable* stripY, LHCb::MuonTileID& tile, double& x, double& deltax, double& y, double& deltay,
                   double& z);

int transform_for_uncrossed_hits(LHCb::MuonTileID& tile, MuonTable* pad, MuonTable* stripX, MuonTable* stripY,
                                 double& x, double& dx, double& y, double& dy, double& z);

#endif //ALLEN_MUONTABLE_H
