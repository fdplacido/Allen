#pragma once

#include "MuonTileID.h"
#include "MuonLayout.h"
#include "MuonDefinitions.cuh"
#include <array>

constexpr unsigned int padGridX[4]{48, 48, 12, 12};
constexpr unsigned int padGridY[4]{8, 8, 8, 8};
constexpr std::array<unsigned int, 16> stripXGridX{48, 48, 48, 48, 48, 48, 48, 48, 12, 12, 12, 12, 12, 12, 12, 12};
constexpr std::array<unsigned int, 16> stripXGridY{1, 2, 2, 2, 1, 2, 2, 2, 8, 2, 2, 2, 8, 2, 2, 2};
constexpr std::array<unsigned int, 16> stripYGridX{8, 4, 2, 2, 8, 4, 2, 2, 12, 4, 2, 2, 12, 4, 2, 2};
constexpr std::array<unsigned int, 16> stripYGridY{8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8};

inline unsigned int getLayoutX(unsigned int i, unsigned int station, unsigned int region) {
  return (i == 0) ? stripXGridX[station * 4 + region] : stripYGridX[station * 4 +region];
}

inline unsigned int getLayoutY(unsigned int i, unsigned int station, unsigned int region) {
  return (i == 0) ? stripXGridY[station * 4 + region] : stripYGridY[station * 4 + region];
}

struct MuonTable {
  int gridX[16]{}, gridY[16]{};
  unsigned int offset[16]{};
  float sizeX[16]{}, sizeY[16]{};
  std::vector<std::array<float, 3>> points[4];
};

inline void calcTilePos(MuonTable* pad, Muon::MuonTileID& tile, double& x, double& deltax, double& y, double& deltay,
                 double& z);

inline void calcStripXPos(MuonTable* stripX, Muon::MuonTileID& tile, double& x, double& deltax, double& y,
    double& deltay, double& z);

inline void calcStripYPos(MuonTable* stripY, Muon::MuonTileID& tile, double& x, double& deltax, double& y,
    double& deltay, double& z);

inline int transform_for_uncrossed_hits(Muon::MuonTileID& tile, MuonTable* pad, MuonTable* stripX, MuonTable* stripY,
                                 double& x, double& dx, double& y, double& dy, double& z);

