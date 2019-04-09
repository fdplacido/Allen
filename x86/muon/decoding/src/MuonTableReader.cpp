#include "MuonTableReader.h"

void MuonTableReader::read(const char* raw_input, MuonTable* pad, MuonTable* stripX, MuonTable* stripY) {
  MuonTable* muonTables[3] = {pad, stripX, stripY};
  for (MuonTable* muonTable : muonTables) {
    size_t gridXSize;
    std::copy_n((size_t*) raw_input, 1, &gridXSize);
    raw_input += sizeof(size_t);
    std::copy_n((int*) raw_input, gridXSize, muonTable -> gridX);
    raw_input += sizeof(int) * gridXSize;

    size_t gridYSize;
    std::copy_n((size_t*) raw_input, 1, &gridYSize);
    raw_input += sizeof(size_t);
    std::copy_n((int*) raw_input, gridYSize, muonTable -> gridY);
    raw_input += sizeof(int) * gridYSize;

    size_t sizeXSize;
    std::copy_n((size_t *) raw_input, 1, &sizeXSize);
    raw_input += sizeof(size_t);
    std::copy_n((float*) raw_input, sizeXSize, muonTable -> sizeX);
    raw_input += sizeof(float) * sizeXSize;

    size_t sizeYSize;
    std::copy_n((size_t*) raw_input, 1, &sizeYSize);
    raw_input += sizeof(size_t);
    std::copy_n((float*) raw_input, sizeYSize, muonTable -> sizeY);
    raw_input += sizeof(float) * sizeYSize;

    size_t offsetSize;
    std::copy_n((size_t*) raw_input, 1, &offsetSize);
    raw_input += sizeof(size_t);
    std::copy_n((unsigned int*) raw_input, offsetSize, muonTable -> offset);
    raw_input += sizeof(unsigned int) * offsetSize;

    size_t tableSize;
    std::copy_n((size_t*) raw_input, 1, &tableSize);
    raw_input += sizeof(size_t);
    for (int i = 0; i < tableSize; i++) {
      size_t stationTableSize;
      std::copy_n((size_t*) raw_input, 1, &stationTableSize);
      raw_input += sizeof(size_t);
      (muonTable -> points)[i].resize(stationTableSize);
      for (int j = 0; j < stationTableSize; j++) {
        float point[3];
        std::copy_n((float *) raw_input, 3, point);
        raw_input += sizeof(float) * 3;
        for (int k = 0; k < 3; k++) {
          (muonTable->points)[i][j][k] = point[k];
        }
      }
    }
  }
}

