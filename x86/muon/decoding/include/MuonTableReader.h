#pragma once

#include "MuonTable.h"
#include <algorithm>

class MuonTableReader {
public:
  void read(const char* raw_input, MuonTable* pad, MuonTable* stripX, MuonTable* stripY);
};
