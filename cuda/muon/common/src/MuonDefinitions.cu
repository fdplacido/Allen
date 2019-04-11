#include "MuonDefinitions.cuh"

void Muon::HitsSoA::addAtIndex(size_t index, int tile, float x, float dx, float y, float dy, float z, float dz,
  int uncrossed, unsigned int time, int delta_time, int cluster_size, int region) {
    this -> tile[index] = tile;
    this -> x[index] = x;
    this -> dx[index] = dx;
    this -> y[index] = y;
    this -> dy[index] = dy;
    this -> z[index] = z;
    this -> dz[index] = dz;
    this -> uncrossed[index] = uncrossed;
    this -> time[index] = time;
    this -> delta_time[index] = delta_time;
    this -> cluster_size[index] = cluster_size;
    this -> region[index] = region;
}
