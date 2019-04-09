/*
void commonMuonHitsToHitsSoA(std::array<std::array<CommonMuonHits, 4>, 4>& commonMuonHits, Muon::HitsSoA* muon_station_hits) {
  unsigned int i = 0;
  unsigned int currentStation = 0;
  unsigned int hitsPerStation = 0;
  unsigned int offset = 0;
  for (auto commonMuonHitsPerStation: commonMuonHits) {
    hitsPerStation = 0;
    for (auto commonMuonHitsPerRegion: commonMuonHitsPerStation) {
      for (auto commonMuonHit: commonMuonHitsPerRegion) {
        muon_station_hits->tile[i] = commonMuonHit.tile().id();
        muon_station_hits->x[i] = commonMuonHit.x();
        muon_station_hits->dx[i] = commonMuonHit.dx();
        muon_station_hits->y[i] = commonMuonHit.y();
        muon_station_hits->dy[i] = commonMuonHit.dy();
        muon_station_hits->z[i] = commonMuonHit.z();
        muon_station_hits->dz[i] = commonMuonHit.dz();
        muon_station_hits->uncrossed[i] = commonMuonHit.uncrossed();
        muon_station_hits->time[i] = commonMuonHit.time();
        muon_station_hits->delta_time[i] = commonMuonHit.deltaTime();
        muon_station_hits->cluster_size[i] = commonMuonHit.clusterSize();
        muon_station_hits->region[i] = commonMuonHit.tile().region();
        i++;
      }
    }
    muon_station_hits->number_of_hits_per_station[currentStation] = hitsPerStation;
    muon_station_hits->station_offsets[currentStation] = offset;
    currentStation++;
    offset += hitsPerStation;
  }
}
*/