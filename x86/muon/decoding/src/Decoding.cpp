#include "Decoding.h"
void MuonRawToHitsDecoder::f() {

}
/*
void MuonRawToHitsDecoder::decode(gsl::span<char> events, gsl::span<unsigned int> offsets, std::vector<Muon::HitsSoA>& muon_hits_events) {
  std::string file_name_muon_table = "../input/muon/muon_table.bin";
  int number_of_events_requested = 5;
  MuonTable pad = MuonTable();
  MuonTable stripX = MuonTable();
  MuonTable stripY = MuonTable();
  char muon_table_raw_input[1200000];
  memset(muon_table_raw_input, 0, sizeof(muon_table_raw_input));
  std::ifstream muon_table_file(file_name_muon_table, std::ios::binary);
  muon_table_file.read(muon_table_raw_input, sizeof(muon_table_raw_input));
  muon_table_file.close();
  MuonTableReader muonTableReader = MuonTableReader();
  muonTableReader.read(muon_table_raw_input, &pad, &stripX, &stripY);
  MuonRawToHits muonRawToHits = MuonRawToHits(&pad, &stripX, &stripY);
  for (size_t i = 0; i < muon_hits_events.size(); i++) {
    gsl::span<char> rawEventSpan = {events.begin() + offsets[i], offsets[i + 1] - offsets[i]};
    Muon::MuonRawEvent rawEvent = Muon::MuonRawEvent((const char *) rawEventSpan.begin());
    muonRawToHits(rawEvent, &muon_hits_events[i]);
  }
}
*/