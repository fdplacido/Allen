#pragma once

#include <cmath>
#include <array>
#include <vector>
#include <algorithm>
#include <fstream>

#include <cassert>

#include "SciFiDefinitions.cuh"

// access hits from a layer
// first zone number: y < 0
// second zone number: y > 0
void get_offset_and_n_hits_for_layer(
  const int first_zone, 
  const SciFi::HitCount& scifi_hit_count, 
  const float y, 
  int& n_hits, 
  int& zone_offset) {

  assert( first_zone < SciFi::Constants::n_zones-1 );

  if ( y < 0 ) {
    n_hits = scifi_hit_count.zone_number_of_hits(first_zone);
    zone_offset = scifi_hit_count.zone_offset(first_zone);
  } else {
    n_hits = scifi_hit_count.zone_number_of_hits(first_zone+1);
    zone_offset = scifi_hit_count.zone_offset(first_zone+1);
  }
   
  
}
