#pragma once

#include <stdint.h>
#include <ostream>

#include "SciFiDefinitions.cuh"

namespace SciFi {

  //----------------------------------------------------------------------
  // Struct for hit information. For now don't use a HitBase.
  struct Hit {
    float x0;
    float z0;
    float endPointY;
    uint32_t channel;
    
    // Cluster reference:
    //   raw bank: 8 bits
    //   element (it): 8 bits
    //   Condition 1-2-3: 2 bits
    //   Condition 2.1-2.2: 1 bit
    //   Condition 2.1: log2(n+1) - 8 bits
    uint32_t assembled_datatype;
    
    friend std::ostream& operator<<(std::ostream& stream, const Hit& hit) {
      stream << "SciFi hit {"
             << hit.x0 << ", "
             << hit.z0 << ", "
             << hit.channel << ", "
             << hit.assembled_datatype
             << "}";
      
      return stream;
    } 
  };

  
  /**
   * @brief Offset and number of hits of each layer.
   */
  struct HitCount {
    uint* mat_offsets;
    uint* n_hits_mats; 

    __device__ __host__
    void typecast_before_prefix_sum(
      uint* base_pointer,
      const uint event_number
    ){
      n_hits_mats = base_pointer + event_number * SciFi::Constants::n_mats;
    }
    
    __device__ __host__
    void typecast_after_prefix_sum(
      uint* base_pointer,
      const uint event_number,
      const uint number_of_events
    ){
      mat_offsets = base_pointer + event_number * SciFi::Constants::n_mats;
      n_hits_mats = base_pointer + number_of_events * SciFi::Constants::n_mats + 1 + event_number * SciFi::Constants::n_mats;
    }
    
    __device__ __host__
    uint mat_offset(const uint mat_number) const {
      return mat_offsets[mat_number];
    }
    
    __device__ __host__
    uint mat_number_of_hits(const uint mat_number) const {
      return mat_offsets[mat_number+1] - mat_offsets[mat_number];
    }

    __device__ __host__
    uint zone_offset(const uint zone_number) const {
      // TODO: Make this a constant
      constexpr uint32_t first_corrected_unique_mat_in_zone[] = {
      0, 40, 80, 120, 160, 200, 240, 280, 320, 360, 400, 440, 480, 520, 560, 600, 640,
      688, 736, 784, 832, 880, 928, 976, 1024};
      return mat_offsets[first_corrected_unique_mat_in_zone[zone_number]];
    }

    __device__ __host__
    uint zone_number_of_hits(const uint zone_number) const {
      return zone_offset(zone_number + 1) - zone_offset(zone_number);
    }
    
    __device__ __host__
    uint event_number_of_hits() const {
      return mat_offsets[SciFi::Constants::n_mats] - mat_offsets[0];
    }
    
    __device__ __host__
    uint event_offset() const {
      return mat_offsets[0];
    } 
  };

  struct BaseHits {
    float* x0;
    float* z0;
    float* m_endPointY;
    uint32_t* channel; 
    uint32_t* assembled_datatype;   
    
    const SciFiGeometry* geom;
    const float *dev_inv_clus_res;
    
    __device__ __host__ float w(uint32_t index) const {
      assert(pseudoSize(index) < 9 && "Wrong pseudo size.");
      float werrX = dev_inv_clus_res[pseudoSize(index)];
      return werrX * werrX;
    };
    
    __device__ __host__ float dxdy(uint32_t index) const {
      return geom->dxdy[mat(index)];
    };
    
    __device__ __host__ float dzdy(uint32_t index) const {
      return geom->dzdy[mat(index)];
    };

    __device__ __host__ float endPointY(uint32_t index) const {
      const SciFiChannelID id(channel[index]);
      float uFromChannel = geom->uBegin[mat(index)] + (2 * id.channel() + 1 + fraction(index)) * geom->halfChannelPitch[mat(index)];
      if( id.die() ) uFromChannel += geom->dieGap[mat(index)];
      return geom->mirrorPointY[mat(index)] + geom->ddxY[mat(index)] * uFromChannel;
    }

    __device__ __host__ float yMin(uint32_t index) const {
      const SciFiChannelID id(channel[index]);
      return m_endPointY[index] + id.isBottom() * geom->globaldy[mat(index)];
    };

    __device__ __host__ float yMax(uint32_t index) const {
      const SciFiChannelID id(channel[index]);
      return m_endPointY[index] + !id.isBottom() * geom->globaldy[mat(index)];
    };

    __device__ __host__ uint32_t LHCbID(uint32_t index) const {
      return (10u << 28) + channel[index];
    }; 

    __device__ __host__ uint32_t mat(uint32_t index) const {
      return assembled_datatype[index] & 0x7ff;
    };

    __device__ __host__ uint32_t pseudoSize(uint32_t index) const{
      return (assembled_datatype[index] >> 11) & 0xf;
    };

    __device__ __host__ uint32_t planeCode(uint32_t index) const {
      return (assembled_datatype[index] >> 15) & 0x1f;
    };

    __device__ __host__ uint32_t fraction(uint32_t index) const{
      return (assembled_datatype[index] >> 20) & 0x1;
    };
  
  }; 

  struct Hits : BaseHits {
    // Cluster reference:
    //   raw bank: 8 bits
    //   element (it): 8 bits
    //   Condition 1-2-3: 2 bits
    //   Condition 2.1-2.2: 1 bit
    //   Condition 2.1: log2(n+1) - 8 bits
    uint32_t* cluster_reference;   
    
    __device__ __host__
    Hits(uint* base,
      const uint32_t total_number_of_hits,
      const SciFiGeometry* param_geom,
      const float* param_dev_inv_clus_res) 
    {
      geom = param_geom;
      x0 = reinterpret_cast<float*>(base);
      z0 = reinterpret_cast<float*>(base + total_number_of_hits);
      m_endPointY = reinterpret_cast<float*>(base + 2 * total_number_of_hits);
      channel = reinterpret_cast<uint32_t*>(base + 3 * total_number_of_hits);
      assembled_datatype = reinterpret_cast<uint32_t*>(base + 4 * total_number_of_hits);
      cluster_reference = reinterpret_cast<uint32_t*>(base + 5 * total_number_of_hits);
      dev_inv_clus_res = param_dev_inv_clus_res; 
    } 
  };

  /**
   * Track object used for storing tracks
   */
  struct TrackHits {
    short hits[SciFi::Constants::max_track_size];
    float qop;
    unsigned short hitsNum = 0;
    float chi2;
    unsigned int UTTrackIndex; // Index of velo-UT track

    __host__ __device__ void addHit(unsigned int idx){
      assert(hitsNum < SciFi::Constants::max_track_size);
      hits[hitsNum++] = idx;
    }
  };
  
}