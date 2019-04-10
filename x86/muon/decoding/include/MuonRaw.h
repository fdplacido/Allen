#pragma once

#include <stdint.h>
#include <gsl-lite.hpp>

namespace Muon {
  struct MuonRawBank {
    uint32_t sourceID;
    uint16_t* data;
    uint16_t* last;
    size_t m_size;

    MuonRawBank(const char* raw_bank, const char* end) {
      const char* p = raw_bank;
      sourceID = *((uint32_t*) p);
      p += sizeof(uint32_t);
      data = (uint16_t*) p;
      last = (uint16_t*) end;
      m_size = last - data;
    }

    size_t size() const {
      return m_size;
    }

    /// Begin iterator
    template <typename T> T* begin(){    return (T*)data;         }

    /// End iterator
    template <typename T> T* end()  {    return ((T*)data) + size()/sizeof(T);  }

    /// Begin iterator over const iteration
    template <typename T> const T* begin() const {return (T*)data;}

    /// End iterator of const iteration
    template <typename T>  const T* end() const  {return ((T*)data) + size()/sizeof(T);  }

    /// return a range of 'const T'
    template <typename T>
    gsl::span<const T> range() const {
      return {begin<T>(), end<T>()};
    }
  };

  struct MuonRawEvent {
    uint32_t number_of_raw_banks;
    uint32_t* raw_bank_offset;
    char* payload;

    MuonRawEvent(const char* event) {
      const char* p = event;
      number_of_raw_banks = *((uint32_t*) p);
      p += sizeof(uint32_t);
      raw_bank_offset = (uint32_t*) p;
      p += (number_of_raw_banks + 1) * sizeof(uint32_t);
      payload = (char*) p;
    }

    MuonRawBank getMuonBank(const uint32_t index) const {
      MuonRawBank bank(payload + raw_bank_offset[index], payload + raw_bank_offset[index + 1]);
      return bank;
    }
  };
}
