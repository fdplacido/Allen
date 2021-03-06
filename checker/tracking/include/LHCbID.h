/** @file LHCbID.h
 *
 * @brief encapsulate an LHCbID
 *
 * @author Rainer Schwemmer
 * @author Daniel Campora
 * @author Manuel Schiller
 * @date 2018-02-18
 */

#pragma once

#include <array>
#include <cstdint>
#include <cmath>
#include <vector>

/// encapsulate an LHCbID
class LHCbID {
private:
  uint32_t m_id;

public:
  constexpr LHCbID(uint32_t id) : m_id(id) {}

  enum channelIDtype { Velo = 0x8, FT = 0xa, UT = 0xb };

  LHCbID() = default;
  LHCbID(const LHCbID& other) = default;
  LHCbID(LHCbID&& other) = default;
  LHCbID& operator=(const LHCbID& other) = default;
  LHCbID& operator=(LHCbID&& other) = default;

  /// convert back to integer
  constexpr operator uint32_t() const noexcept { return m_id; }
  /// ordering of LHCbIDs
  bool operator==(const LHCbID& other) const noexcept
  {
    // dcampora: Relaxed check:
    // if (m_id == other.m_id) {
    //     return true;
    // } else if ((m_id & 0xFFFF0000) == (other.m_id & 0xFFFF0000)) {
    //     const uint row = m_id & 0xFF;
    //     const uint col = (m_id >> 8) & 0xFF;
    //     const uint other_row = other.m_id & 0xFF;
    //     const uint other_col = (other.m_id >> 8) & 0xFF;
    //     const float distance = (row - other_row) * (row - other_row) +
    //                            (col - other_col) * (col - other_col);
    //     if (std::sqrt(distance) < 2.0) {
    //       return true;
    //     }
    // }
    // return false;
    return m_id == other.m_id;
  }
  /// ordering of LHCbIDs
  bool operator!=(const LHCbID& other) const noexcept { return !this->operator==(other); }
  /// ordering of LHCbIDs
  constexpr bool operator<(const LHCbID& other) const noexcept { return m_id < other.m_id; }
  /// ordering of LHCbIDs
  constexpr bool operator<=(const LHCbID& other) const noexcept { return m_id <= other.m_id; }
  /// ordering of LHCbIDs
  constexpr bool operator>(const LHCbID& other) const noexcept { return m_id > other.m_id; }
  /// ordering of LHCbIDs
  constexpr bool operator>=(const LHCbID& other) const noexcept { return m_id >= other.m_id; }

  // get subdetector
  inline uint32_t detectorType() const { return (unsigned int) ((m_id & 0xF0000000L) >> 28); }
  bool isVelo() const { return (Velo == detectorType()); };
  bool isUT() const { return (UT == detectorType()); };
  bool isSciFi() const { return (FT == detectorType()); };
};

typedef std::vector<LHCbID> LHCbIDs;
