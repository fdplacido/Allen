#pragma once

#include <map>
#include <string>

#include <Dumpers/IUpdater.h>

int allen(std::map<std::string, std::string> options, Allen::NonEventData::IUpdater* updater);

namespace {
  constexpr size_t n_io = 1;
  constexpr size_t n_mon = 1;
  constexpr size_t max_stream_threads = 1024;
} // namespace
