#include <string>

#include <CudaCommon.h>
#include <Common.h>
#include <Consumers.h>

namespace {
  using std::string;
  using std::to_string;
} // namespace

Consumers::Beamline::Beamline(gsl::span<float>& dev_beamline) : m_dev_beamline {dev_beamline} {}

void Consumers::Beamline::consume(std::vector<char> const& data)
{
  if (m_dev_beamline.get().empty()) {
    // Allocate space
    float* p = nullptr;
    cudaCheck(cudaMalloc((void**) &p, data.size()));
    m_dev_beamline.get() = {p, data.size() / sizeof(float)};
  }
  else if (data.size() != m_dev_beamline.get().size()) {
    throw StrException {string {"sizes don't match: "} + to_string(m_dev_beamline.get().size()) + " " +
                        to_string(data.size())};
  }

  cudaCheck(cudaMemcpy(m_dev_beamline.get().data(), data.data(), data.size(), cudaMemcpyHostToDevice));
}
