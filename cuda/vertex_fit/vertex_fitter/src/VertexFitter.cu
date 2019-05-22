#include "VertexFitter.cuh"
#include "ParKalmanMath.cuh"
#include "ParKalmanDefinitions.cuh"
#include <gsl-lite.hpp>

namespace VertexFit {

  //----------------------------------------------------------------------
  // Point of closest approach. Reimplementation from TrackVertexUtils.
  __device__ bool poca(
    const ParKalmanFilter::FittedTrack& trackA,
    const ParKalmanFilter::FittedTrack& trackB,
    float& x,
    float& y,
    float& z)
  {
    float zA = trackA.z;
    float xA = trackA.state[0];
    float yA = trackA.state[1];
    float txA = trackA.state[2];
    float tyA = trackA.state[3];
    float zB = trackB.z;
    float xB = trackB.state[0];
    float yB = trackB.state[1];
    float txB = trackB.state[2];
    float tyB = trackB.state[3];
    float secondAA = txA * txA + tyA * tyA + 1.0;
    float secondBB = txB * txB + tyB * tyB + 1.0;
    float secondAB = -txA * txB - tyA * tyB - 1.0;
    float det = secondAA * secondBB - secondAB * secondAB;
    if (std::abs(det) > 0) {
      float secondinvAA = secondBB / det;
      float secondinvBB = secondAA / det;
      float secondinvAB = -secondAB / det;
      float firstA = txA * (xA - xB) + tyA * (yA - yB) + (zA - zB);
      float firstB = -txB * (xA - xB) - tyB * (yA - yB) - (zA - zB);
      float muA = -(secondinvAA * firstA + secondinvAB * firstB);
      float muB = -(secondinvBB * firstB + secondinvAB * firstA);
      x = 0.5 * (xA + muA * txA + xB + muB * txB);
      y = 0.5 * (yA + muA * tyA + yB + muB * tyB);
      z = 0.5 * (zA + muA + zB + muB);
      return true;
    }
    return false;
  }

  //----------------------------------------------------------------------
  // Add the contribution from one track to the vertex weight
  // matrix. NB this assumes (x, tx) and (y, ty) are uncorrelated.
  __device__ float addToDerivatives(
    const ParKalmanFilter::FittedTrack& track,
    const float& x,
    const float& y,
    const float& z,
    float& halfDChi2_0,
    float& halfDChi2_1,
    float& halfDChi2_2,
    float& halfD2Chi2_00,
    float& halfD2Chi2_11,
    float& halfD2Chi2_20,
    float& halfD2Chi2_21,
    float& halfD2Chi2_22)
  {
    float dz = z - track.z;
    float rX = track.state[0] + dz * track.state[2] - x;
    float rY = track.state[1] + dz * track.state[3] - y;
    float cov00 = track.cov(0, 0) + dz * dz * track.cov(2, 2) + 2 * dz * track.cov(2, 0);
    float cov11 = track.cov(1, 1) + dz * dz * track.cov(3, 3) + 2 * dz * track.cov(3, 1);
    float invcov00 = 1. / cov00;
    float invcov11 = 1. / cov11;
    halfDChi2_0 += invcov00 * rX;
    halfDChi2_1 += invcov11 * rY;
    halfDChi2_2 += -(invcov00 * rX * track.state[2] + invcov11 * rY * track.state[3]);
    halfD2Chi2_00 += invcov00;
    halfD2Chi2_11 += invcov11;
    halfD2Chi2_20 += -invcov00 * track.state[2];
    halfD2Chi2_21 += -invcov11 * track.state[3];
    halfD2Chi2_22 += invcov00 * track.state[0] * track.state[0] + invcov11 * track.state[1] * track.state[1];
    return invcov00 * rX * rX + invcov11 * rY * rY;
  }

  //----------------------------------------------------------------------
  // Correct the POCA to find the vertex position.
  __device__ float solve(
    float& x,
    float& y,
    float& z,
    float& cov00,
    float& cov11,
    float& cov20,
    float& cov21,
    float& cov22,
    const float& halfDChi2_0,
    const float& halfDChi2_1,
    const float& halfDChi2_2,
    const float& halfD2Chi2_00,
    const float& halfD2Chi2_11,
    const float& halfD2Chi2_20,
    const float& halfD2Chi2_21,
    const float& halfD2Chi2_22)
  {
    const float det = halfD2Chi2_00 * halfD2Chi2_11 * halfD2Chi2_22
      - halfD2Chi2_00 * halfD2Chi2_21 * halfD2Chi2_21
      - halfD2Chi2_11 * halfD2Chi2_20 * halfD2Chi2_20;
    const float invdet = 1. / det;
    cov00 = (halfD2Chi2_11 * halfD2Chi2_22 - halfD2Chi2_21 * halfD2Chi2_21) * invdet;
    cov11 = (halfD2Chi2_00 * halfD2Chi2_22 - halfD2Chi2_20 * halfD2Chi2_20) * invdet;
    cov20 = -halfD2Chi2_11 * halfD2Chi2_20 * invdet;
    cov21 = -halfD2Chi2_00 * halfD2Chi2_21 * invdet;
    cov22 = halfD2Chi2_00 * halfD2Chi2_11 * invdet;
    x += halfDChi2_0 * cov00 + halfDChi2_2 * cov20;
    y += halfDChi2_1 * cov11 + halfDChi2_2 * cov21;
    z += halfDChi2_0 * cov20 + halfDChi2_1 * cov21 + halfDChi2_2 * cov22;
    return -1 * (halfDChi2_0 * (halfDChi2_0 * cov00 + halfDChi2_2 * cov20)
                 + halfDChi2_1 * (halfDChi2_1 * cov11 + halfDChi2_2 * cov21)
                 + halfDChi2_2 * (halfDChi2_0 * cov20 + halfDChi2_1 * cov21 + halfDChi2_2 * cov22));
  }

  //----------------------------------------------------------------------
  // Perform a vertex fit assuming x and y are uncorrelated.
  __device__ bool doFit(
    const ParKalmanFilter::FittedTrack& trackA,
    const ParKalmanFilter::FittedTrack& trackB,
    Vertex& vertex)
  {
    if (!poca(trackA, trackB, vertex.x, vertex.y, vertex.z)) {
      return false;
    };
    float vertexweight00 = 0.f;
    float vertexweight11 = 0.f;
    float vertexweight20 = 0.f;
    float vertexweight21 = 0.f;
    float vertexweight22 = 0.f;
    float halfDChi2_0 = 0.f;
    float halfDChi2_1 = 0.f;
    float halfDChi2_2 = 0.f;
    vertex.chi2 = addToDerivatives(
      trackA,
      vertex.x,
      vertex.y,
      vertex.z,
      halfDChi2_0,
      halfDChi2_1,
      halfDChi2_2,
      vertexweight00,
      vertexweight11,
      vertexweight20,
      vertexweight21,
      vertexweight22);
    vertex.chi2 += addToDerivatives(
      trackB,
      vertex.x,
      vertex.y,
      vertex.z,
      halfDChi2_0,
      halfDChi2_1,
      halfDChi2_2,
      vertexweight00,
      vertexweight11,
      vertexweight20,
      vertexweight21,
      vertexweight22);
    vertex.chi2 += solve(
      vertex.x,
      vertex.y,
      vertex.z,
      vertex.cov00,
      vertex.cov11,
      vertex.cov20,
      vertex.cov21,
      vertex.cov22,
      halfDChi2_0,
      halfDChi2_1,
      halfDChi2_2,
      vertexweight00,
      vertexweight11,
      vertexweight20,
      vertexweight21,
      vertexweight22);
    return true;
  }

  __device__ void fill_extra_info(
    Vertex& sv,
    const PV::Vertex& pv,
    const ParKalmanFilter::FittedTrack& trackA,
    const ParKalmanFilter::FittedTrack& trackB)
  {
    // SV momentum.
    sv.px = trackA.px() + trackB.px();
    sv.py = trackA.py() + trackB.py();
    sv.pz = trackA.pz() + trackB.pz();

    // Sum of track pT.
    sv.sumpt = trackA.pt() + trackB.pt();

    // Number of tracks with ip chi2 < 16.
    sv.ntrks16 = (trackA.ipChi2 < 16) + (trackB.ipChi2 < 16);
    
    // Get PV-SV separation.
    const float dx = sv.x - pv.position.x;
    const float dy = sv.y - pv.position.y;
    const float dz = sv.z - pv.position.z;
    const float fd = std::sqrt(dx * dx + dy * dy + dz * dz);

    // Get covariance and FD chi2.
    const float cov00 = sv.cov00 + pv.cov00;
    const float cov11 = sv.cov11 + pv.cov11;
    const float cov20 = sv.cov20 + pv.cov20;
    const float cov21 = sv.cov21 + pv.cov21;
    const float cov22 = sv.cov22 + pv.cov22;
    const float invdet = 1. / (cov00 * cov11 * cov22
                               - cov00 * cov21 * cov21
                               - cov11 * cov20 * cov20);
    const float invcov00 = (cov11 * cov22 - cov21 * cov21) * invdet;
    const float invcov11 = (cov00 * cov22 - cov20 * cov20) * invdet;
    const float invcov20 = -cov11 * cov20 * invdet;
    const float invcov21 = cov00 * cov21 * invdet;
    const float invcov22 = cov00 * cov22 * invdet;
    sv.fdchi2 = invcov00 * dx * dx + invcov11 * dy * dy + invcov22 * dz * dz
      + 2 * invcov20 * dx * dz + 2 * invcov21 * dy * dz;

    // PV-SV eta.
    sv.eta = std::atanh(dz / fd);
    
    // Corrected mass.
    const float mpi = 139.57;
    const float px = trackA.px() + trackB.px();
    const float py = trackA.py() + trackB.py();
    const float pz = trackA.pz() + trackB.pz();
    const float mvis2 = 2 * mpi * mpi
      + 2 * (std::sqrt((trackA.p() * trackA.p() + mpi * mpi) *
                       (trackB.p() * trackB.p() + mpi * mpi))
             - trackA.px() * trackB.px()
             - trackA.py() * trackB.py()
             - trackA.pz() * trackB.pz());
    const float pperp2 = ((py * dz - dy * pz) * (py * dz - dy * pz) +
                          (pz * dx - dz * px) * (pz * dx - dz * px) +
                          (px * dy - dx * py) * (px * dy - dx * py)) / fd / fd;
    sv.mcor = std::sqrt(mvis2 + pperp2) + std::sqrt(pperp2);
  }
  
}

__global__ void fit_secondary_vertices(
  const ParKalmanFilter::FittedTrack* dev_kf_tracks,
  int* dev_n_scifi_tracks,
  uint* dev_scifi_track_hit_number,
  char* dev_scifi_consolidated_hits,
  float* dev_scifi_qop,
  MiniState* dev_scifi_states,
  uint* dev_ut_indices,
  PV::Vertex* dev_multi_fit_vertices,
  uint* dev_number_of_multi_fit_vertices,
  char* dev_kalman_pv_ipchi2,
  uint* dev_sv_offsets,
  VertexFit::Vertex* dev_secondary_vertices)
{
  const uint number_of_events = gridDim.x;
  const uint event_number = blockIdx.x;
  const uint sv_offset = dev_sv_offsets[event_number];
  
  // Consolidated SciFi tracks.
  const SciFi::Consolidated::Tracks scifi_tracks {
    (uint*) dev_n_scifi_tracks,
      (uint*) dev_scifi_track_hit_number,
      (float*) dev_scifi_qop,
      (MiniState*) dev_scifi_states,
      (uint*) dev_ut_indices,
      event_number,
      number_of_events};
  const uint event_tracks_offset = scifi_tracks.tracks_offset(event_number);
  const uint n_scifi_tracks = scifi_tracks.number_of_tracks(event_number);
  
  // Track-PV association table.
  const Associate::Consolidated::Table kalman_pv_ipchi2 {dev_kalman_pv_ipchi2, scifi_tracks.total_number_of_tracks};
  const auto pv_table = kalman_pv_ipchi2.event_table(scifi_tracks, event_number);

  // Kalman fitted tracks.
  const ParKalmanFilter::FittedTrack* event_tracks = dev_kf_tracks + event_tracks_offset;

  // Primary vertices.
  gsl::span<PV::Vertex const> vertices {dev_multi_fit_vertices + event_number * PV::max_number_vertices,
      *(dev_number_of_multi_fit_vertices + event_number)};

  // Secondary vertices.
  VertexFit::Vertex* event_secondary_vertices = dev_secondary_vertices + sv_offset;
    
  // Loop over tracks.
  for (int i_track = threadIdx.x; i_track < n_scifi_tracks; i_track += blockDim.x) {

    // Set the fit status for all possible vertices.
    for (int j_track = i_track + 1; j_track < n_scifi_tracks; j_track += 1) {
      uint vertex_idx = (int)n_scifi_tracks * ((int)n_scifi_tracks - 3) / 2
        - ((int)n_scifi_tracks - 1 - i_track) * ((int)n_scifi_tracks - 2 - i_track) / 2 + j_track;
      event_secondary_vertices[vertex_idx].chi2 = -1;
    }
    const ParKalmanFilter::FittedTrack trackA = event_tracks[i_track];
    
    // Preselection on first track.
    if (trackA.pt() < VertexFit::trackMinPt || trackA.ipChi2 < VertexFit::trackMinIPChi2) continue;
    
    // Loop over second track.
    for (int j_track = i_track + 1; j_track < n_scifi_tracks; j_track += 1) {

      // Preselection on second track.
      const ParKalmanFilter::FittedTrack trackB = event_tracks[j_track];
      if (trackB.pt() < VertexFit::trackMinPt || trackB.ipChi2 < VertexFit::trackMinIPChi2) continue;
      
      // Only combine tracks from the same PV.
      // Note: all tracks with IP chi2 > 16 are assigned to vertex -1.
      if (pv_table.pv[i_track] != pv_table.pv[j_track] && (pv_table.pv[i_track] < 1000 || pv_table.pv[j_track] < 1000)) continue;
            
      const int vertex_idx = (int)n_scifi_tracks * ((int)n_scifi_tracks - 3) / 2
        - ((int)n_scifi_tracks - 1 - i_track) * ((int)n_scifi_tracks - 2 - i_track) / 2 + j_track;
      
      // Do the vertex fit.
      doFit(
        trackA,
        trackB,
        event_secondary_vertices[vertex_idx]);

      // Fill extra info.
      int ipv = pv_table.pv[i_track];
      if (ipv >= 1000) ipv -= 1000;
      auto pv = vertices[ipv];
      fill_extra_info(
        event_secondary_vertices[vertex_idx],
        pv,
        trackA,
        trackB);
    }
    
  }
  
}
