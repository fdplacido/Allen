#include "../../../main/include/Common.h"

#include "../../../main/include/Tools.h"


//#include "../../../PatPV/include/PVSeedTool.h"
#include "../../../PatPV/include/AdaptivePV3DFitter.h"
#include "../../../cuda/patPV/include/patPV_Definitions.cuh"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include <algorithm>



struct MCVertex {
  double x;
  double y;
  double z;
  int numberTracks;
};

typedef struct {  
  MCVertex* pMCPV;     // pointer to MC PV 
  int nRecTracks;            // number of reconstructed tracks from this MCPV  
  int nRecBackTracks;        // number of reconstructed backward tracks
  int indexRecPVInfo;        // index to reconstructed PVInfo (-1 if not reco)
  int nCorrectTracks;        // correct tracks belonging to reconstructed PV
  int multClosestMCPV;       // multiplicity of closest reconstructable MCPV
  double distToClosestMCPV;  // distance to closest reconstructible MCPV
  int decayCharm;             // type of mother particle                                                                                                                                           
  int decayBeauty;
  //std::vector<LHCb::MCParticle*> m_mcPartInMCPV;
  //std::vector<LHCb::Track*> m_recTracksInMCPV;
} MCPVInfo;

typedef struct {
 public:
  int nTracks;                    // number of tracks
  int nVeloTracks;                // number of velo tracks in a vertex       
  int nLongTracks;
  double minTrackRD;              //                                                                                                                                        
  double maxTrackRD;              //                                                                             
  double chi2;
  double nDoF;
  double d0;
  double d0nTr;
  double chi2nTr;
  double mind0;
  double maxd0;
  int mother;
  //XYZPoint position;       // position
  double x;
  double y;
  double z;
  XYZPoint positionSigma;  // position sigmas
  int indexMCPVInfo;              // index to MCPVInfo
  Vertex* pRECPV;        // pointer to REC PV
} RecPVInfo;


void match_mc_vertex_by_distance(int ipv, std::vector<RecPVInfo>& rinfo, std::vector<MCPVInfo>& mcpvvec) {

  double mindist = 999999.;
  int indexmc = -1;

  for(int imc = 0; imc < (int) mcpvvec.size(); imc++) {
    if ( mcpvvec[imc].indexRecPVInfo  > -1) continue;
    double dist = fabs(mcpvvec[imc].pMCPV->z -
                       rinfo[ipv].z);
    if(dist < mindist) {
      mindist = dist;
      indexmc = imc;
    }
  }
  if ( indexmc > -1 ) {
    if(mindist < 5.0 * rinfo[ipv].positionSigma.z) {
      rinfo[ipv].indexMCPVInfo = indexmc;
      mcpvvec[indexmc].indexRecPVInfo = ipv;
    }
  }

}

void printRat(std::string mes, int a, int b) {

  double rat = 0.;
  if(b>0) rat = 1.0*a/b;

  // reformat message
  unsigned int len = 20;
  std::string pmes = mes;
  while(pmes.length() < len) {
    pmes+=" ";
  }
  pmes+= " : ";

  std::cout << pmes << " " << rat << std::endl;

}