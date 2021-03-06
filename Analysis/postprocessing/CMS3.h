// -*- C++ -*-
#ifndef CMS3_H
#define CMS3_H
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TBits.h"
#include <vector>
#include <unistd.h>
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

// Generated with file: output_12.root

using namespace std;
class CMS3 {
private:
protected:
  unsigned int index;
  float ele_full5x5_hcalOverEcal_;
  TBranch *ele_full5x5_hcalOverEcal_branch;
  bool ele_full5x5_hcalOverEcal_isLoaded;
  int MC_TrueNumInteractions_;
  TBranch *MC_TrueNumInteractions_branch;
  bool MC_TrueNumInteractions_isLoaded;
  float ele_fbrem_;
  TBranch *ele_fbrem_branch;
  bool ele_fbrem_isLoaded;
  float ele_oldsirir_;
  TBranch *ele_oldsirir_branch;
  bool ele_oldsirir_isLoaded;
  int ele_N_;
  TBranch *ele_N_branch;
  bool ele_N_isLoaded;
  float mc_gen_pt_;
  TBranch *mc_gen_pt_branch;
  bool mc_gen_pt_isLoaded;
  float mc_event_weight_;
  TBranch *mc_event_weight_branch;
  bool mc_event_weight_isLoaded;
  int ele_kfhits_;
  TBranch *ele_kfhits_branch;
  bool ele_kfhits_isLoaded;
  bool ele_trackerDrivenSeed_;
  TBranch *ele_trackerDrivenSeed_branch;
  bool ele_trackerDrivenSeed_isLoaded;
  int mc_ele_matchedFromCB2_;
  TBranch *mc_ele_matchedFromCB2_branch;
  bool mc_ele_matchedFromCB2_isLoaded;
  vector<float> *rhs_e_;
  TBranch *rhs_e_branch;
  bool rhs_e_isLoaded;
  float ele_eelepout_;
  TBranch *ele_eelepout_branch;
  bool ele_eelepout_isLoaded;
  bool ele_isEBEEGap_;
  TBranch *ele_isEBEEGap_branch;
  bool ele_isEBEEGap_isLoaded;
  float ele_oldsigmaietaieta_;
  TBranch *ele_oldsigmaietaieta_branch;
  bool ele_oldsigmaietaieta_isLoaded;
  int nEvent_;
  TBranch *nEvent_branch;
  bool nEvent_isLoaded;
  int ele_nbrem_;
  TBranch *ele_nbrem_branch;
  bool ele_nbrem_isLoaded;
  int ele_index_;
  TBranch *ele_index_branch;
  bool ele_index_isLoaded;
  float scl_phi_;
  TBranch *scl_phi_branch;
  bool scl_phi_isLoaded;
  bool ele_isEBPhiGap_;
  TBranch *ele_isEBPhiGap_branch;
  bool ele_isEBPhiGap_isLoaded;
  bool ele_isEE_;
  TBranch *ele_isEE_branch;
  bool ele_isEE_isLoaded;
  bool ele_isEB_;
  TBranch *ele_isEB_branch;
  bool ele_isEB_isLoaded;
  float ele_SIP_;
  TBranch *ele_SIP_branch;
  bool ele_SIP_isLoaded;
  float ele_sclphiwidth_;
  TBranch *ele_sclphiwidth_branch;
  bool ele_sclphiwidth_isLoaded;
  int mc_ele_matchedFromCB_;
  TBranch *mc_ele_matchedFromCB_branch;
  bool mc_ele_matchedFromCB_isLoaded;
  int ele_ID2_cat_;
  TBranch *ele_ID2_cat_branch;
  bool ele_ID2_cat_isLoaded;
  float scl_E_;
  TBranch *scl_E_branch;
  bool scl_E_isLoaded;
  int ele_eClass_;
  TBranch *ele_eClass_branch;
  bool ele_eClass_isLoaded;
  float mc_gen_eta_;
  TBranch *mc_gen_eta_branch;
  bool mc_gen_eta_isLoaded;
  float ele_sclRawE_;
  TBranch *ele_sclRawE_branch;
  bool ele_sclRawE_isLoaded;
  int PU_N_;
  TBranch *PU_N_branch;
  bool PU_N_isLoaded;
  float ele_kfchi2_;
  TBranch *ele_kfchi2_branch;
  bool ele_kfchi2_isLoaded;
  float ele_hadronicOverEm_;
  TBranch *ele_hadronicOverEm_branch;
  bool ele_hadronicOverEm_isLoaded;
  bool is_signal_;
  TBranch *is_signal_branch;
  bool is_signal_isLoaded;
  float ele_gsfchi2_;
  TBranch *ele_gsfchi2_branch;
  bool ele_gsfchi2_isLoaded;
  float ele_oldcircularity_;
  TBranch *ele_oldcircularity_branch;
  bool ele_oldcircularity_isLoaded;
  vector<int> *rhs_iphi_;
  TBranch *rhs_iphi_branch;
  bool rhs_iphi_isLoaded;
  int ele_sclNclus_;
  TBranch *ele_sclNclus_branch;
  bool ele_sclNclus_isLoaded;
  float ele_IoEmIop_;
  TBranch *ele_IoEmIop_branch;
  bool ele_IoEmIop_isLoaded;
  float ele_pt_;
  TBranch *ele_pt_branch;
  bool ele_pt_isLoaded;
  int seed_iphi_;
  TBranch *seed_iphi_branch;
  bool seed_iphi_isLoaded;
  int vtx_N_;
  TBranch *vtx_N_branch;
  bool vtx_N_isLoaded;
  float ele_ID2_;
  TBranch *ele_ID2_branch;
  bool ele_ID2_isLoaded;
  int ele_valid_hits_;
  TBranch *ele_valid_hits_branch;
  bool ele_valid_hits_isLoaded;
  float ele_olde55_;
  TBranch *ele_olde55_branch;
  bool ele_olde55_isLoaded;
  float ele_oldr9_;
  TBranch *ele_oldr9_branch;
  bool ele_oldr9_isLoaded;
  bool ele_isEERingGap_;
  TBranch *ele_isEERingGap_branch;
  bool ele_isEERingGap_isLoaded;
  float ele_deltaetain_;
  TBranch *ele_deltaetain_branch;
  bool ele_deltaetain_isLoaded;
  int ele_lost_hits_;
  TBranch *ele_lost_hits_branch;
  bool ele_lost_hits_isLoaded;
  int mc_pho_mother_status_;
  TBranch *mc_pho_mother_status_branch;
  bool mc_pho_mother_status_isLoaded;
  float ele_oldhebc_;
  TBranch *ele_oldhebc_branch;
  bool ele_oldhebc_isLoaded;
  float ele_conv_radius_;
  TBranch *ele_conv_radius_branch;
  bool ele_conv_radius_isLoaded;
  float ele_scletawidth_;
  TBranch *ele_scletawidth_branch;
  bool ele_scletawidth_isLoaded;
  float ele_ID1_;
  TBranch *ele_ID1_branch;
  bool ele_ID1_isLoaded;
  float ele_conv_dcot_;
  TBranch *ele_conv_dcot_branch;
  bool ele_conv_dcot_isLoaded;
  int mc_pho_mother_id_;
  TBranch *mc_pho_mother_id_branch;
  bool mc_pho_mother_id_isLoaded;
  vector<int> *rhs_ieta_;
  TBranch *rhs_ieta_branch;
  bool rhs_ieta_isLoaded;
  float scl_eta_;
  TBranch *scl_eta_branch;
  bool scl_eta_isLoaded;
  int nRun_;
  TBranch *nRun_branch;
  bool nRun_isLoaded;
  bool ele_isEBEtaGap_;
  TBranch *ele_isEBEtaGap_branch;
  bool ele_isEBEtaGap_isLoaded;
  int seed_ieta_;
  TBranch *seed_ieta_branch;
  bool seed_ieta_isLoaded;
  int ele_gsfhits_;
  TBranch *ele_gsfhits_branch;
  bool ele_gsfhits_isLoaded;
  float mc_gen_ID_;
  TBranch *mc_gen_ID_branch;
  bool mc_gen_ID_isLoaded;
  float ele_3x3_;
  TBranch *ele_3x3_branch;
  bool ele_3x3_isLoaded;
  float ele_psEoverEraw_;
  TBranch *ele_psEoverEraw_branch;
  bool ele_psEoverEraw_isLoaded;
  int ele_ID2_pass_;
  TBranch *ele_ID2_pass_branch;
  bool ele_ID2_pass_isLoaded;
  float scl_Et_;
  TBranch *scl_Et_branch;
  bool scl_Et_isLoaded;
  float ele_combErr_;
  TBranch *ele_combErr_branch;
  bool ele_combErr_isLoaded;
  float ele_deltaphiin_;
  TBranch *ele_deltaphiin_branch;
  bool ele_deltaphiin_isLoaded;
  float seed_e_;
  TBranch *seed_e_branch;
  bool seed_e_isLoaded;
  float ele_trackErr_;
  TBranch *ele_trackErr_branch;
  bool ele_trackErr_isLoaded;
  float ele_PFcombErr_;
  TBranch *ele_PFcombErr_branch;
  bool ele_PFcombErr_isLoaded;
  int ele_ID1_cat_;
  TBranch *ele_ID1_cat_branch;
  bool ele_ID1_cat_isLoaded;
  float ele_oldsigmaiphiiphi_;
  TBranch *ele_oldsigmaiphiiphi_branch;
  bool ele_oldsigmaiphiiphi_isLoaded;
  float ele_SCfbrem_;
  TBranch *ele_SCfbrem_branch;
  bool ele_SCfbrem_isLoaded;
  int ele_vtxconv_;
  TBranch *ele_vtxconv_branch;
  bool ele_vtxconv_isLoaded;
  float ele_dz_;
  TBranch *ele_dz_branch;
  bool ele_dz_isLoaded;
  int nLumi_;
  TBranch *nLumi_branch;
  bool nLumi_isLoaded;
  int ele_echarge_;
  TBranch *ele_echarge_branch;
  bool ele_echarge_isLoaded;
  int ele_ID1_pass_;
  TBranch *ele_ID1_pass_branch;
  bool ele_ID1_pass_isLoaded;
  bool ele_isEEDeeGap_;
  TBranch *ele_isEEDeeGap_branch;
  bool ele_isEEDeeGap_isLoaded;
  float ele_conversionVertexFitProbability_;
  TBranch *ele_conversionVertexFitProbability_branch;
  bool ele_conversionVertexFitProbability_isLoaded;
  float ele_olde25max_;
  TBranch *ele_olde25max_branch;
  bool ele_olde25max_isLoaded;
  bool ele_ecalDrivenSeed_;
  TBranch *ele_ecalDrivenSeed_branch;
  bool ele_ecalDrivenSeed_isLoaded;
  float ele_eta_;
  TBranch *ele_eta_branch;
  bool ele_eta_isLoaded;
  float ele_oldhe_;
  TBranch *ele_oldhe_branch;
  bool ele_oldhe_isLoaded;
  float ele_olde15_;
  TBranch *ele_olde15_branch;
  bool ele_olde15_isLoaded;
  int mc_ele_matchMother_PDGID_;
  TBranch *mc_ele_matchMother_PDGID_branch;
  bool mc_ele_matchMother_PDGID_isLoaded;
  float ele_ep_;
  TBranch *ele_ep_branch;
  bool ele_ep_isLoaded;
  float ele_deltaetaseed_;
  TBranch *ele_deltaetaseed_branch;
  bool ele_deltaetaseed_isLoaded;
  float ele_conv_dist_;
  TBranch *ele_conv_dist_branch;
  bool ele_conv_dist_isLoaded;
  int ele_expected_inner_hits_;
  TBranch *ele_expected_inner_hits_branch;
  bool ele_expected_inner_hits_isLoaded;
  int mc_ele_matchMother_status_;
  TBranch *mc_ele_matchMother_status_branch;
  bool mc_ele_matchMother_status_isLoaded;
public:
  void Init(TTree *tree);
  void GetEntry(unsigned int idx);
  void LoadAllBranches();
  const float &ele_full5x5_hcalOverEcal();
  const int &MC_TrueNumInteractions();
  const float &ele_fbrem();
  const float &ele_oldsirir();
  const int &ele_N();
  const float &mc_gen_pt();
  const float &mc_event_weight();
  const int &ele_kfhits();
  const bool &ele_trackerDrivenSeed();
  const int &mc_ele_matchedFromCB2();
  const vector<float> &rhs_e();
  const float &ele_eelepout();
  const bool &ele_isEBEEGap();
  const float &ele_oldsigmaietaieta();
  const int &nEvent();
  const int &ele_nbrem();
  const int &ele_index();
  const float &scl_phi();
  const bool &ele_isEBPhiGap();
  const bool &ele_isEE();
  const bool &ele_isEB();
  const float &ele_SIP();
  const float &ele_sclphiwidth();
  const int &mc_ele_matchedFromCB();
  const int &ele_ID2_cat();
  const float &scl_E();
  const int &ele_eClass();
  const float &mc_gen_eta();
  const float &ele_sclRawE();
  const int &PU_N();
  const float &ele_kfchi2();
  const float &ele_hadronicOverEm();
  const bool &is_signal();
  const float &ele_gsfchi2();
  const float &ele_oldcircularity();
  const vector<int> &rhs_iphi();
  const int &ele_sclNclus();
  const float &ele_IoEmIop();
  const float &ele_pt();
  const int &seed_iphi();
  const int &vtx_N();
  const float &ele_ID2();
  const int &ele_valid_hits();
  const float &ele_olde55();
  const float &ele_oldr9();
  const bool &ele_isEERingGap();
  const float &ele_deltaetain();
  const int &ele_lost_hits();
  const int &mc_pho_mother_status();
  const float &ele_oldhebc();
  const float &ele_conv_radius();
  const float &ele_scletawidth();
  const float &ele_ID1();
  const float &ele_conv_dcot();
  const int &mc_pho_mother_id();
  const vector<int> &rhs_ieta();
  const float &scl_eta();
  const int &nRun();
  const bool &ele_isEBEtaGap();
  const int &seed_ieta();
  const int &ele_gsfhits();
  const float &mc_gen_ID();
  const float &ele_3x3();
  const float &ele_psEoverEraw();
  const int &ele_ID2_pass();
  const float &scl_Et();
  const float &ele_combErr();
  const float &ele_deltaphiin();
  const float &seed_e();
  const float &ele_trackErr();
  const float &ele_PFcombErr();
  const int &ele_ID1_cat();
  const float &ele_oldsigmaiphiiphi();
  const float &ele_SCfbrem();
  const int &ele_vtxconv();
  const float &ele_dz();
  const int &nLumi();
  const int &ele_echarge();
  const int &ele_ID1_pass();
  const bool &ele_isEEDeeGap();
  const float &ele_conversionVertexFitProbability();
  const float &ele_olde25max();
  const bool &ele_ecalDrivenSeed();
  const float &ele_eta();
  const float &ele_oldhe();
  const float &ele_olde15();
  const int &mc_ele_matchMother_PDGID();
  const float &ele_ep();
  const float &ele_deltaetaseed();
  const float &ele_conv_dist();
  const int &ele_expected_inner_hits();
  const int &mc_ele_matchMother_status();
  static void progress( int nEventsTotal, int nEventsChain );
};

#ifndef __CINT__
extern CMS3 cms3;
#endif

namespace tas {

  const float &ele_full5x5_hcalOverEcal();
  const int &MC_TrueNumInteractions();
  const float &ele_fbrem();
  const float &ele_oldsirir();
  const int &ele_N();
  const float &mc_gen_pt();
  const float &mc_event_weight();
  const int &ele_kfhits();
  const bool &ele_trackerDrivenSeed();
  const int &mc_ele_matchedFromCB2();
  const vector<float> &rhs_e();
  const float &ele_eelepout();
  const bool &ele_isEBEEGap();
  const float &ele_oldsigmaietaieta();
  const int &nEvent();
  const int &ele_nbrem();
  const int &ele_index();
  const float &scl_phi();
  const bool &ele_isEBPhiGap();
  const bool &ele_isEE();
  const bool &ele_isEB();
  const float &ele_SIP();
  const float &ele_sclphiwidth();
  const int &mc_ele_matchedFromCB();
  const int &ele_ID2_cat();
  const float &scl_E();
  const int &ele_eClass();
  const float &mc_gen_eta();
  const float &ele_sclRawE();
  const int &PU_N();
  const float &ele_kfchi2();
  const float &ele_hadronicOverEm();
  const bool &is_signal();
  const float &ele_gsfchi2();
  const float &ele_oldcircularity();
  const vector<int> &rhs_iphi();
  const int &ele_sclNclus();
  const float &ele_IoEmIop();
  const float &ele_pt();
  const int &seed_iphi();
  const int &vtx_N();
  const float &ele_ID2();
  const int &ele_valid_hits();
  const float &ele_olde55();
  const float &ele_oldr9();
  const bool &ele_isEERingGap();
  const float &ele_deltaetain();
  const int &ele_lost_hits();
  const int &mc_pho_mother_status();
  const float &ele_oldhebc();
  const float &ele_conv_radius();
  const float &ele_scletawidth();
  const float &ele_ID1();
  const float &ele_conv_dcot();
  const int &mc_pho_mother_id();
  const vector<int> &rhs_ieta();
  const float &scl_eta();
  const int &nRun();
  const bool &ele_isEBEtaGap();
  const int &seed_ieta();
  const int &ele_gsfhits();
  const float &mc_gen_ID();
  const float &ele_3x3();
  const float &ele_psEoverEraw();
  const int &ele_ID2_pass();
  const float &scl_Et();
  const float &ele_combErr();
  const float &ele_deltaphiin();
  const float &seed_e();
  const float &ele_trackErr();
  const float &ele_PFcombErr();
  const int &ele_ID1_cat();
  const float &ele_oldsigmaiphiiphi();
  const float &ele_SCfbrem();
  const int &ele_vtxconv();
  const float &ele_dz();
  const int &nLumi();
  const int &ele_echarge();
  const int &ele_ID1_pass();
  const bool &ele_isEEDeeGap();
  const float &ele_conversionVertexFitProbability();
  const float &ele_olde25max();
  const bool &ele_ecalDrivenSeed();
  const float &ele_eta();
  const float &ele_oldhe();
  const float &ele_olde15();
  const int &mc_ele_matchMother_PDGID();
  const float &ele_ep();
  const float &ele_deltaetaseed();
  const float &ele_conv_dist();
  const int &ele_expected_inner_hits();
  const int &mc_ele_matchMother_status();
}
#endif
