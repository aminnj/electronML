#include "CMS3.h"
CMS3 cms3;

void CMS3::Init(TTree *tree) {
  tree->SetMakeClass(1);
  ele_full5x5_hcalOverEcal_branch = 0;
  if (tree->GetBranch("ele_full5x5_hcalOverEcal") != 0) {
    ele_full5x5_hcalOverEcal_branch = tree->GetBranch("ele_full5x5_hcalOverEcal");
    if (ele_full5x5_hcalOverEcal_branch) { ele_full5x5_hcalOverEcal_branch->SetAddress(&ele_full5x5_hcalOverEcal_); }
  }
  MC_TrueNumInteractions_branch = 0;
  if (tree->GetBranch("MC_TrueNumInteractions") != 0) {
    MC_TrueNumInteractions_branch = tree->GetBranch("MC_TrueNumInteractions");
    if (MC_TrueNumInteractions_branch) { MC_TrueNumInteractions_branch->SetAddress(&MC_TrueNumInteractions_); }
  }
  ele_fbrem_branch = 0;
  if (tree->GetBranch("ele_fbrem") != 0) {
    ele_fbrem_branch = tree->GetBranch("ele_fbrem");
    if (ele_fbrem_branch) { ele_fbrem_branch->SetAddress(&ele_fbrem_); }
  }
  ele_oldsirir_branch = 0;
  if (tree->GetBranch("ele_oldsirir") != 0) {
    ele_oldsirir_branch = tree->GetBranch("ele_oldsirir");
    if (ele_oldsirir_branch) { ele_oldsirir_branch->SetAddress(&ele_oldsirir_); }
  }
  ele_N_branch = 0;
  if (tree->GetBranch("ele_N") != 0) {
    ele_N_branch = tree->GetBranch("ele_N");
    if (ele_N_branch) { ele_N_branch->SetAddress(&ele_N_); }
  }
  mc_gen_pt_branch = 0;
  if (tree->GetBranch("mc_gen_pt") != 0) {
    mc_gen_pt_branch = tree->GetBranch("mc_gen_pt");
    if (mc_gen_pt_branch) { mc_gen_pt_branch->SetAddress(&mc_gen_pt_); }
  }
  mc_event_weight_branch = 0;
  if (tree->GetBranch("mc_event_weight") != 0) {
    mc_event_weight_branch = tree->GetBranch("mc_event_weight");
    if (mc_event_weight_branch) { mc_event_weight_branch->SetAddress(&mc_event_weight_); }
  }
  ele_kfhits_branch = 0;
  if (tree->GetBranch("ele_kfhits") != 0) {
    ele_kfhits_branch = tree->GetBranch("ele_kfhits");
    if (ele_kfhits_branch) { ele_kfhits_branch->SetAddress(&ele_kfhits_); }
  }
  ele_trackerDrivenSeed_branch = 0;
  if (tree->GetBranch("ele_trackerDrivenSeed") != 0) {
    ele_trackerDrivenSeed_branch = tree->GetBranch("ele_trackerDrivenSeed");
    if (ele_trackerDrivenSeed_branch) { ele_trackerDrivenSeed_branch->SetAddress(&ele_trackerDrivenSeed_); }
  }
  mc_ele_matchedFromCB2_branch = 0;
  if (tree->GetBranch("mc_ele_matchedFromCB2") != 0) {
    mc_ele_matchedFromCB2_branch = tree->GetBranch("mc_ele_matchedFromCB2");
    if (mc_ele_matchedFromCB2_branch) { mc_ele_matchedFromCB2_branch->SetAddress(&mc_ele_matchedFromCB2_); }
  }
  rhs_e_branch = 0;
  if (tree->GetBranch("rhs_e") != 0) {
    rhs_e_branch = tree->GetBranch("rhs_e");
    if (rhs_e_branch) { rhs_e_branch->SetAddress(&rhs_e_); }
  }
  ele_eelepout_branch = 0;
  if (tree->GetBranch("ele_eelepout") != 0) {
    ele_eelepout_branch = tree->GetBranch("ele_eelepout");
    if (ele_eelepout_branch) { ele_eelepout_branch->SetAddress(&ele_eelepout_); }
  }
  ele_isEBEEGap_branch = 0;
  if (tree->GetBranch("ele_isEBEEGap") != 0) {
    ele_isEBEEGap_branch = tree->GetBranch("ele_isEBEEGap");
    if (ele_isEBEEGap_branch) { ele_isEBEEGap_branch->SetAddress(&ele_isEBEEGap_); }
  }
  ele_oldsigmaietaieta_branch = 0;
  if (tree->GetBranch("ele_oldsigmaietaieta") != 0) {
    ele_oldsigmaietaieta_branch = tree->GetBranch("ele_oldsigmaietaieta");
    if (ele_oldsigmaietaieta_branch) { ele_oldsigmaietaieta_branch->SetAddress(&ele_oldsigmaietaieta_); }
  }
  nEvent_branch = 0;
  if (tree->GetBranch("nEvent") != 0) {
    nEvent_branch = tree->GetBranch("nEvent");
    if (nEvent_branch) { nEvent_branch->SetAddress(&nEvent_); }
  }
  ele_nbrem_branch = 0;
  if (tree->GetBranch("ele_nbrem") != 0) {
    ele_nbrem_branch = tree->GetBranch("ele_nbrem");
    if (ele_nbrem_branch) { ele_nbrem_branch->SetAddress(&ele_nbrem_); }
  }
  ele_index_branch = 0;
  if (tree->GetBranch("ele_index") != 0) {
    ele_index_branch = tree->GetBranch("ele_index");
    if (ele_index_branch) { ele_index_branch->SetAddress(&ele_index_); }
  }
  scl_phi_branch = 0;
  if (tree->GetBranch("scl_phi") != 0) {
    scl_phi_branch = tree->GetBranch("scl_phi");
    if (scl_phi_branch) { scl_phi_branch->SetAddress(&scl_phi_); }
  }
  ele_isEBPhiGap_branch = 0;
  if (tree->GetBranch("ele_isEBPhiGap") != 0) {
    ele_isEBPhiGap_branch = tree->GetBranch("ele_isEBPhiGap");
    if (ele_isEBPhiGap_branch) { ele_isEBPhiGap_branch->SetAddress(&ele_isEBPhiGap_); }
  }
  ele_isEE_branch = 0;
  if (tree->GetBranch("ele_isEE") != 0) {
    ele_isEE_branch = tree->GetBranch("ele_isEE");
    if (ele_isEE_branch) { ele_isEE_branch->SetAddress(&ele_isEE_); }
  }
  ele_isEB_branch = 0;
  if (tree->GetBranch("ele_isEB") != 0) {
    ele_isEB_branch = tree->GetBranch("ele_isEB");
    if (ele_isEB_branch) { ele_isEB_branch->SetAddress(&ele_isEB_); }
  }
  ele_SIP_branch = 0;
  if (tree->GetBranch("ele_SIP") != 0) {
    ele_SIP_branch = tree->GetBranch("ele_SIP");
    if (ele_SIP_branch) { ele_SIP_branch->SetAddress(&ele_SIP_); }
  }
  ele_sclphiwidth_branch = 0;
  if (tree->GetBranch("ele_sclphiwidth") != 0) {
    ele_sclphiwidth_branch = tree->GetBranch("ele_sclphiwidth");
    if (ele_sclphiwidth_branch) { ele_sclphiwidth_branch->SetAddress(&ele_sclphiwidth_); }
  }
  mc_ele_matchedFromCB_branch = 0;
  if (tree->GetBranch("mc_ele_matchedFromCB") != 0) {
    mc_ele_matchedFromCB_branch = tree->GetBranch("mc_ele_matchedFromCB");
    if (mc_ele_matchedFromCB_branch) { mc_ele_matchedFromCB_branch->SetAddress(&mc_ele_matchedFromCB_); }
  }
  ele_ID2_cat_branch = 0;
  if (tree->GetBranch("ele_ID2_cat") != 0) {
    ele_ID2_cat_branch = tree->GetBranch("ele_ID2_cat");
    if (ele_ID2_cat_branch) { ele_ID2_cat_branch->SetAddress(&ele_ID2_cat_); }
  }
  scl_E_branch = 0;
  if (tree->GetBranch("scl_E") != 0) {
    scl_E_branch = tree->GetBranch("scl_E");
    if (scl_E_branch) { scl_E_branch->SetAddress(&scl_E_); }
  }
  ele_eClass_branch = 0;
  if (tree->GetBranch("ele_eClass") != 0) {
    ele_eClass_branch = tree->GetBranch("ele_eClass");
    if (ele_eClass_branch) { ele_eClass_branch->SetAddress(&ele_eClass_); }
  }
  mc_gen_eta_branch = 0;
  if (tree->GetBranch("mc_gen_eta") != 0) {
    mc_gen_eta_branch = tree->GetBranch("mc_gen_eta");
    if (mc_gen_eta_branch) { mc_gen_eta_branch->SetAddress(&mc_gen_eta_); }
  }
  ele_sclRawE_branch = 0;
  if (tree->GetBranch("ele_sclRawE") != 0) {
    ele_sclRawE_branch = tree->GetBranch("ele_sclRawE");
    if (ele_sclRawE_branch) { ele_sclRawE_branch->SetAddress(&ele_sclRawE_); }
  }
  PU_N_branch = 0;
  if (tree->GetBranch("PU_N") != 0) {
    PU_N_branch = tree->GetBranch("PU_N");
    if (PU_N_branch) { PU_N_branch->SetAddress(&PU_N_); }
  }
  ele_kfchi2_branch = 0;
  if (tree->GetBranch("ele_kfchi2") != 0) {
    ele_kfchi2_branch = tree->GetBranch("ele_kfchi2");
    if (ele_kfchi2_branch) { ele_kfchi2_branch->SetAddress(&ele_kfchi2_); }
  }
  ele_hadronicOverEm_branch = 0;
  if (tree->GetBranch("ele_hadronicOverEm") != 0) {
    ele_hadronicOverEm_branch = tree->GetBranch("ele_hadronicOverEm");
    if (ele_hadronicOverEm_branch) { ele_hadronicOverEm_branch->SetAddress(&ele_hadronicOverEm_); }
  }
  is_signal_branch = 0;
  if (tree->GetBranch("is_signal") != 0) {
    is_signal_branch = tree->GetBranch("is_signal");
    if (is_signal_branch) { is_signal_branch->SetAddress(&is_signal_); }
  }
  ele_gsfchi2_branch = 0;
  if (tree->GetBranch("ele_gsfchi2") != 0) {
    ele_gsfchi2_branch = tree->GetBranch("ele_gsfchi2");
    if (ele_gsfchi2_branch) { ele_gsfchi2_branch->SetAddress(&ele_gsfchi2_); }
  }
  ele_oldcircularity_branch = 0;
  if (tree->GetBranch("ele_oldcircularity") != 0) {
    ele_oldcircularity_branch = tree->GetBranch("ele_oldcircularity");
    if (ele_oldcircularity_branch) { ele_oldcircularity_branch->SetAddress(&ele_oldcircularity_); }
  }
  rhs_iphi_branch = 0;
  if (tree->GetBranch("rhs_iphi") != 0) {
    rhs_iphi_branch = tree->GetBranch("rhs_iphi");
    if (rhs_iphi_branch) { rhs_iphi_branch->SetAddress(&rhs_iphi_); }
  }
  ele_sclNclus_branch = 0;
  if (tree->GetBranch("ele_sclNclus") != 0) {
    ele_sclNclus_branch = tree->GetBranch("ele_sclNclus");
    if (ele_sclNclus_branch) { ele_sclNclus_branch->SetAddress(&ele_sclNclus_); }
  }
  ele_IoEmIop_branch = 0;
  if (tree->GetBranch("ele_IoEmIop") != 0) {
    ele_IoEmIop_branch = tree->GetBranch("ele_IoEmIop");
    if (ele_IoEmIop_branch) { ele_IoEmIop_branch->SetAddress(&ele_IoEmIop_); }
  }
  ele_pt_branch = 0;
  if (tree->GetBranch("ele_pt") != 0) {
    ele_pt_branch = tree->GetBranch("ele_pt");
    if (ele_pt_branch) { ele_pt_branch->SetAddress(&ele_pt_); }
  }
  seed_iphi_branch = 0;
  if (tree->GetBranch("seed_iphi") != 0) {
    seed_iphi_branch = tree->GetBranch("seed_iphi");
    if (seed_iphi_branch) { seed_iphi_branch->SetAddress(&seed_iphi_); }
  }
  vtx_N_branch = 0;
  if (tree->GetBranch("vtx_N") != 0) {
    vtx_N_branch = tree->GetBranch("vtx_N");
    if (vtx_N_branch) { vtx_N_branch->SetAddress(&vtx_N_); }
  }
  ele_ID2_branch = 0;
  if (tree->GetBranch("ele_ID2") != 0) {
    ele_ID2_branch = tree->GetBranch("ele_ID2");
    if (ele_ID2_branch) { ele_ID2_branch->SetAddress(&ele_ID2_); }
  }
  ele_valid_hits_branch = 0;
  if (tree->GetBranch("ele_valid_hits") != 0) {
    ele_valid_hits_branch = tree->GetBranch("ele_valid_hits");
    if (ele_valid_hits_branch) { ele_valid_hits_branch->SetAddress(&ele_valid_hits_); }
  }
  ele_olde55_branch = 0;
  if (tree->GetBranch("ele_olde55") != 0) {
    ele_olde55_branch = tree->GetBranch("ele_olde55");
    if (ele_olde55_branch) { ele_olde55_branch->SetAddress(&ele_olde55_); }
  }
  ele_oldr9_branch = 0;
  if (tree->GetBranch("ele_oldr9") != 0) {
    ele_oldr9_branch = tree->GetBranch("ele_oldr9");
    if (ele_oldr9_branch) { ele_oldr9_branch->SetAddress(&ele_oldr9_); }
  }
  ele_isEERingGap_branch = 0;
  if (tree->GetBranch("ele_isEERingGap") != 0) {
    ele_isEERingGap_branch = tree->GetBranch("ele_isEERingGap");
    if (ele_isEERingGap_branch) { ele_isEERingGap_branch->SetAddress(&ele_isEERingGap_); }
  }
  ele_deltaetain_branch = 0;
  if (tree->GetBranch("ele_deltaetain") != 0) {
    ele_deltaetain_branch = tree->GetBranch("ele_deltaetain");
    if (ele_deltaetain_branch) { ele_deltaetain_branch->SetAddress(&ele_deltaetain_); }
  }
  ele_lost_hits_branch = 0;
  if (tree->GetBranch("ele_lost_hits") != 0) {
    ele_lost_hits_branch = tree->GetBranch("ele_lost_hits");
    if (ele_lost_hits_branch) { ele_lost_hits_branch->SetAddress(&ele_lost_hits_); }
  }
  mc_pho_mother_status_branch = 0;
  if (tree->GetBranch("mc_pho_mother_status") != 0) {
    mc_pho_mother_status_branch = tree->GetBranch("mc_pho_mother_status");
    if (mc_pho_mother_status_branch) { mc_pho_mother_status_branch->SetAddress(&mc_pho_mother_status_); }
  }
  ele_oldhebc_branch = 0;
  if (tree->GetBranch("ele_oldhebc") != 0) {
    ele_oldhebc_branch = tree->GetBranch("ele_oldhebc");
    if (ele_oldhebc_branch) { ele_oldhebc_branch->SetAddress(&ele_oldhebc_); }
  }
  ele_conv_radius_branch = 0;
  if (tree->GetBranch("ele_conv_radius") != 0) {
    ele_conv_radius_branch = tree->GetBranch("ele_conv_radius");
    if (ele_conv_radius_branch) { ele_conv_radius_branch->SetAddress(&ele_conv_radius_); }
  }
  ele_scletawidth_branch = 0;
  if (tree->GetBranch("ele_scletawidth") != 0) {
    ele_scletawidth_branch = tree->GetBranch("ele_scletawidth");
    if (ele_scletawidth_branch) { ele_scletawidth_branch->SetAddress(&ele_scletawidth_); }
  }
  ele_ID1_branch = 0;
  if (tree->GetBranch("ele_ID1") != 0) {
    ele_ID1_branch = tree->GetBranch("ele_ID1");
    if (ele_ID1_branch) { ele_ID1_branch->SetAddress(&ele_ID1_); }
  }
  ele_conv_dcot_branch = 0;
  if (tree->GetBranch("ele_conv_dcot") != 0) {
    ele_conv_dcot_branch = tree->GetBranch("ele_conv_dcot");
    if (ele_conv_dcot_branch) { ele_conv_dcot_branch->SetAddress(&ele_conv_dcot_); }
  }
  mc_pho_mother_id_branch = 0;
  if (tree->GetBranch("mc_pho_mother_id") != 0) {
    mc_pho_mother_id_branch = tree->GetBranch("mc_pho_mother_id");
    if (mc_pho_mother_id_branch) { mc_pho_mother_id_branch->SetAddress(&mc_pho_mother_id_); }
  }
  rhs_ieta_branch = 0;
  if (tree->GetBranch("rhs_ieta") != 0) {
    rhs_ieta_branch = tree->GetBranch("rhs_ieta");
    if (rhs_ieta_branch) { rhs_ieta_branch->SetAddress(&rhs_ieta_); }
  }
  scl_eta_branch = 0;
  if (tree->GetBranch("scl_eta") != 0) {
    scl_eta_branch = tree->GetBranch("scl_eta");
    if (scl_eta_branch) { scl_eta_branch->SetAddress(&scl_eta_); }
  }
  nRun_branch = 0;
  if (tree->GetBranch("nRun") != 0) {
    nRun_branch = tree->GetBranch("nRun");
    if (nRun_branch) { nRun_branch->SetAddress(&nRun_); }
  }
  ele_isEBEtaGap_branch = 0;
  if (tree->GetBranch("ele_isEBEtaGap") != 0) {
    ele_isEBEtaGap_branch = tree->GetBranch("ele_isEBEtaGap");
    if (ele_isEBEtaGap_branch) { ele_isEBEtaGap_branch->SetAddress(&ele_isEBEtaGap_); }
  }
  seed_ieta_branch = 0;
  if (tree->GetBranch("seed_ieta") != 0) {
    seed_ieta_branch = tree->GetBranch("seed_ieta");
    if (seed_ieta_branch) { seed_ieta_branch->SetAddress(&seed_ieta_); }
  }
  ele_gsfhits_branch = 0;
  if (tree->GetBranch("ele_gsfhits") != 0) {
    ele_gsfhits_branch = tree->GetBranch("ele_gsfhits");
    if (ele_gsfhits_branch) { ele_gsfhits_branch->SetAddress(&ele_gsfhits_); }
  }
  mc_gen_ID_branch = 0;
  if (tree->GetBranch("mc_gen_ID") != 0) {
    mc_gen_ID_branch = tree->GetBranch("mc_gen_ID");
    if (mc_gen_ID_branch) { mc_gen_ID_branch->SetAddress(&mc_gen_ID_); }
  }
  ele_3x3_branch = 0;
  if (tree->GetBranch("ele_3x3") != 0) {
    ele_3x3_branch = tree->GetBranch("ele_3x3");
    if (ele_3x3_branch) { ele_3x3_branch->SetAddress(&ele_3x3_); }
  }
  ele_psEoverEraw_branch = 0;
  if (tree->GetBranch("ele_psEoverEraw") != 0) {
    ele_psEoverEraw_branch = tree->GetBranch("ele_psEoverEraw");
    if (ele_psEoverEraw_branch) { ele_psEoverEraw_branch->SetAddress(&ele_psEoverEraw_); }
  }
  ele_ID2_pass_branch = 0;
  if (tree->GetBranch("ele_ID2_pass") != 0) {
    ele_ID2_pass_branch = tree->GetBranch("ele_ID2_pass");
    if (ele_ID2_pass_branch) { ele_ID2_pass_branch->SetAddress(&ele_ID2_pass_); }
  }
  scl_Et_branch = 0;
  if (tree->GetBranch("scl_Et") != 0) {
    scl_Et_branch = tree->GetBranch("scl_Et");
    if (scl_Et_branch) { scl_Et_branch->SetAddress(&scl_Et_); }
  }
  ele_combErr_branch = 0;
  if (tree->GetBranch("ele_combErr") != 0) {
    ele_combErr_branch = tree->GetBranch("ele_combErr");
    if (ele_combErr_branch) { ele_combErr_branch->SetAddress(&ele_combErr_); }
  }
  ele_deltaphiin_branch = 0;
  if (tree->GetBranch("ele_deltaphiin") != 0) {
    ele_deltaphiin_branch = tree->GetBranch("ele_deltaphiin");
    if (ele_deltaphiin_branch) { ele_deltaphiin_branch->SetAddress(&ele_deltaphiin_); }
  }
  seed_e_branch = 0;
  if (tree->GetBranch("seed_e") != 0) {
    seed_e_branch = tree->GetBranch("seed_e");
    if (seed_e_branch) { seed_e_branch->SetAddress(&seed_e_); }
  }
  ele_trackErr_branch = 0;
  if (tree->GetBranch("ele_trackErr") != 0) {
    ele_trackErr_branch = tree->GetBranch("ele_trackErr");
    if (ele_trackErr_branch) { ele_trackErr_branch->SetAddress(&ele_trackErr_); }
  }
  ele_PFcombErr_branch = 0;
  if (tree->GetBranch("ele_PFcombErr") != 0) {
    ele_PFcombErr_branch = tree->GetBranch("ele_PFcombErr");
    if (ele_PFcombErr_branch) { ele_PFcombErr_branch->SetAddress(&ele_PFcombErr_); }
  }
  ele_ID1_cat_branch = 0;
  if (tree->GetBranch("ele_ID1_cat") != 0) {
    ele_ID1_cat_branch = tree->GetBranch("ele_ID1_cat");
    if (ele_ID1_cat_branch) { ele_ID1_cat_branch->SetAddress(&ele_ID1_cat_); }
  }
  ele_oldsigmaiphiiphi_branch = 0;
  if (tree->GetBranch("ele_oldsigmaiphiiphi") != 0) {
    ele_oldsigmaiphiiphi_branch = tree->GetBranch("ele_oldsigmaiphiiphi");
    if (ele_oldsigmaiphiiphi_branch) { ele_oldsigmaiphiiphi_branch->SetAddress(&ele_oldsigmaiphiiphi_); }
  }
  ele_SCfbrem_branch = 0;
  if (tree->GetBranch("ele_SCfbrem") != 0) {
    ele_SCfbrem_branch = tree->GetBranch("ele_SCfbrem");
    if (ele_SCfbrem_branch) { ele_SCfbrem_branch->SetAddress(&ele_SCfbrem_); }
  }
  ele_vtxconv_branch = 0;
  if (tree->GetBranch("ele_vtxconv") != 0) {
    ele_vtxconv_branch = tree->GetBranch("ele_vtxconv");
    if (ele_vtxconv_branch) { ele_vtxconv_branch->SetAddress(&ele_vtxconv_); }
  }
  ele_dz_branch = 0;
  if (tree->GetBranch("ele_dz") != 0) {
    ele_dz_branch = tree->GetBranch("ele_dz");
    if (ele_dz_branch) { ele_dz_branch->SetAddress(&ele_dz_); }
  }
  nLumi_branch = 0;
  if (tree->GetBranch("nLumi") != 0) {
    nLumi_branch = tree->GetBranch("nLumi");
    if (nLumi_branch) { nLumi_branch->SetAddress(&nLumi_); }
  }
  ele_echarge_branch = 0;
  if (tree->GetBranch("ele_echarge") != 0) {
    ele_echarge_branch = tree->GetBranch("ele_echarge");
    if (ele_echarge_branch) { ele_echarge_branch->SetAddress(&ele_echarge_); }
  }
  ele_ID1_pass_branch = 0;
  if (tree->GetBranch("ele_ID1_pass") != 0) {
    ele_ID1_pass_branch = tree->GetBranch("ele_ID1_pass");
    if (ele_ID1_pass_branch) { ele_ID1_pass_branch->SetAddress(&ele_ID1_pass_); }
  }
  ele_isEEDeeGap_branch = 0;
  if (tree->GetBranch("ele_isEEDeeGap") != 0) {
    ele_isEEDeeGap_branch = tree->GetBranch("ele_isEEDeeGap");
    if (ele_isEEDeeGap_branch) { ele_isEEDeeGap_branch->SetAddress(&ele_isEEDeeGap_); }
  }
  ele_conversionVertexFitProbability_branch = 0;
  if (tree->GetBranch("ele_conversionVertexFitProbability") != 0) {
    ele_conversionVertexFitProbability_branch = tree->GetBranch("ele_conversionVertexFitProbability");
    if (ele_conversionVertexFitProbability_branch) { ele_conversionVertexFitProbability_branch->SetAddress(&ele_conversionVertexFitProbability_); }
  }
  ele_olde25max_branch = 0;
  if (tree->GetBranch("ele_olde25max") != 0) {
    ele_olde25max_branch = tree->GetBranch("ele_olde25max");
    if (ele_olde25max_branch) { ele_olde25max_branch->SetAddress(&ele_olde25max_); }
  }
  ele_ecalDrivenSeed_branch = 0;
  if (tree->GetBranch("ele_ecalDrivenSeed") != 0) {
    ele_ecalDrivenSeed_branch = tree->GetBranch("ele_ecalDrivenSeed");
    if (ele_ecalDrivenSeed_branch) { ele_ecalDrivenSeed_branch->SetAddress(&ele_ecalDrivenSeed_); }
  }
  ele_eta_branch = 0;
  if (tree->GetBranch("ele_eta") != 0) {
    ele_eta_branch = tree->GetBranch("ele_eta");
    if (ele_eta_branch) { ele_eta_branch->SetAddress(&ele_eta_); }
  }
  ele_oldhe_branch = 0;
  if (tree->GetBranch("ele_oldhe") != 0) {
    ele_oldhe_branch = tree->GetBranch("ele_oldhe");
    if (ele_oldhe_branch) { ele_oldhe_branch->SetAddress(&ele_oldhe_); }
  }
  ele_olde15_branch = 0;
  if (tree->GetBranch("ele_olde15") != 0) {
    ele_olde15_branch = tree->GetBranch("ele_olde15");
    if (ele_olde15_branch) { ele_olde15_branch->SetAddress(&ele_olde15_); }
  }
  mc_ele_matchMother_PDGID_branch = 0;
  if (tree->GetBranch("mc_ele_matchMother_PDGID") != 0) {
    mc_ele_matchMother_PDGID_branch = tree->GetBranch("mc_ele_matchMother_PDGID");
    if (mc_ele_matchMother_PDGID_branch) { mc_ele_matchMother_PDGID_branch->SetAddress(&mc_ele_matchMother_PDGID_); }
  }
  ele_ep_branch = 0;
  if (tree->GetBranch("ele_ep") != 0) {
    ele_ep_branch = tree->GetBranch("ele_ep");
    if (ele_ep_branch) { ele_ep_branch->SetAddress(&ele_ep_); }
  }
  ele_deltaetaseed_branch = 0;
  if (tree->GetBranch("ele_deltaetaseed") != 0) {
    ele_deltaetaseed_branch = tree->GetBranch("ele_deltaetaseed");
    if (ele_deltaetaseed_branch) { ele_deltaetaseed_branch->SetAddress(&ele_deltaetaseed_); }
  }
  ele_conv_dist_branch = 0;
  if (tree->GetBranch("ele_conv_dist") != 0) {
    ele_conv_dist_branch = tree->GetBranch("ele_conv_dist");
    if (ele_conv_dist_branch) { ele_conv_dist_branch->SetAddress(&ele_conv_dist_); }
  }
  ele_expected_inner_hits_branch = 0;
  if (tree->GetBranch("ele_expected_inner_hits") != 0) {
    ele_expected_inner_hits_branch = tree->GetBranch("ele_expected_inner_hits");
    if (ele_expected_inner_hits_branch) { ele_expected_inner_hits_branch->SetAddress(&ele_expected_inner_hits_); }
  }
  mc_ele_matchMother_status_branch = 0;
  if (tree->GetBranch("mc_ele_matchMother_status") != 0) {
    mc_ele_matchMother_status_branch = tree->GetBranch("mc_ele_matchMother_status");
    if (mc_ele_matchMother_status_branch) { mc_ele_matchMother_status_branch->SetAddress(&mc_ele_matchMother_status_); }
  }
  tree->SetMakeClass(0);
}
void CMS3::GetEntry(unsigned int idx) {
  index = idx;
  ele_full5x5_hcalOverEcal_isLoaded = false;
  MC_TrueNumInteractions_isLoaded = false;
  ele_fbrem_isLoaded = false;
  ele_oldsirir_isLoaded = false;
  ele_N_isLoaded = false;
  mc_gen_pt_isLoaded = false;
  mc_event_weight_isLoaded = false;
  ele_kfhits_isLoaded = false;
  ele_trackerDrivenSeed_isLoaded = false;
  mc_ele_matchedFromCB2_isLoaded = false;
  rhs_e_isLoaded = false;
  ele_eelepout_isLoaded = false;
  ele_isEBEEGap_isLoaded = false;
  ele_oldsigmaietaieta_isLoaded = false;
  nEvent_isLoaded = false;
  ele_nbrem_isLoaded = false;
  ele_index_isLoaded = false;
  scl_phi_isLoaded = false;
  ele_isEBPhiGap_isLoaded = false;
  ele_isEE_isLoaded = false;
  ele_isEB_isLoaded = false;
  ele_SIP_isLoaded = false;
  ele_sclphiwidth_isLoaded = false;
  mc_ele_matchedFromCB_isLoaded = false;
  ele_ID2_cat_isLoaded = false;
  scl_E_isLoaded = false;
  ele_eClass_isLoaded = false;
  mc_gen_eta_isLoaded = false;
  ele_sclRawE_isLoaded = false;
  PU_N_isLoaded = false;
  ele_kfchi2_isLoaded = false;
  ele_hadronicOverEm_isLoaded = false;
  is_signal_isLoaded = false;
  ele_gsfchi2_isLoaded = false;
  ele_oldcircularity_isLoaded = false;
  rhs_iphi_isLoaded = false;
  ele_sclNclus_isLoaded = false;
  ele_IoEmIop_isLoaded = false;
  ele_pt_isLoaded = false;
  seed_iphi_isLoaded = false;
  vtx_N_isLoaded = false;
  ele_ID2_isLoaded = false;
  ele_valid_hits_isLoaded = false;
  ele_olde55_isLoaded = false;
  ele_oldr9_isLoaded = false;
  ele_isEERingGap_isLoaded = false;
  ele_deltaetain_isLoaded = false;
  ele_lost_hits_isLoaded = false;
  mc_pho_mother_status_isLoaded = false;
  ele_oldhebc_isLoaded = false;
  ele_conv_radius_isLoaded = false;
  ele_scletawidth_isLoaded = false;
  ele_ID1_isLoaded = false;
  ele_conv_dcot_isLoaded = false;
  mc_pho_mother_id_isLoaded = false;
  rhs_ieta_isLoaded = false;
  scl_eta_isLoaded = false;
  nRun_isLoaded = false;
  ele_isEBEtaGap_isLoaded = false;
  seed_ieta_isLoaded = false;
  ele_gsfhits_isLoaded = false;
  mc_gen_ID_isLoaded = false;
  ele_3x3_isLoaded = false;
  ele_psEoverEraw_isLoaded = false;
  ele_ID2_pass_isLoaded = false;
  scl_Et_isLoaded = false;
  ele_combErr_isLoaded = false;
  ele_deltaphiin_isLoaded = false;
  seed_e_isLoaded = false;
  ele_trackErr_isLoaded = false;
  ele_PFcombErr_isLoaded = false;
  ele_ID1_cat_isLoaded = false;
  ele_oldsigmaiphiiphi_isLoaded = false;
  ele_SCfbrem_isLoaded = false;
  ele_vtxconv_isLoaded = false;
  ele_dz_isLoaded = false;
  nLumi_isLoaded = false;
  ele_echarge_isLoaded = false;
  ele_ID1_pass_isLoaded = false;
  ele_isEEDeeGap_isLoaded = false;
  ele_conversionVertexFitProbability_isLoaded = false;
  ele_olde25max_isLoaded = false;
  ele_ecalDrivenSeed_isLoaded = false;
  ele_eta_isLoaded = false;
  ele_oldhe_isLoaded = false;
  ele_olde15_isLoaded = false;
  mc_ele_matchMother_PDGID_isLoaded = false;
  ele_ep_isLoaded = false;
  ele_deltaetaseed_isLoaded = false;
  ele_conv_dist_isLoaded = false;
  ele_expected_inner_hits_isLoaded = false;
  mc_ele_matchMother_status_isLoaded = false;
}
void CMS3::LoadAllBranches() {
  if (ele_full5x5_hcalOverEcal_branch != 0) ele_full5x5_hcalOverEcal();
  if (MC_TrueNumInteractions_branch != 0) MC_TrueNumInteractions();
  if (ele_fbrem_branch != 0) ele_fbrem();
  if (ele_oldsirir_branch != 0) ele_oldsirir();
  if (ele_N_branch != 0) ele_N();
  if (mc_gen_pt_branch != 0) mc_gen_pt();
  if (mc_event_weight_branch != 0) mc_event_weight();
  if (ele_kfhits_branch != 0) ele_kfhits();
  if (ele_trackerDrivenSeed_branch != 0) ele_trackerDrivenSeed();
  if (mc_ele_matchedFromCB2_branch != 0) mc_ele_matchedFromCB2();
  if (rhs_e_branch != 0) rhs_e();
  if (ele_eelepout_branch != 0) ele_eelepout();
  if (ele_isEBEEGap_branch != 0) ele_isEBEEGap();
  if (ele_oldsigmaietaieta_branch != 0) ele_oldsigmaietaieta();
  if (nEvent_branch != 0) nEvent();
  if (ele_nbrem_branch != 0) ele_nbrem();
  if (ele_index_branch != 0) ele_index();
  if (scl_phi_branch != 0) scl_phi();
  if (ele_isEBPhiGap_branch != 0) ele_isEBPhiGap();
  if (ele_isEE_branch != 0) ele_isEE();
  if (ele_isEB_branch != 0) ele_isEB();
  if (ele_SIP_branch != 0) ele_SIP();
  if (ele_sclphiwidth_branch != 0) ele_sclphiwidth();
  if (mc_ele_matchedFromCB_branch != 0) mc_ele_matchedFromCB();
  if (ele_ID2_cat_branch != 0) ele_ID2_cat();
  if (scl_E_branch != 0) scl_E();
  if (ele_eClass_branch != 0) ele_eClass();
  if (mc_gen_eta_branch != 0) mc_gen_eta();
  if (ele_sclRawE_branch != 0) ele_sclRawE();
  if (PU_N_branch != 0) PU_N();
  if (ele_kfchi2_branch != 0) ele_kfchi2();
  if (ele_hadronicOverEm_branch != 0) ele_hadronicOverEm();
  if (is_signal_branch != 0) is_signal();
  if (ele_gsfchi2_branch != 0) ele_gsfchi2();
  if (ele_oldcircularity_branch != 0) ele_oldcircularity();
  if (rhs_iphi_branch != 0) rhs_iphi();
  if (ele_sclNclus_branch != 0) ele_sclNclus();
  if (ele_IoEmIop_branch != 0) ele_IoEmIop();
  if (ele_pt_branch != 0) ele_pt();
  if (seed_iphi_branch != 0) seed_iphi();
  if (vtx_N_branch != 0) vtx_N();
  if (ele_ID2_branch != 0) ele_ID2();
  if (ele_valid_hits_branch != 0) ele_valid_hits();
  if (ele_olde55_branch != 0) ele_olde55();
  if (ele_oldr9_branch != 0) ele_oldr9();
  if (ele_isEERingGap_branch != 0) ele_isEERingGap();
  if (ele_deltaetain_branch != 0) ele_deltaetain();
  if (ele_lost_hits_branch != 0) ele_lost_hits();
  if (mc_pho_mother_status_branch != 0) mc_pho_mother_status();
  if (ele_oldhebc_branch != 0) ele_oldhebc();
  if (ele_conv_radius_branch != 0) ele_conv_radius();
  if (ele_scletawidth_branch != 0) ele_scletawidth();
  if (ele_ID1_branch != 0) ele_ID1();
  if (ele_conv_dcot_branch != 0) ele_conv_dcot();
  if (mc_pho_mother_id_branch != 0) mc_pho_mother_id();
  if (rhs_ieta_branch != 0) rhs_ieta();
  if (scl_eta_branch != 0) scl_eta();
  if (nRun_branch != 0) nRun();
  if (ele_isEBEtaGap_branch != 0) ele_isEBEtaGap();
  if (seed_ieta_branch != 0) seed_ieta();
  if (ele_gsfhits_branch != 0) ele_gsfhits();
  if (mc_gen_ID_branch != 0) mc_gen_ID();
  if (ele_3x3_branch != 0) ele_3x3();
  if (ele_psEoverEraw_branch != 0) ele_psEoverEraw();
  if (ele_ID2_pass_branch != 0) ele_ID2_pass();
  if (scl_Et_branch != 0) scl_Et();
  if (ele_combErr_branch != 0) ele_combErr();
  if (ele_deltaphiin_branch != 0) ele_deltaphiin();
  if (seed_e_branch != 0) seed_e();
  if (ele_trackErr_branch != 0) ele_trackErr();
  if (ele_PFcombErr_branch != 0) ele_PFcombErr();
  if (ele_ID1_cat_branch != 0) ele_ID1_cat();
  if (ele_oldsigmaiphiiphi_branch != 0) ele_oldsigmaiphiiphi();
  if (ele_SCfbrem_branch != 0) ele_SCfbrem();
  if (ele_vtxconv_branch != 0) ele_vtxconv();
  if (ele_dz_branch != 0) ele_dz();
  if (nLumi_branch != 0) nLumi();
  if (ele_echarge_branch != 0) ele_echarge();
  if (ele_ID1_pass_branch != 0) ele_ID1_pass();
  if (ele_isEEDeeGap_branch != 0) ele_isEEDeeGap();
  if (ele_conversionVertexFitProbability_branch != 0) ele_conversionVertexFitProbability();
  if (ele_olde25max_branch != 0) ele_olde25max();
  if (ele_ecalDrivenSeed_branch != 0) ele_ecalDrivenSeed();
  if (ele_eta_branch != 0) ele_eta();
  if (ele_oldhe_branch != 0) ele_oldhe();
  if (ele_olde15_branch != 0) ele_olde15();
  if (mc_ele_matchMother_PDGID_branch != 0) mc_ele_matchMother_PDGID();
  if (ele_ep_branch != 0) ele_ep();
  if (ele_deltaetaseed_branch != 0) ele_deltaetaseed();
  if (ele_conv_dist_branch != 0) ele_conv_dist();
  if (ele_expected_inner_hits_branch != 0) ele_expected_inner_hits();
  if (mc_ele_matchMother_status_branch != 0) mc_ele_matchMother_status();
}
const float &CMS3::ele_full5x5_hcalOverEcal() {
  if (not ele_full5x5_hcalOverEcal_isLoaded) {
    if (ele_full5x5_hcalOverEcal_branch != 0) {
      ele_full5x5_hcalOverEcal_branch->GetEntry(index);
    } else {
      printf("branch ele_full5x5_hcalOverEcal_branch does not exist!\n");
      exit(1);
    }
    ele_full5x5_hcalOverEcal_isLoaded = true;
  }
  return ele_full5x5_hcalOverEcal_;
}
const int &CMS3::MC_TrueNumInteractions() {
  if (not MC_TrueNumInteractions_isLoaded) {
    if (MC_TrueNumInteractions_branch != 0) {
      MC_TrueNumInteractions_branch->GetEntry(index);
    } else {
      printf("branch MC_TrueNumInteractions_branch does not exist!\n");
      exit(1);
    }
    MC_TrueNumInteractions_isLoaded = true;
  }
  return MC_TrueNumInteractions_;
}
const float &CMS3::ele_fbrem() {
  if (not ele_fbrem_isLoaded) {
    if (ele_fbrem_branch != 0) {
      ele_fbrem_branch->GetEntry(index);
    } else {
      printf("branch ele_fbrem_branch does not exist!\n");
      exit(1);
    }
    ele_fbrem_isLoaded = true;
  }
  return ele_fbrem_;
}
const float &CMS3::ele_oldsirir() {
  if (not ele_oldsirir_isLoaded) {
    if (ele_oldsirir_branch != 0) {
      ele_oldsirir_branch->GetEntry(index);
    } else {
      printf("branch ele_oldsirir_branch does not exist!\n");
      exit(1);
    }
    ele_oldsirir_isLoaded = true;
  }
  return ele_oldsirir_;
}
const int &CMS3::ele_N() {
  if (not ele_N_isLoaded) {
    if (ele_N_branch != 0) {
      ele_N_branch->GetEntry(index);
    } else {
      printf("branch ele_N_branch does not exist!\n");
      exit(1);
    }
    ele_N_isLoaded = true;
  }
  return ele_N_;
}
const float &CMS3::mc_gen_pt() {
  if (not mc_gen_pt_isLoaded) {
    if (mc_gen_pt_branch != 0) {
      mc_gen_pt_branch->GetEntry(index);
    } else {
      printf("branch mc_gen_pt_branch does not exist!\n");
      exit(1);
    }
    mc_gen_pt_isLoaded = true;
  }
  return mc_gen_pt_;
}
const float &CMS3::mc_event_weight() {
  if (not mc_event_weight_isLoaded) {
    if (mc_event_weight_branch != 0) {
      mc_event_weight_branch->GetEntry(index);
    } else {
      printf("branch mc_event_weight_branch does not exist!\n");
      exit(1);
    }
    mc_event_weight_isLoaded = true;
  }
  return mc_event_weight_;
}
const int &CMS3::ele_kfhits() {
  if (not ele_kfhits_isLoaded) {
    if (ele_kfhits_branch != 0) {
      ele_kfhits_branch->GetEntry(index);
    } else {
      printf("branch ele_kfhits_branch does not exist!\n");
      exit(1);
    }
    ele_kfhits_isLoaded = true;
  }
  return ele_kfhits_;
}
const bool &CMS3::ele_trackerDrivenSeed() {
  if (not ele_trackerDrivenSeed_isLoaded) {
    if (ele_trackerDrivenSeed_branch != 0) {
      ele_trackerDrivenSeed_branch->GetEntry(index);
    } else {
      printf("branch ele_trackerDrivenSeed_branch does not exist!\n");
      exit(1);
    }
    ele_trackerDrivenSeed_isLoaded = true;
  }
  return ele_trackerDrivenSeed_;
}
const int &CMS3::mc_ele_matchedFromCB2() {
  if (not mc_ele_matchedFromCB2_isLoaded) {
    if (mc_ele_matchedFromCB2_branch != 0) {
      mc_ele_matchedFromCB2_branch->GetEntry(index);
    } else {
      printf("branch mc_ele_matchedFromCB2_branch does not exist!\n");
      exit(1);
    }
    mc_ele_matchedFromCB2_isLoaded = true;
  }
  return mc_ele_matchedFromCB2_;
}
const vector<float> &CMS3::rhs_e() {
  if (not rhs_e_isLoaded) {
    if (rhs_e_branch != 0) {
      rhs_e_branch->GetEntry(index);
    } else {
      printf("branch rhs_e_branch does not exist!\n");
      exit(1);
    }
    rhs_e_isLoaded = true;
  }
  return rhs_e_;
}
const float &CMS3::ele_eelepout() {
  if (not ele_eelepout_isLoaded) {
    if (ele_eelepout_branch != 0) {
      ele_eelepout_branch->GetEntry(index);
    } else {
      printf("branch ele_eelepout_branch does not exist!\n");
      exit(1);
    }
    ele_eelepout_isLoaded = true;
  }
  return ele_eelepout_;
}
const bool &CMS3::ele_isEBEEGap() {
  if (not ele_isEBEEGap_isLoaded) {
    if (ele_isEBEEGap_branch != 0) {
      ele_isEBEEGap_branch->GetEntry(index);
    } else {
      printf("branch ele_isEBEEGap_branch does not exist!\n");
      exit(1);
    }
    ele_isEBEEGap_isLoaded = true;
  }
  return ele_isEBEEGap_;
}
const float &CMS3::ele_oldsigmaietaieta() {
  if (not ele_oldsigmaietaieta_isLoaded) {
    if (ele_oldsigmaietaieta_branch != 0) {
      ele_oldsigmaietaieta_branch->GetEntry(index);
    } else {
      printf("branch ele_oldsigmaietaieta_branch does not exist!\n");
      exit(1);
    }
    ele_oldsigmaietaieta_isLoaded = true;
  }
  return ele_oldsigmaietaieta_;
}
const int &CMS3::nEvent() {
  if (not nEvent_isLoaded) {
    if (nEvent_branch != 0) {
      nEvent_branch->GetEntry(index);
    } else {
      printf("branch nEvent_branch does not exist!\n");
      exit(1);
    }
    nEvent_isLoaded = true;
  }
  return nEvent_;
}
const int &CMS3::ele_nbrem() {
  if (not ele_nbrem_isLoaded) {
    if (ele_nbrem_branch != 0) {
      ele_nbrem_branch->GetEntry(index);
    } else {
      printf("branch ele_nbrem_branch does not exist!\n");
      exit(1);
    }
    ele_nbrem_isLoaded = true;
  }
  return ele_nbrem_;
}
const int &CMS3::ele_index() {
  if (not ele_index_isLoaded) {
    if (ele_index_branch != 0) {
      ele_index_branch->GetEntry(index);
    } else {
      printf("branch ele_index_branch does not exist!\n");
      exit(1);
    }
    ele_index_isLoaded = true;
  }
  return ele_index_;
}
const float &CMS3::scl_phi() {
  if (not scl_phi_isLoaded) {
    if (scl_phi_branch != 0) {
      scl_phi_branch->GetEntry(index);
    } else {
      printf("branch scl_phi_branch does not exist!\n");
      exit(1);
    }
    scl_phi_isLoaded = true;
  }
  return scl_phi_;
}
const bool &CMS3::ele_isEBPhiGap() {
  if (not ele_isEBPhiGap_isLoaded) {
    if (ele_isEBPhiGap_branch != 0) {
      ele_isEBPhiGap_branch->GetEntry(index);
    } else {
      printf("branch ele_isEBPhiGap_branch does not exist!\n");
      exit(1);
    }
    ele_isEBPhiGap_isLoaded = true;
  }
  return ele_isEBPhiGap_;
}
const bool &CMS3::ele_isEE() {
  if (not ele_isEE_isLoaded) {
    if (ele_isEE_branch != 0) {
      ele_isEE_branch->GetEntry(index);
    } else {
      printf("branch ele_isEE_branch does not exist!\n");
      exit(1);
    }
    ele_isEE_isLoaded = true;
  }
  return ele_isEE_;
}
const bool &CMS3::ele_isEB() {
  if (not ele_isEB_isLoaded) {
    if (ele_isEB_branch != 0) {
      ele_isEB_branch->GetEntry(index);
    } else {
      printf("branch ele_isEB_branch does not exist!\n");
      exit(1);
    }
    ele_isEB_isLoaded = true;
  }
  return ele_isEB_;
}
const float &CMS3::ele_SIP() {
  if (not ele_SIP_isLoaded) {
    if (ele_SIP_branch != 0) {
      ele_SIP_branch->GetEntry(index);
    } else {
      printf("branch ele_SIP_branch does not exist!\n");
      exit(1);
    }
    ele_SIP_isLoaded = true;
  }
  return ele_SIP_;
}
const float &CMS3::ele_sclphiwidth() {
  if (not ele_sclphiwidth_isLoaded) {
    if (ele_sclphiwidth_branch != 0) {
      ele_sclphiwidth_branch->GetEntry(index);
    } else {
      printf("branch ele_sclphiwidth_branch does not exist!\n");
      exit(1);
    }
    ele_sclphiwidth_isLoaded = true;
  }
  return ele_sclphiwidth_;
}
const int &CMS3::mc_ele_matchedFromCB() {
  if (not mc_ele_matchedFromCB_isLoaded) {
    if (mc_ele_matchedFromCB_branch != 0) {
      mc_ele_matchedFromCB_branch->GetEntry(index);
    } else {
      printf("branch mc_ele_matchedFromCB_branch does not exist!\n");
      exit(1);
    }
    mc_ele_matchedFromCB_isLoaded = true;
  }
  return mc_ele_matchedFromCB_;
}
const int &CMS3::ele_ID2_cat() {
  if (not ele_ID2_cat_isLoaded) {
    if (ele_ID2_cat_branch != 0) {
      ele_ID2_cat_branch->GetEntry(index);
    } else {
      printf("branch ele_ID2_cat_branch does not exist!\n");
      exit(1);
    }
    ele_ID2_cat_isLoaded = true;
  }
  return ele_ID2_cat_;
}
const float &CMS3::scl_E() {
  if (not scl_E_isLoaded) {
    if (scl_E_branch != 0) {
      scl_E_branch->GetEntry(index);
    } else {
      printf("branch scl_E_branch does not exist!\n");
      exit(1);
    }
    scl_E_isLoaded = true;
  }
  return scl_E_;
}
const int &CMS3::ele_eClass() {
  if (not ele_eClass_isLoaded) {
    if (ele_eClass_branch != 0) {
      ele_eClass_branch->GetEntry(index);
    } else {
      printf("branch ele_eClass_branch does not exist!\n");
      exit(1);
    }
    ele_eClass_isLoaded = true;
  }
  return ele_eClass_;
}
const float &CMS3::mc_gen_eta() {
  if (not mc_gen_eta_isLoaded) {
    if (mc_gen_eta_branch != 0) {
      mc_gen_eta_branch->GetEntry(index);
    } else {
      printf("branch mc_gen_eta_branch does not exist!\n");
      exit(1);
    }
    mc_gen_eta_isLoaded = true;
  }
  return mc_gen_eta_;
}
const float &CMS3::ele_sclRawE() {
  if (not ele_sclRawE_isLoaded) {
    if (ele_sclRawE_branch != 0) {
      ele_sclRawE_branch->GetEntry(index);
    } else {
      printf("branch ele_sclRawE_branch does not exist!\n");
      exit(1);
    }
    ele_sclRawE_isLoaded = true;
  }
  return ele_sclRawE_;
}
const int &CMS3::PU_N() {
  if (not PU_N_isLoaded) {
    if (PU_N_branch != 0) {
      PU_N_branch->GetEntry(index);
    } else {
      printf("branch PU_N_branch does not exist!\n");
      exit(1);
    }
    PU_N_isLoaded = true;
  }
  return PU_N_;
}
const float &CMS3::ele_kfchi2() {
  if (not ele_kfchi2_isLoaded) {
    if (ele_kfchi2_branch != 0) {
      ele_kfchi2_branch->GetEntry(index);
    } else {
      printf("branch ele_kfchi2_branch does not exist!\n");
      exit(1);
    }
    ele_kfchi2_isLoaded = true;
  }
  return ele_kfchi2_;
}
const float &CMS3::ele_hadronicOverEm() {
  if (not ele_hadronicOverEm_isLoaded) {
    if (ele_hadronicOverEm_branch != 0) {
      ele_hadronicOverEm_branch->GetEntry(index);
    } else {
      printf("branch ele_hadronicOverEm_branch does not exist!\n");
      exit(1);
    }
    ele_hadronicOverEm_isLoaded = true;
  }
  return ele_hadronicOverEm_;
}
const bool &CMS3::is_signal() {
  if (not is_signal_isLoaded) {
    if (is_signal_branch != 0) {
      is_signal_branch->GetEntry(index);
    } else {
      printf("branch is_signal_branch does not exist!\n");
      exit(1);
    }
    is_signal_isLoaded = true;
  }
  return is_signal_;
}
const float &CMS3::ele_gsfchi2() {
  if (not ele_gsfchi2_isLoaded) {
    if (ele_gsfchi2_branch != 0) {
      ele_gsfchi2_branch->GetEntry(index);
    } else {
      printf("branch ele_gsfchi2_branch does not exist!\n");
      exit(1);
    }
    ele_gsfchi2_isLoaded = true;
  }
  return ele_gsfchi2_;
}
const float &CMS3::ele_oldcircularity() {
  if (not ele_oldcircularity_isLoaded) {
    if (ele_oldcircularity_branch != 0) {
      ele_oldcircularity_branch->GetEntry(index);
    } else {
      printf("branch ele_oldcircularity_branch does not exist!\n");
      exit(1);
    }
    ele_oldcircularity_isLoaded = true;
  }
  return ele_oldcircularity_;
}
const vector<int> &CMS3::rhs_iphi() {
  if (not rhs_iphi_isLoaded) {
    if (rhs_iphi_branch != 0) {
      rhs_iphi_branch->GetEntry(index);
    } else {
      printf("branch rhs_iphi_branch does not exist!\n");
      exit(1);
    }
    rhs_iphi_isLoaded = true;
  }
  return rhs_iphi_;
}
const int &CMS3::ele_sclNclus() {
  if (not ele_sclNclus_isLoaded) {
    if (ele_sclNclus_branch != 0) {
      ele_sclNclus_branch->GetEntry(index);
    } else {
      printf("branch ele_sclNclus_branch does not exist!\n");
      exit(1);
    }
    ele_sclNclus_isLoaded = true;
  }
  return ele_sclNclus_;
}
const float &CMS3::ele_IoEmIop() {
  if (not ele_IoEmIop_isLoaded) {
    if (ele_IoEmIop_branch != 0) {
      ele_IoEmIop_branch->GetEntry(index);
    } else {
      printf("branch ele_IoEmIop_branch does not exist!\n");
      exit(1);
    }
    ele_IoEmIop_isLoaded = true;
  }
  return ele_IoEmIop_;
}
const float &CMS3::ele_pt() {
  if (not ele_pt_isLoaded) {
    if (ele_pt_branch != 0) {
      ele_pt_branch->GetEntry(index);
    } else {
      printf("branch ele_pt_branch does not exist!\n");
      exit(1);
    }
    ele_pt_isLoaded = true;
  }
  return ele_pt_;
}
const int &CMS3::seed_iphi() {
  if (not seed_iphi_isLoaded) {
    if (seed_iphi_branch != 0) {
      seed_iphi_branch->GetEntry(index);
    } else {
      printf("branch seed_iphi_branch does not exist!\n");
      exit(1);
    }
    seed_iphi_isLoaded = true;
  }
  return seed_iphi_;
}
const int &CMS3::vtx_N() {
  if (not vtx_N_isLoaded) {
    if (vtx_N_branch != 0) {
      vtx_N_branch->GetEntry(index);
    } else {
      printf("branch vtx_N_branch does not exist!\n");
      exit(1);
    }
    vtx_N_isLoaded = true;
  }
  return vtx_N_;
}
const float &CMS3::ele_ID2() {
  if (not ele_ID2_isLoaded) {
    if (ele_ID2_branch != 0) {
      ele_ID2_branch->GetEntry(index);
    } else {
      printf("branch ele_ID2_branch does not exist!\n");
      exit(1);
    }
    ele_ID2_isLoaded = true;
  }
  return ele_ID2_;
}
const int &CMS3::ele_valid_hits() {
  if (not ele_valid_hits_isLoaded) {
    if (ele_valid_hits_branch != 0) {
      ele_valid_hits_branch->GetEntry(index);
    } else {
      printf("branch ele_valid_hits_branch does not exist!\n");
      exit(1);
    }
    ele_valid_hits_isLoaded = true;
  }
  return ele_valid_hits_;
}
const float &CMS3::ele_olde55() {
  if (not ele_olde55_isLoaded) {
    if (ele_olde55_branch != 0) {
      ele_olde55_branch->GetEntry(index);
    } else {
      printf("branch ele_olde55_branch does not exist!\n");
      exit(1);
    }
    ele_olde55_isLoaded = true;
  }
  return ele_olde55_;
}
const float &CMS3::ele_oldr9() {
  if (not ele_oldr9_isLoaded) {
    if (ele_oldr9_branch != 0) {
      ele_oldr9_branch->GetEntry(index);
    } else {
      printf("branch ele_oldr9_branch does not exist!\n");
      exit(1);
    }
    ele_oldr9_isLoaded = true;
  }
  return ele_oldr9_;
}
const bool &CMS3::ele_isEERingGap() {
  if (not ele_isEERingGap_isLoaded) {
    if (ele_isEERingGap_branch != 0) {
      ele_isEERingGap_branch->GetEntry(index);
    } else {
      printf("branch ele_isEERingGap_branch does not exist!\n");
      exit(1);
    }
    ele_isEERingGap_isLoaded = true;
  }
  return ele_isEERingGap_;
}
const float &CMS3::ele_deltaetain() {
  if (not ele_deltaetain_isLoaded) {
    if (ele_deltaetain_branch != 0) {
      ele_deltaetain_branch->GetEntry(index);
    } else {
      printf("branch ele_deltaetain_branch does not exist!\n");
      exit(1);
    }
    ele_deltaetain_isLoaded = true;
  }
  return ele_deltaetain_;
}
const int &CMS3::ele_lost_hits() {
  if (not ele_lost_hits_isLoaded) {
    if (ele_lost_hits_branch != 0) {
      ele_lost_hits_branch->GetEntry(index);
    } else {
      printf("branch ele_lost_hits_branch does not exist!\n");
      exit(1);
    }
    ele_lost_hits_isLoaded = true;
  }
  return ele_lost_hits_;
}
const int &CMS3::mc_pho_mother_status() {
  if (not mc_pho_mother_status_isLoaded) {
    if (mc_pho_mother_status_branch != 0) {
      mc_pho_mother_status_branch->GetEntry(index);
    } else {
      printf("branch mc_pho_mother_status_branch does not exist!\n");
      exit(1);
    }
    mc_pho_mother_status_isLoaded = true;
  }
  return mc_pho_mother_status_;
}
const float &CMS3::ele_oldhebc() {
  if (not ele_oldhebc_isLoaded) {
    if (ele_oldhebc_branch != 0) {
      ele_oldhebc_branch->GetEntry(index);
    } else {
      printf("branch ele_oldhebc_branch does not exist!\n");
      exit(1);
    }
    ele_oldhebc_isLoaded = true;
  }
  return ele_oldhebc_;
}
const float &CMS3::ele_conv_radius() {
  if (not ele_conv_radius_isLoaded) {
    if (ele_conv_radius_branch != 0) {
      ele_conv_radius_branch->GetEntry(index);
    } else {
      printf("branch ele_conv_radius_branch does not exist!\n");
      exit(1);
    }
    ele_conv_radius_isLoaded = true;
  }
  return ele_conv_radius_;
}
const float &CMS3::ele_scletawidth() {
  if (not ele_scletawidth_isLoaded) {
    if (ele_scletawidth_branch != 0) {
      ele_scletawidth_branch->GetEntry(index);
    } else {
      printf("branch ele_scletawidth_branch does not exist!\n");
      exit(1);
    }
    ele_scletawidth_isLoaded = true;
  }
  return ele_scletawidth_;
}
const float &CMS3::ele_ID1() {
  if (not ele_ID1_isLoaded) {
    if (ele_ID1_branch != 0) {
      ele_ID1_branch->GetEntry(index);
    } else {
      printf("branch ele_ID1_branch does not exist!\n");
      exit(1);
    }
    ele_ID1_isLoaded = true;
  }
  return ele_ID1_;
}
const float &CMS3::ele_conv_dcot() {
  if (not ele_conv_dcot_isLoaded) {
    if (ele_conv_dcot_branch != 0) {
      ele_conv_dcot_branch->GetEntry(index);
    } else {
      printf("branch ele_conv_dcot_branch does not exist!\n");
      exit(1);
    }
    ele_conv_dcot_isLoaded = true;
  }
  return ele_conv_dcot_;
}
const int &CMS3::mc_pho_mother_id() {
  if (not mc_pho_mother_id_isLoaded) {
    if (mc_pho_mother_id_branch != 0) {
      mc_pho_mother_id_branch->GetEntry(index);
    } else {
      printf("branch mc_pho_mother_id_branch does not exist!\n");
      exit(1);
    }
    mc_pho_mother_id_isLoaded = true;
  }
  return mc_pho_mother_id_;
}
const vector<int> &CMS3::rhs_ieta() {
  if (not rhs_ieta_isLoaded) {
    if (rhs_ieta_branch != 0) {
      rhs_ieta_branch->GetEntry(index);
    } else {
      printf("branch rhs_ieta_branch does not exist!\n");
      exit(1);
    }
    rhs_ieta_isLoaded = true;
  }
  return rhs_ieta_;
}
const float &CMS3::scl_eta() {
  if (not scl_eta_isLoaded) {
    if (scl_eta_branch != 0) {
      scl_eta_branch->GetEntry(index);
    } else {
      printf("branch scl_eta_branch does not exist!\n");
      exit(1);
    }
    scl_eta_isLoaded = true;
  }
  return scl_eta_;
}
const int &CMS3::nRun() {
  if (not nRun_isLoaded) {
    if (nRun_branch != 0) {
      nRun_branch->GetEntry(index);
    } else {
      printf("branch nRun_branch does not exist!\n");
      exit(1);
    }
    nRun_isLoaded = true;
  }
  return nRun_;
}
const bool &CMS3::ele_isEBEtaGap() {
  if (not ele_isEBEtaGap_isLoaded) {
    if (ele_isEBEtaGap_branch != 0) {
      ele_isEBEtaGap_branch->GetEntry(index);
    } else {
      printf("branch ele_isEBEtaGap_branch does not exist!\n");
      exit(1);
    }
    ele_isEBEtaGap_isLoaded = true;
  }
  return ele_isEBEtaGap_;
}
const int &CMS3::seed_ieta() {
  if (not seed_ieta_isLoaded) {
    if (seed_ieta_branch != 0) {
      seed_ieta_branch->GetEntry(index);
    } else {
      printf("branch seed_ieta_branch does not exist!\n");
      exit(1);
    }
    seed_ieta_isLoaded = true;
  }
  return seed_ieta_;
}
const int &CMS3::ele_gsfhits() {
  if (not ele_gsfhits_isLoaded) {
    if (ele_gsfhits_branch != 0) {
      ele_gsfhits_branch->GetEntry(index);
    } else {
      printf("branch ele_gsfhits_branch does not exist!\n");
      exit(1);
    }
    ele_gsfhits_isLoaded = true;
  }
  return ele_gsfhits_;
}
const float &CMS3::mc_gen_ID() {
  if (not mc_gen_ID_isLoaded) {
    if (mc_gen_ID_branch != 0) {
      mc_gen_ID_branch->GetEntry(index);
    } else {
      printf("branch mc_gen_ID_branch does not exist!\n");
      exit(1);
    }
    mc_gen_ID_isLoaded = true;
  }
  return mc_gen_ID_;
}
const float &CMS3::ele_3x3() {
  if (not ele_3x3_isLoaded) {
    if (ele_3x3_branch != 0) {
      ele_3x3_branch->GetEntry(index);
    } else {
      printf("branch ele_3x3_branch does not exist!\n");
      exit(1);
    }
    ele_3x3_isLoaded = true;
  }
  return ele_3x3_;
}
const float &CMS3::ele_psEoverEraw() {
  if (not ele_psEoverEraw_isLoaded) {
    if (ele_psEoverEraw_branch != 0) {
      ele_psEoverEraw_branch->GetEntry(index);
    } else {
      printf("branch ele_psEoverEraw_branch does not exist!\n");
      exit(1);
    }
    ele_psEoverEraw_isLoaded = true;
  }
  return ele_psEoverEraw_;
}
const int &CMS3::ele_ID2_pass() {
  if (not ele_ID2_pass_isLoaded) {
    if (ele_ID2_pass_branch != 0) {
      ele_ID2_pass_branch->GetEntry(index);
    } else {
      printf("branch ele_ID2_pass_branch does not exist!\n");
      exit(1);
    }
    ele_ID2_pass_isLoaded = true;
  }
  return ele_ID2_pass_;
}
const float &CMS3::scl_Et() {
  if (not scl_Et_isLoaded) {
    if (scl_Et_branch != 0) {
      scl_Et_branch->GetEntry(index);
    } else {
      printf("branch scl_Et_branch does not exist!\n");
      exit(1);
    }
    scl_Et_isLoaded = true;
  }
  return scl_Et_;
}
const float &CMS3::ele_combErr() {
  if (not ele_combErr_isLoaded) {
    if (ele_combErr_branch != 0) {
      ele_combErr_branch->GetEntry(index);
    } else {
      printf("branch ele_combErr_branch does not exist!\n");
      exit(1);
    }
    ele_combErr_isLoaded = true;
  }
  return ele_combErr_;
}
const float &CMS3::ele_deltaphiin() {
  if (not ele_deltaphiin_isLoaded) {
    if (ele_deltaphiin_branch != 0) {
      ele_deltaphiin_branch->GetEntry(index);
    } else {
      printf("branch ele_deltaphiin_branch does not exist!\n");
      exit(1);
    }
    ele_deltaphiin_isLoaded = true;
  }
  return ele_deltaphiin_;
}
const float &CMS3::seed_e() {
  if (not seed_e_isLoaded) {
    if (seed_e_branch != 0) {
      seed_e_branch->GetEntry(index);
    } else {
      printf("branch seed_e_branch does not exist!\n");
      exit(1);
    }
    seed_e_isLoaded = true;
  }
  return seed_e_;
}
const float &CMS3::ele_trackErr() {
  if (not ele_trackErr_isLoaded) {
    if (ele_trackErr_branch != 0) {
      ele_trackErr_branch->GetEntry(index);
    } else {
      printf("branch ele_trackErr_branch does not exist!\n");
      exit(1);
    }
    ele_trackErr_isLoaded = true;
  }
  return ele_trackErr_;
}
const float &CMS3::ele_PFcombErr() {
  if (not ele_PFcombErr_isLoaded) {
    if (ele_PFcombErr_branch != 0) {
      ele_PFcombErr_branch->GetEntry(index);
    } else {
      printf("branch ele_PFcombErr_branch does not exist!\n");
      exit(1);
    }
    ele_PFcombErr_isLoaded = true;
  }
  return ele_PFcombErr_;
}
const int &CMS3::ele_ID1_cat() {
  if (not ele_ID1_cat_isLoaded) {
    if (ele_ID1_cat_branch != 0) {
      ele_ID1_cat_branch->GetEntry(index);
    } else {
      printf("branch ele_ID1_cat_branch does not exist!\n");
      exit(1);
    }
    ele_ID1_cat_isLoaded = true;
  }
  return ele_ID1_cat_;
}
const float &CMS3::ele_oldsigmaiphiiphi() {
  if (not ele_oldsigmaiphiiphi_isLoaded) {
    if (ele_oldsigmaiphiiphi_branch != 0) {
      ele_oldsigmaiphiiphi_branch->GetEntry(index);
    } else {
      printf("branch ele_oldsigmaiphiiphi_branch does not exist!\n");
      exit(1);
    }
    ele_oldsigmaiphiiphi_isLoaded = true;
  }
  return ele_oldsigmaiphiiphi_;
}
const float &CMS3::ele_SCfbrem() {
  if (not ele_SCfbrem_isLoaded) {
    if (ele_SCfbrem_branch != 0) {
      ele_SCfbrem_branch->GetEntry(index);
    } else {
      printf("branch ele_SCfbrem_branch does not exist!\n");
      exit(1);
    }
    ele_SCfbrem_isLoaded = true;
  }
  return ele_SCfbrem_;
}
const int &CMS3::ele_vtxconv() {
  if (not ele_vtxconv_isLoaded) {
    if (ele_vtxconv_branch != 0) {
      ele_vtxconv_branch->GetEntry(index);
    } else {
      printf("branch ele_vtxconv_branch does not exist!\n");
      exit(1);
    }
    ele_vtxconv_isLoaded = true;
  }
  return ele_vtxconv_;
}
const float &CMS3::ele_dz() {
  if (not ele_dz_isLoaded) {
    if (ele_dz_branch != 0) {
      ele_dz_branch->GetEntry(index);
    } else {
      printf("branch ele_dz_branch does not exist!\n");
      exit(1);
    }
    ele_dz_isLoaded = true;
  }
  return ele_dz_;
}
const int &CMS3::nLumi() {
  if (not nLumi_isLoaded) {
    if (nLumi_branch != 0) {
      nLumi_branch->GetEntry(index);
    } else {
      printf("branch nLumi_branch does not exist!\n");
      exit(1);
    }
    nLumi_isLoaded = true;
  }
  return nLumi_;
}
const int &CMS3::ele_echarge() {
  if (not ele_echarge_isLoaded) {
    if (ele_echarge_branch != 0) {
      ele_echarge_branch->GetEntry(index);
    } else {
      printf("branch ele_echarge_branch does not exist!\n");
      exit(1);
    }
    ele_echarge_isLoaded = true;
  }
  return ele_echarge_;
}
const int &CMS3::ele_ID1_pass() {
  if (not ele_ID1_pass_isLoaded) {
    if (ele_ID1_pass_branch != 0) {
      ele_ID1_pass_branch->GetEntry(index);
    } else {
      printf("branch ele_ID1_pass_branch does not exist!\n");
      exit(1);
    }
    ele_ID1_pass_isLoaded = true;
  }
  return ele_ID1_pass_;
}
const bool &CMS3::ele_isEEDeeGap() {
  if (not ele_isEEDeeGap_isLoaded) {
    if (ele_isEEDeeGap_branch != 0) {
      ele_isEEDeeGap_branch->GetEntry(index);
    } else {
      printf("branch ele_isEEDeeGap_branch does not exist!\n");
      exit(1);
    }
    ele_isEEDeeGap_isLoaded = true;
  }
  return ele_isEEDeeGap_;
}
const float &CMS3::ele_conversionVertexFitProbability() {
  if (not ele_conversionVertexFitProbability_isLoaded) {
    if (ele_conversionVertexFitProbability_branch != 0) {
      ele_conversionVertexFitProbability_branch->GetEntry(index);
    } else {
      printf("branch ele_conversionVertexFitProbability_branch does not exist!\n");
      exit(1);
    }
    ele_conversionVertexFitProbability_isLoaded = true;
  }
  return ele_conversionVertexFitProbability_;
}
const float &CMS3::ele_olde25max() {
  if (not ele_olde25max_isLoaded) {
    if (ele_olde25max_branch != 0) {
      ele_olde25max_branch->GetEntry(index);
    } else {
      printf("branch ele_olde25max_branch does not exist!\n");
      exit(1);
    }
    ele_olde25max_isLoaded = true;
  }
  return ele_olde25max_;
}
const bool &CMS3::ele_ecalDrivenSeed() {
  if (not ele_ecalDrivenSeed_isLoaded) {
    if (ele_ecalDrivenSeed_branch != 0) {
      ele_ecalDrivenSeed_branch->GetEntry(index);
    } else {
      printf("branch ele_ecalDrivenSeed_branch does not exist!\n");
      exit(1);
    }
    ele_ecalDrivenSeed_isLoaded = true;
  }
  return ele_ecalDrivenSeed_;
}
const float &CMS3::ele_eta() {
  if (not ele_eta_isLoaded) {
    if (ele_eta_branch != 0) {
      ele_eta_branch->GetEntry(index);
    } else {
      printf("branch ele_eta_branch does not exist!\n");
      exit(1);
    }
    ele_eta_isLoaded = true;
  }
  return ele_eta_;
}
const float &CMS3::ele_oldhe() {
  if (not ele_oldhe_isLoaded) {
    if (ele_oldhe_branch != 0) {
      ele_oldhe_branch->GetEntry(index);
    } else {
      printf("branch ele_oldhe_branch does not exist!\n");
      exit(1);
    }
    ele_oldhe_isLoaded = true;
  }
  return ele_oldhe_;
}
const float &CMS3::ele_olde15() {
  if (not ele_olde15_isLoaded) {
    if (ele_olde15_branch != 0) {
      ele_olde15_branch->GetEntry(index);
    } else {
      printf("branch ele_olde15_branch does not exist!\n");
      exit(1);
    }
    ele_olde15_isLoaded = true;
  }
  return ele_olde15_;
}
const int &CMS3::mc_ele_matchMother_PDGID() {
  if (not mc_ele_matchMother_PDGID_isLoaded) {
    if (mc_ele_matchMother_PDGID_branch != 0) {
      mc_ele_matchMother_PDGID_branch->GetEntry(index);
    } else {
      printf("branch mc_ele_matchMother_PDGID_branch does not exist!\n");
      exit(1);
    }
    mc_ele_matchMother_PDGID_isLoaded = true;
  }
  return mc_ele_matchMother_PDGID_;
}
const float &CMS3::ele_ep() {
  if (not ele_ep_isLoaded) {
    if (ele_ep_branch != 0) {
      ele_ep_branch->GetEntry(index);
    } else {
      printf("branch ele_ep_branch does not exist!\n");
      exit(1);
    }
    ele_ep_isLoaded = true;
  }
  return ele_ep_;
}
const float &CMS3::ele_deltaetaseed() {
  if (not ele_deltaetaseed_isLoaded) {
    if (ele_deltaetaseed_branch != 0) {
      ele_deltaetaseed_branch->GetEntry(index);
    } else {
      printf("branch ele_deltaetaseed_branch does not exist!\n");
      exit(1);
    }
    ele_deltaetaseed_isLoaded = true;
  }
  return ele_deltaetaseed_;
}
const float &CMS3::ele_conv_dist() {
  if (not ele_conv_dist_isLoaded) {
    if (ele_conv_dist_branch != 0) {
      ele_conv_dist_branch->GetEntry(index);
    } else {
      printf("branch ele_conv_dist_branch does not exist!\n");
      exit(1);
    }
    ele_conv_dist_isLoaded = true;
  }
  return ele_conv_dist_;
}
const int &CMS3::ele_expected_inner_hits() {
  if (not ele_expected_inner_hits_isLoaded) {
    if (ele_expected_inner_hits_branch != 0) {
      ele_expected_inner_hits_branch->GetEntry(index);
    } else {
      printf("branch ele_expected_inner_hits_branch does not exist!\n");
      exit(1);
    }
    ele_expected_inner_hits_isLoaded = true;
  }
  return ele_expected_inner_hits_;
}
const int &CMS3::mc_ele_matchMother_status() {
  if (not mc_ele_matchMother_status_isLoaded) {
    if (mc_ele_matchMother_status_branch != 0) {
      mc_ele_matchMother_status_branch->GetEntry(index);
    } else {
      printf("branch mc_ele_matchMother_status_branch does not exist!\n");
      exit(1);
    }
    mc_ele_matchMother_status_isLoaded = true;
  }
  return mc_ele_matchMother_status_;
}
void CMS3::progress( int nEventsTotal, int nEventsChain ){
  int period = 1000;
  if(nEventsTotal%1000 == 0) {
    if (isatty(1)) {
      if( ( nEventsChain - nEventsTotal ) > period ){
        float frac = (float)nEventsTotal/(nEventsChain*0.01);
        printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
               "\033[0m\033[32m <---\033[0m\015", frac);
        fflush(stdout);
      }
      else {
        printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
               "\033[0m\033[32m <---\033[0m\015", 100.);
        cout << endl;
      }
    }
  }
}
namespace tas {
  const float &ele_full5x5_hcalOverEcal() { return cms3.ele_full5x5_hcalOverEcal(); }
  const int &MC_TrueNumInteractions() { return cms3.MC_TrueNumInteractions(); }
  const float &ele_fbrem() { return cms3.ele_fbrem(); }
  const float &ele_oldsirir() { return cms3.ele_oldsirir(); }
  const int &ele_N() { return cms3.ele_N(); }
  const float &mc_gen_pt() { return cms3.mc_gen_pt(); }
  const float &mc_event_weight() { return cms3.mc_event_weight(); }
  const int &ele_kfhits() { return cms3.ele_kfhits(); }
  const bool &ele_trackerDrivenSeed() { return cms3.ele_trackerDrivenSeed(); }
  const int &mc_ele_matchedFromCB2() { return cms3.mc_ele_matchedFromCB2(); }
  const vector<float> &rhs_e() { return cms3.rhs_e(); }
  const float &ele_eelepout() { return cms3.ele_eelepout(); }
  const bool &ele_isEBEEGap() { return cms3.ele_isEBEEGap(); }
  const float &ele_oldsigmaietaieta() { return cms3.ele_oldsigmaietaieta(); }
  const int &nEvent() { return cms3.nEvent(); }
  const int &ele_nbrem() { return cms3.ele_nbrem(); }
  const int &ele_index() { return cms3.ele_index(); }
  const float &scl_phi() { return cms3.scl_phi(); }
  const bool &ele_isEBPhiGap() { return cms3.ele_isEBPhiGap(); }
  const bool &ele_isEE() { return cms3.ele_isEE(); }
  const bool &ele_isEB() { return cms3.ele_isEB(); }
  const float &ele_SIP() { return cms3.ele_SIP(); }
  const float &ele_sclphiwidth() { return cms3.ele_sclphiwidth(); }
  const int &mc_ele_matchedFromCB() { return cms3.mc_ele_matchedFromCB(); }
  const int &ele_ID2_cat() { return cms3.ele_ID2_cat(); }
  const float &scl_E() { return cms3.scl_E(); }
  const int &ele_eClass() { return cms3.ele_eClass(); }
  const float &mc_gen_eta() { return cms3.mc_gen_eta(); }
  const float &ele_sclRawE() { return cms3.ele_sclRawE(); }
  const int &PU_N() { return cms3.PU_N(); }
  const float &ele_kfchi2() { return cms3.ele_kfchi2(); }
  const float &ele_hadronicOverEm() { return cms3.ele_hadronicOverEm(); }
  const bool &is_signal() { return cms3.is_signal(); }
  const float &ele_gsfchi2() { return cms3.ele_gsfchi2(); }
  const float &ele_oldcircularity() { return cms3.ele_oldcircularity(); }
  const vector<int> &rhs_iphi() { return cms3.rhs_iphi(); }
  const int &ele_sclNclus() { return cms3.ele_sclNclus(); }
  const float &ele_IoEmIop() { return cms3.ele_IoEmIop(); }
  const float &ele_pt() { return cms3.ele_pt(); }
  const int &seed_iphi() { return cms3.seed_iphi(); }
  const int &vtx_N() { return cms3.vtx_N(); }
  const float &ele_ID2() { return cms3.ele_ID2(); }
  const int &ele_valid_hits() { return cms3.ele_valid_hits(); }
  const float &ele_olde55() { return cms3.ele_olde55(); }
  const float &ele_oldr9() { return cms3.ele_oldr9(); }
  const bool &ele_isEERingGap() { return cms3.ele_isEERingGap(); }
  const float &ele_deltaetain() { return cms3.ele_deltaetain(); }
  const int &ele_lost_hits() { return cms3.ele_lost_hits(); }
  const int &mc_pho_mother_status() { return cms3.mc_pho_mother_status(); }
  const float &ele_oldhebc() { return cms3.ele_oldhebc(); }
  const float &ele_conv_radius() { return cms3.ele_conv_radius(); }
  const float &ele_scletawidth() { return cms3.ele_scletawidth(); }
  const float &ele_ID1() { return cms3.ele_ID1(); }
  const float &ele_conv_dcot() { return cms3.ele_conv_dcot(); }
  const int &mc_pho_mother_id() { return cms3.mc_pho_mother_id(); }
  const vector<int> &rhs_ieta() { return cms3.rhs_ieta(); }
  const float &scl_eta() { return cms3.scl_eta(); }
  const int &nRun() { return cms3.nRun(); }
  const bool &ele_isEBEtaGap() { return cms3.ele_isEBEtaGap(); }
  const int &seed_ieta() { return cms3.seed_ieta(); }
  const int &ele_gsfhits() { return cms3.ele_gsfhits(); }
  const float &mc_gen_ID() { return cms3.mc_gen_ID(); }
  const float &ele_3x3() { return cms3.ele_3x3(); }
  const float &ele_psEoverEraw() { return cms3.ele_psEoverEraw(); }
  const int &ele_ID2_pass() { return cms3.ele_ID2_pass(); }
  const float &scl_Et() { return cms3.scl_Et(); }
  const float &ele_combErr() { return cms3.ele_combErr(); }
  const float &ele_deltaphiin() { return cms3.ele_deltaphiin(); }
  const float &seed_e() { return cms3.seed_e(); }
  const float &ele_trackErr() { return cms3.ele_trackErr(); }
  const float &ele_PFcombErr() { return cms3.ele_PFcombErr(); }
  const int &ele_ID1_cat() { return cms3.ele_ID1_cat(); }
  const float &ele_oldsigmaiphiiphi() { return cms3.ele_oldsigmaiphiiphi(); }
  const float &ele_SCfbrem() { return cms3.ele_SCfbrem(); }
  const int &ele_vtxconv() { return cms3.ele_vtxconv(); }
  const float &ele_dz() { return cms3.ele_dz(); }
  const int &nLumi() { return cms3.nLumi(); }
  const int &ele_echarge() { return cms3.ele_echarge(); }
  const int &ele_ID1_pass() { return cms3.ele_ID1_pass(); }
  const bool &ele_isEEDeeGap() { return cms3.ele_isEEDeeGap(); }
  const float &ele_conversionVertexFitProbability() { return cms3.ele_conversionVertexFitProbability(); }
  const float &ele_olde25max() { return cms3.ele_olde25max(); }
  const bool &ele_ecalDrivenSeed() { return cms3.ele_ecalDrivenSeed(); }
  const float &ele_eta() { return cms3.ele_eta(); }
  const float &ele_oldhe() { return cms3.ele_oldhe(); }
  const float &ele_olde15() { return cms3.ele_olde15(); }
  const int &mc_ele_matchMother_PDGID() { return cms3.mc_ele_matchMother_PDGID(); }
  const float &ele_ep() { return cms3.ele_ep(); }
  const float &ele_deltaetaseed() { return cms3.ele_deltaetaseed(); }
  const float &ele_conv_dist() { return cms3.ele_conv_dist(); }
  const int &ele_expected_inner_hits() { return cms3.ele_expected_inner_hits(); }
  const int &mc_ele_matchMother_status() { return cms3.mc_ele_matchMother_status(); }
}
