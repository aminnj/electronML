import ROOT as r
import numpy as np
from tqdm import tqdm


ch = r.TChain("tree")
ch.Add("output_total.root")
# for i in range(1,10):
#     ch.Add("/hadoop/cms/store/user/namin/ProjectMetis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2_MINIAODSIM_v5/output_{}.root".format(i))

xdata = []
mvadata = []
ydata = []

IMAGE_WIDTH = 29 # x - phi
IMAGE_HEIGHT = 15 # y - eta
# print rhs_ieta, rhs_iphi

def get_image(iphis,ietas,energies):
    image_half_phi = (IMAGE_WIDTH - 1)//2
    image_half_eta = (IMAGE_HEIGHT - 1)//2

    image1 = np.zeros((IMAGE_HEIGHT,IMAGE_WIDTH))
    for iphi,ieta,energy in zip(iphis,ietas,energies):
        if abs(ieta) > image_half_eta: continue
        if abs(iphi) > image_half_phi: continue
        image1[(ieta+image_half_eta,iphi+image_half_phi)] = energy

    # for L2 norm
    # norm = np.sqrt(np.sum(self.image*self.image))
    # self.image /= norm
    # self.image = np.c_[self.image1,self.image2]
    return image1

def print_image(image):
    buff = ""
    for row in image:
        for col in row:
            if col > 0.1:
                buff += "#"
            elif col > 0.01:
                buff += "*"
            elif col > 0.001:
                buff += "o"
            elif col > 0.0001:
                buff += "."
            else:
                buff += "`"
            buff += "  "
        buff += "\n"
    print buff


ievent = 0
do_extra = True
ch.SetBranchStatus("*",0)
for var in [
        "ele_eta",
        "ele_pt",
        "mc_ele_matchedFromCB",
        "ele_ID1",
        "seed_ieta",
        "seed_iphi",
        "rhs_ieta",
        "rhs_iphi",
        "rhs_e",
        
        "ele_kfhits",
        "ele_oldsigmaietaieta",
        "ele_oldsigmaiphiiphi",
        "ele_oldcircularity",
        "ele_oldr9",
        "ele_scletawidth",
        "ele_sclphiwidth",
        "ele_oldhe",
        "ele_psEoverEraw",
        "ele_kfchi2",
        "ele_gsfchi2",
        "ele_fbrem",
        "ele_ep",
        "ele_eelepout",
        "ele_IoEmIop",
        "ele_deltaetain",
        "ele_deltaphiin",
        "ele_deltaetaseed",
        "scl_eta",
        "ele_gsfhits",
        "ele_expected_inner_hits",
        "ele_conversionVertexFitProbability",

        ]:

    ch.SetBranchStatus(var,1)
for event in tqdm(ch, total=ch.GetEntries()):
    eta = event.ele_eta

    ievent += 1
    # if ievent > 2000000: break;

    if abs(eta) > 1.4: continue # FIXME

    pt = event.ele_pt
    if pt < 5.: continue
    mva = event.ele_ID1
    ismatch = -1


    theircode = event.mc_ele_matchedFromCB
    if theircode == 0: # bkg unmatch
        ismatch = 0
        matchtype = 0
    elif theircode == 1: # sig
        ismatch = 1
        matchtype = 2
    elif theircode == 3: # bkg from b,c
        ismatch = 0
        matchtype = 1
    else:
        print "shouldn't get here with theircode =",theircode

    if do_extra:
        scl_eta_ = event.scl_eta

        ele_kfhits_ = event.ele_kfhits
        ele_oldsigmaietaieta_ = event.ele_oldsigmaietaieta
        ele_oldsigmaiphiiphi_ = event.ele_oldsigmaiphiiphi
        ele_oldcircularity_ = max(min(event.ele_oldcircularity, 2.), -1.)
        ele_oldr9_ = min(event.ele_oldr9, 5.)
        ele_scletawidth_ = event.ele_scletawidth
        ele_sclphiwidth_ = event.ele_sclphiwidth
        ele_oldhe_ = event.ele_oldhe
        ele_psEoverEraw_ = event.ele_psEoverEraw
        ele_kfchi2_ = min(event.ele_kfchi2, 10.)
        ele_chi2_hits_ = min(event.ele_gsfchi2, 200.)
        ele_fbrem_ = max(event.ele_fbrem, -1.)
        ele_ep_ = min(event.ele_ep, 20.)
        ele_eelepout_ = min(event.ele_eelepout, 20.)
        ele_IoEmIop_ = event.ele_IoEmIop
        ele_deltaetain_ = min(abs(event.ele_deltaetain), 0.06)
        ele_deltaphiin_ = min(abs(event.ele_deltaphiin), 0.6)
        ele_deltaetaseed_ = min(abs(event.ele_deltaetaseed), 0.2)
        ele_gsfhits_ = event.ele_gsfhits
        ele_expectedMissingInnerHits_ = event.ele_expected_inner_hits
        ele_convVtxFitProbability_ = event.ele_conversionVertexFitProbability
        if abs(ele_oldsigmaiphiiphi_) > 1e9: ele_oldsigmaiphiiphi_ = 0.

    seed_ieta = event.seed_ieta
    seed_iphi = event.seed_iphi
    seed_e = event.seed_iphi

    if seed_e < 1e-6: continue

    # rhs_iphi = np.array(event.rhs_iphi)
    # rhs_ieta = np.array(event.rhs_ieta)
    # rhs_iphi[rhs_iphi < -180] += 360 # correct edge effects
    # rhs_iphi[rhs_iphi > 180] -= 360 # correct edge effects
    # rhs_e = np.array(event.rhs_e)
    # rhs_iphi -= seed_iphi
    # rhs_ieta -= seed_ieta
    # rhs_e /= np.max(rhs_e)
    # image = get_image(rhs_iphi,rhs_ieta,rhs_e)
    # xdata.append(image)


    mvarow = [
        ele_oldsigmaietaieta_, # cluster shape
        ele_oldsigmaiphiiphi_, # cluster shape
        ele_oldcircularity_, # cluster shape
        ele_oldr9_, # cluster shape
        ele_scletawidth_, # cluster shape
        ele_sclphiwidth_, # cluster shape
        ele_kfhits_,
        ele_oldhe_,
        ele_psEoverEraw_,
        ele_kfchi2_,
        ele_chi2_hits_,
        ele_fbrem_,
        ele_ep_,
        ele_eelepout_,
        ele_IoEmIop_,
        ele_deltaetain_,
        ele_deltaphiin_,
        ele_deltaetaseed_,
        ele_gsfhits_,
        ele_expectedMissingInnerHits_,
        ele_convVtxFitProbability_,
        ]

    mvadata.append(mvarow)
    ydata.append([matchtype, pt, eta, mva])


np.array(mvadata, dtype=np.float32).dump("dump_mvadata.npa")
# np.array(xdata, dtype=np.float32).dump("dump_xdata.npa")
np.array(ydata, dtype=np.float32).dump("dump_ydata.npa")
