import numba
import numpy as np
import uproot
from tqdm import tqdm
import os
import glob
from multiprocessing import Pool as ThreadPool


IMAGE_WIDTH = 29  # x - phi
IMAGE_HEIGHT = 15  # y - eta
@numba.jit(nopython=True)
def get_image(iphis,ietas,energies, seed_iphi, seed_ieta, charge):
    # center on seed
    iphis -= seed_iphi
    ietas -= seed_ieta
    # correct for wrap around effects in phi
    iphis[iphis < -180] += 360
    iphis[iphis > 180] -= 360

    # flip negatively charged e, so that bend goes in same direction
    if charge < 0: iphis *= -1

    image_half_phi = (IMAGE_WIDTH - 1)//2
    image_half_eta = (IMAGE_HEIGHT - 1)//2
    image1 = np.zeros((IMAGE_HEIGHT,IMAGE_WIDTH))
    for iphi,ieta,energy in zip(iphis,ietas,energies):
        if abs(ieta) > image_half_eta: continue
        if abs(iphi) > image_half_phi: continue
        image1[(ieta+image_half_eta,iphi+image_half_phi)] = energy
    return image1

def clamp(arr, lower=None, upper=None):
    if lower is not None:
        arr[arr < lower] = lower
    if upper is not None:
        arr[arr > upper] = upper
    return arr

def convert_file(fname, outputdir="outputs/", prefix="", suffix="", treename="tree"):
    if not os.path.exists(fname):
        raise Exception("File does not exist: {}".format(fname))

    if not os.path.exists(outputdir):
        os.system("mkdir -p {}".format(outputdir))

    f = uproot.open(fname)
    t = f[treename]

    ele_oldsigmaietaieta_, ele_oldsigmaiphiiphi_, \
        ele_oldcircularity_, ele_oldr9_, \
        ele_scletawidth_, ele_sclphiwidth_, \
        ele_kfhits_, ele_oldhe_, \
        ele_psEoverEraw_, ele_kfchi2_, \
        ele_chi2_hits_, ele_fbrem_, \
        ele_ep_, ele_eelepout_, \
        ele_IoEmIop_, ele_deltaetain_, \
        ele_deltaphiin_, ele_deltaetaseed_, \
        ele_gsfhits_, ele_expectedMissingInnerHits_, \
        ele_convVtxFitProbability_, \
        mvas, seed_ietas, seed_iphis, seed_es, \
        charges, rhs_iphis, rhs_ietas, rhs_es, \
        pts, etas, truenumints, theircode = \
        t.arrays([
                "ele_oldsigmaietaieta", "ele_oldsigmaiphiiphi",
                "ele_oldcircularity", "ele_oldr9",
                "ele_scletawidth", "ele_sclphiwidth",
                "ele_kfchi2", "ele_oldhe",
                "ele_psEoverEraw", "ele_kfchi2",
                "ele_gsfchi2", "ele_fbrem",
                "ele_ep", "ele_eelepout",
                "ele_IoEmIop", "ele_deltaetain",
                "ele_deltaphiin", "ele_deltaetaseed",
                "ele_gsfhits", "ele_expected_inner_hits",
                "ele_conversionVertexFitProbability",
                "ele_ID1", "seed_ieta", "seed_iphi",
                "seed_e", "ele_echarge", "rhs_iphi",
                "rhs_ieta", "rhs_e", "ele_pt",
                "ele_eta", "MC_TrueNumInteractions",
                "mc_ele_matchedFromCB",
                ], outputtype=tuple)

    # 0, 1, 3 = bkg unmatch, sig, bkg from b/c
    matchtype = np.copy(theircode)
    # remap from their convention to mine
    matchtype[matchtype == 0] = 0
    matchtype[matchtype == 1] = 2
    matchtype[matchtype == 3] = 1

    ele_deltaetain_ = np.abs(ele_deltaetain_)
    ele_deltaphiin_ = np.abs(ele_deltaphiin_)
    ele_deltaetaseed_ = np.abs(ele_deltaetaseed_)
    ele_oldcircularity_ = clamp(ele_oldcircularity_, lower=-1, upper=2.)
    ele_oldr9_ = clamp(ele_oldr9_, upper=5.)
    ele_kfchi2_ = clamp(ele_kfchi2_, upper=10.)
    ele_chi2_hits_ = clamp(ele_chi2_hits_, upper=200.)
    ele_fbrem_ = clamp(ele_fbrem_, lower=-1.)
    ele_ep_ = clamp(ele_ep_, upper=20.)
    ele_eelepout_ = clamp(ele_eelepout_, upper=20.)
    ele_deltaetain_ = clamp(ele_deltaetain_, upper=0.06)
    ele_deltaphiin_ = clamp(ele_deltaphiin_, upper=0.6)
    ele_deltaetaseed_ = clamp(ele_deltaetaseed_, upper=0.2)
    ele_oldsigmaiphiiphi_[np.abs(ele_oldsigmaiphiiphi_) > 1e9] = 0.

    mvamat = np.c_[
            ele_oldsigmaietaieta_, ele_oldsigmaiphiiphi_,
            ele_oldcircularity_, ele_oldr9_,
            ele_scletawidth_, ele_sclphiwidth_,
            ele_kfhits_, ele_oldhe_,
            ele_psEoverEraw_, ele_kfchi2_,
            ele_chi2_hits_, ele_fbrem_,
            ele_ep_, ele_eelepout_,
            ele_IoEmIop_, ele_deltaetain_,
            ele_deltaphiin_, ele_deltaetaseed_,
            ele_gsfhits_, ele_expectedMissingInnerHits_,
            ele_convVtxFitProbability_,
            ]

    isgood = (pts > 5.) & (pts < 2.0e3) & (seed_es > 1.e-6) & (np.abs(etas) < 1.4)

    images = []
    for good, rhs_e, rhs_iphi, rhs_ieta, rhs_e, \
            seed_e, seed_iphi, seed_ieta, charge \
            in zip(isgood, rhs_es, rhs_iphis, rhs_ietas, rhs_es,
                   seed_es, seed_iphis, seed_ietas, charges):

        if not good: continue

        image = get_image(
                rhs_iphi,rhs_ieta,rhs_e,
                seed_iphi, seed_ieta,
                charge
                )
        images.append(image)

    x_data = np.array(images)
    mva_data = mvamat[isgood]
    y_data = np.c_[
            matchtype,
            pts,
            etas,
            mvas,
            seed_iphis,
            seed_ietas,
            truenumints,
            ][isgood]

    idx = int(fname.rsplit(".",1)[0].split("_")[-1])
    outfname = "{}/{}data{}_{}.npz".format(outputdir,prefix,suffix,idx)
    np.savez_compressed(outfname,
                        mva_data=np.array(mva_data,dtype=np.float32),
                        x_data=np.array(x_data,dtype=np.float32),
                        y_data=np.array(y_data,dtype=np.float32),
                        )

    return len(x_data)

if __name__ == "__main__":

    pool = ThreadPool(10)
    thedir = "/hadoop/cms/store/user/namin/elemva/ProjectMetis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2_MINIAODSIM_v10/"
    fnames = glob.glob(thedir + "/*.root")
    for N in tqdm(pool.imap_unordered(convert_file,fnames),total=len(fnames)):
        pass
        # print ">>> Processed {} events (so far)".format(N)
    pool.close()
    pool.join()
