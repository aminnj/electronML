import glob
import os
import numpy as np

def load_data(inputdir="outputs/", prefix="flip_", nfiles=-1, nevents=-1):
    list_xmva_data = []
    list_x_data = []
    list_y_data = []
    tot_events = 0
    nfilesopen = 0
    ifiles = map(lambda x: int(x.rsplit(".",1)[0].rsplit("_",1)[1]), glob.glob(inputdir+"*xdata*npa"))
    for ifile in ifiles:
        fname_mva = "{}/{}dump_mvadata_{}.npa".format(inputdir,prefix,ifile)
        fname_x = "{}/{}dump_xdata_{}.npa".format(inputdir,prefix,ifile)
        fname_y = "{}/{}dump_ydata_{}.npa".format(inputdir,prefix,ifile)
        if not os.path.isfile(fname_mva): continue
        if not os.path.isfile(fname_x): continue
        if not os.path.isfile(fname_y): continue
        print("Loading {},{},{}".format(fname_mva,fname_x,fname_y))
        xf_mva = np.load(fname_mva)
        xf_x = np.load(fname_x)
        xf_y = np.load(fname_y)
        list_xmva_data.append(xf_mva)
        list_x_data.append(xf_x)
        list_y_data.append(xf_y)
        tot_events += len(xf_y)
        nfilesopen += 1
        if nfiles > 0 and nfilesopen >= nfiles:
            break
        if nevents > 0 and tot_events >= nevents:
            break
    xmva_data = np.concatenate(list_xmva_data)
    x_data = np.concatenate(list_x_data)
    y_data = np.concatenate(list_y_data)
    if nevents > 0:
        xmva_data = xmva_data[:nevents]
        x_data = x_data[:nevents]
        y_data = y_data[:nevents]
    print("Loaded {} files with {} events".format(nfilesopen,tot_events))
    return xmva_data, x_data, y_data

def get_ptweight_info(pt,eta,truth,do_print=False):
    """
    takes 1d array of pTs, etas, 1d array of truth info (signal=1, bkg=0)
    returns 3 arrays of bin edges for pt reweighting, scale factors (= signal/background)
    """
    binedges = np.array(range(0,80,5)+range(80,200,20)+range(200,300,50)+range(300,500,100)+range(500,1000+1,500))

    etarange1 = (np.abs(eta) > 0.)    & (np.abs(eta) < 0.8)
    etarange2 = (np.abs(eta) > 0.8)   & (np.abs(eta) < 1.479)
    etarange3 = (np.abs(eta) > 1.479) & (np.abs(eta) < 2.5)

    contents_sig = np.histogram2d(pt[(truth == 1)], np.abs(eta[(truth == 1)]), bins=[binedges, np.array([0.,0.8,1.479,2.5])], normed=True)[0]
    contents_bg = np.histogram2d(pt[(truth == 0)], np.abs(eta[(truth == 0)]), bins=[binedges, np.array([0.,0.8,1.479,2.5])], normed=True)[0]
    etasfs = 1.0*contents_sig/contents_bg

    sf1 = etasfs[:,0]
    sf2 = etasfs[:,1]
    sf3 = etasfs[:,2]

    if do_print:
        for triplet in zip(binedges[:-1],sf1,sf2,sf3):
            print("{}\t{:.4f}\t{:.4f}\t{:.4f}".format(*triplet))

    return binedges, (sf1,sf2,sf3)

def get_ptweight(pt,eta,truth,binedges,sfs):
    """
    takes 1d array of pTs, etas, 1d array of truth info (signal=1, bkg=0),
    1d array of binedges, scale factors, and errors (from get_ptweight_info)
    returns array of weights (1 for signal)
    """
    sfs1, sfs2, sfs3 = sfs

    etarange1 = (np.abs(eta) > 0.)    & (np.abs(eta) < 0.8)
    etarange2 = (np.abs(eta) > 0.8)   & (np.abs(eta) < 1.479)
    etarange3 = (np.abs(eta) > 1.479) & (np.abs(eta) < 2.5)

    weights = np.ones(len(pt))
    for edge,sf1,sf2,sf3 in zip(binedges,sfs1,sfs2,sfs3):
        # print edge, sf, err
        which = (pt>=edge) & (truth==0) & etarange1
        weights[which] = sf1
        which = (pt>=edge) & (truth==0) & etarange2
        weights[which] = sf2
        which = (pt>=edge) & (truth==0) & etarange3
        weights[which] = sf3
    return weights


def train_test_split(*args,**kwargs):
    train_size = 1.-kwargs.get("test_size", 0.5)
    for arg in args:
        n_total = arg.shape[0]
        n_train = int(train_size*n_total)
        train = arg[:n_train]
        test = arg[n_train-n_total:]
        yield train
        yield test
