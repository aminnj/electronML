import time
import os
import glob
import numpy as np
from sklearn.metrics import roc_curve

def load_data(inputdir="outputs/", nfiles=-1, nevents=-1, skipfiles=0,selector=None,return_track_info=False,return_time_info=False):
    """
    Load data from `inputdir`, up to a maximum of `nfiles` files OR
    `nevents` events (whichever comes first), skipping the first
    `skipfiles` files
    `selector` is an optional lambda function applied to y_data_ returning an array
    of bools to only store selected x,y,mva data
    """
    list_mva_data = []
    list_track_data = []
    list_time_data = []
    list_chi2_data = []
    list_x_data = []
    list_y_data = []
    tot_events = 0
    nfilesopen = 0
    fnames = sorted(glob.glob("{}/*.npz".format(inputdir)))
    fnames = fnames[skipfiles:]
    for fname in fnames:
        # idx = int(x.rsplit(".",1)[0].rsplit("_",1)[1])
        data = np.load(fname)
        x_data_ = data["x_data"]
        y_data_ = data["y_data"]
        mva_data_ = data["mva_data"]
        if return_track_info:
            track_data_ = data["track_data"]
        if return_time_info:
            time_data_ = data["time_data"]
            chi2_data_ = data["chi2_data"]
        if selector:
            sel = selector(y_data_)
            x_data_ = x_data_[sel]
            y_data_ = y_data_[sel]
            mva_data_ = mva_data_[sel]
            if return_track_info:
                track_data_ = track_data_[sel]
            if return_time_info:
                time_data_ = time_data_[sel]
                chi2_data_ = chi2_data_[sel]
        list_mva_data.append(mva_data_)
        list_x_data.append(x_data_)
        list_y_data.append(y_data_)
        if return_track_info:
            list_track_data.append(track_data_)
        if return_time_info:
            list_time_data.append(time_data_)
            list_chi2_data.append(chi2_data_)
        tot_events += len(y_data_)
        nfilesopen += 1
        if nfiles > 0 and nfilesopen >= nfiles:
            break
        if nevents > 0 and tot_events >= nevents:
            break
        print("Loaded {} ({} events total)".format(fname, tot_events))
    mva_data = np.concatenate(list_mva_data)
    x_data = np.concatenate(list_x_data)
    y_data = np.concatenate(list_y_data)
    if return_track_info:
        track_data = np.concatenate(list_track_data)
    if return_time_info:
        time_data = np.concatenate(list_time_data)
        chi2_data = np.concatenate(list_chi2_data)
    if nevents > 0:
        mva_data = mva_data[:int(nevents)]
        x_data = x_data[:int(nevents)]
        y_data = y_data[:int(nevents)]
        if return_track_info:
            track_data = track_data[:int(nevents)]
        if return_time_info:
            time_data = time_data[:int(nevents)]
            chi2_data = chi2_data[:int(nevents)]
    print("Loaded {} files with {} events".format(nfilesopen,len(y_data)))
    if return_track_info:
        if return_time_info:
            return x_data, y_data, mva_data, track_data, time_data, chi2_data
        else:
            return x_data, y_data, mva_data, track_data
    else:
        return x_data, y_data, mva_data

def train_test_split(*args,**kwargs):
    train_size = 1.-kwargs.get("test_size", 0.5)
    for arg in args:
        n_total = arg.shape[0]
        n_train = int(train_size*n_total)
        train = arg[:n_train]
        test = arg[n_train-n_total:]
        yield train
        yield test

def culled_indices(thresholds, precision=0.0002):
    """
    input a list of discriminant values and a precision
    and this returns a list of indices on the values for which
    the difference wrt previous values is >= precision
    """
    prev = -999.
    to_keep = []
    for ithr,thr in enumerate(thresholds):
        diff = abs(thr - prev)
        if diff > precision:
            to_keep.append(ithr)
            prev = thr
        else:
            continue
    return np.array(to_keep)

def sigeff_for_given_bkgeff(y_true, y_pred, bkgeff, sel=None):
    """
    For a given background efficiency, this returns the
    signal efficiency, with an optional selection on the input vectors
    """
    yt,yp = y_true, y_pred
    if sel is not None:
        yt = yt[sel]
        yp = yp[sel]
    fpr, tpr, _ = roc_curve(yt,yp)
    idx = np.where(fpr > bkgeff)[0][0]
    return tpr[idx]

def get_sigeffs_for_bins(y_true,y_pred,selector,bins,bkgeff):
    """
    Takes true class labels, predictions, a list of bin edges, and
    a given background efficiency. For each bin (of the selector), computes the
    signal efficiency for that background efficiency. Returns a 3-column matrix
    with columns [bin center, signal efficiency, weight], where weight
    is sqrt(number of samples in bin)
    """
    to_ret = []
    for low,high in zip(bins[:-1],bins[1:]):
        sel = ((selector > low) & (selector < high))
        nsamps = sel.sum()
        if nsamps <= 10: continue
        sigeff = sigeff_for_given_bkgeff(y_true, y_pred, bkgeff=bkgeff, sel=sel)
        to_ret.append( [0.5*(low+high), sigeff, np.sqrt(nsamps)] )
    return np.array(to_ret)

def get_pteta_weights(pts,etas,is_signal,ptbins=None,etabins=None):
    """
    Takes arrays of pts, etas, is_signal bools, and optional pt/eta bins
    and returns an array of weights matching pts.shape with 1 for signal,
    reweighting bkg to signal
    """
    if ptbins is None:
        ptbins = np.array(range(5,78,2)+range(80,200,5)+range(200,300,10)+range(300,500,100)+range(500,1000+1,500)+[4000])
    if etabins is None:
        etabins = np.linspace(-1.4,1.4,20)
    counts_sig, _,_ = np.histogram2d(pts[is_signal],etas[is_signal],bins=[ptbins,etabins])
    counts_bkg, _,_ = np.histogram2d(pts[~is_signal],etas[~is_signal],bins=[ptbins,etabins])
    ptidxs = np.digitize(pts[~is_signal],ptbins)-1
    etaidxs = np.digitize(etas[~is_signal],etabins)-1
    weights = np.ones(len(pts))
    weights[~is_signal] = (1.0*counts_sig[ptidxs,etaidxs]/counts_bkg[ptidxs,etaidxs])
    return weights

