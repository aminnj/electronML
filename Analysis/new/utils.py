import glob
import numpy as np

def load_data(inputdir="outputs/", nfiles=-1, nevents=-1, skipfiles=0):
    """
    Load data from `inputdir`, up to a maximum of `nfiles` files OR
    `nevents` events (whichever comes first), skipping the first
    `skipfiles` files
    """
    list_mva_data = []
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
        list_mva_data.append(mva_data_)
        list_x_data.append(x_data_)
        list_y_data.append(y_data_)
        tot_events += len(y_data_)
        nfilesopen += 1
        if nfiles > 0 and nfilesopen >= nfiles:
            break
        if nevents > 0 and tot_events >= nevents:
            break
        print("Loaded {}".format(fname))
    mva_data = np.concatenate(list_mva_data)
    x_data = np.concatenate(list_x_data)
    y_data = np.concatenate(list_y_data)
    if nevents > 0:
        mva_data = mva_data[:int(nevents)]
        x_data = x_data[:int(nevents)]
        y_data = y_data[:int(nevents)]
    print("Loaded {} files with {} events".format(nfilesopen,len(y_data)))
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
