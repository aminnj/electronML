
'''Trains a simple convnet on the MNIST dataset.

Gets to 99.25% test accuracy after 12 epochs
(there is still a lot of margin for parameter tuning).
16 seconds per epoch on a GRID K520 GPU.
'''

from __future__ import print_function
import numpy as np
np.random.seed(1337)  # for reproducibility
import sys

from keras import backend as K
K.set_image_dim_ordering('th')

from keras.datasets import mnist
from keras.models import Sequential
from keras.layers.core import Dense, Dropout, Activation, Flatten
from keras.layers.convolutional import Convolution2D, MaxPooling2D
from keras.utils import np_utils

# from sklearn.model_selection import train_test_split
# from sklearn.metrics import roc_auc_score
# from sklearn.metrics import roc_curve

batch_size = 512
# batch_size = 1024
num_classes = 2
epochs = 1

# input image dimensions
img_rows, img_cols = 9, 23
# number of convolutional filters to use
nb_filters = 32
# size of pooling area for max pooling
nb_pool = 2
# convolution kernel size
nb_conv = 3

# the data, shuffled and split between train and test sets
# (x_train, y_train), (x_test, y_test) = mnist.load_data()
x_data = np.load("dump_xdata.npa")
y_data = np.load("dump_ydata.npa")

# print(x_data)
# print(y_data)

# print(np.sum(np.sum(x_data,axis=2),axis=1))

truth, extra = y_data[:,0], y_data[:,range(1,y_data.shape[1])]
pt = extra[:,0]
eta = extra[:,1]
mva = extra[:,2]

# matchtype == 2 is our signal, so screw everything else
truth[truth != 2] = 0
truth[truth == 2] = 1

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


binedges, sfs = get_ptweight_info(pt,eta,truth,do_print=True)
weights = get_ptweight(pt,eta,truth,binedges,sfs)
print(weights)
print("got the pt/eta weights")
# print weights, weighterrs

print(len(truth))
print(sum(truth==1))
print(sum(truth==0))


def train_test_split(*args,**kwargs):
    test_size = kwargs.get("test_size", 0.5)
    for arg in args:
        n_total = arg.shape[0]
        n_train = int(test_size*n_total)
        train = arg[:n_train]
        test = arg[n_train-n_total:]
        yield train
        yield test

x_train, x_test, y_train, y_test, extra_train, extra_test, weights_train, weights_test = train_test_split(x_data, truth, extra, weights, test_size=0.7, random_state=43)

x_train = x_train.reshape(x_train.shape[0], 1, img_rows, img_cols)
x_test = x_test.reshape(x_test.shape[0], 1, img_rows, img_cols)

# sys.exit()

x_train = x_train.astype('float32')
x_test = x_test.astype('float32')
print('x_train shape:', x_train.shape)
print(x_train.shape[0], 'train samples')
print(x_test.shape[0], 'test samples')

# convert class vectors to binary class matrices
y_train = np_utils.to_categorical(y_train, num_classes)
y_test = np_utils.to_categorical(y_test, num_classes)

model = Sequential()
model.add(Convolution2D(nb_filters, nb_conv, nb_conv,
                        border_mode='valid',
                        input_shape=(1, img_rows, img_cols)))
model.add(Activation('relu'))
model.add(Convolution2D(nb_filters, nb_conv, nb_conv))
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(nb_pool, nb_pool)))
model.add(Dropout(0.25))

model.add(Flatten())
model.add(Dense(128))
model.add(Activation('relu'))
model.add(Dropout(0.5))
model.add(Dense(num_classes))
model.add(Activation('softmax'))


model.compile(loss='categorical_crossentropy',
              optimizer='adadelta',
              metrics=['accuracy'])

model.fit(x_train, y_train,
          batch_size=batch_size,
          nb_epoch=epochs,
          verbose=1,
          sample_weight=weights_train,
          validation_data=(x_test, y_test, weights_test))
score = model.evaluate(x_test, y_test, verbose=0)
print("predicting")
y_pred = model.predict(x_test)
print('Test loss:', score[0])
print('Test accuracy:', score[1])

#      y_test, y_pred, pt, mva
todump = np.c_[ y_test[:,1], y_pred[:,1], extra_test[:,0], extra_test[:,2] ]
np.array(todump, dtype=np.float32).dump("todump.npa")
# np.array(y_test, dtype=np.float32).dump("dump_ytest.npa")
# np.array(y_pred, dtype=np.float32).dump("dump_ypred.npa")
print("AUC total",roc_auc_score(y_test[:,1],y_pred[:,1]))
