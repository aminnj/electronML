
'''Trains a simple convnet on the MNIST dataset.

Gets to 99.25% test accuracy after 12 epochs
(there is still a lot of margin for parameter tuning).
16 seconds per epoch on a GRID K520 GPU.
'''

from __future__ import print_function
import numpy as np
np.random.seed(1337)  # for reproducibility
import sys
import os
import socket

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"

from keras import backend as K
if "ucsd" in socket.gethostname():
    use_tf = False
    K.set_image_dim_ordering('th')
else:
    use_tf = True

from keras.datasets import mnist
# https://keras.io/getting-started/faq/#how-can-i-save-a-keras-model
from keras.models import Sequential, load_model 
from keras.layers.core import Dense, Dropout, Activation, Flatten
from keras.layers.convolutional import Convolution2D, MaxPooling2D
from keras.utils import np_utils

import xgboost as xgb

import utils

from sklearn.metrics import roc_auc_score

batch_size = 512
num_classes = 2
# epochs = 15
epochs = 10
img_rows, img_cols = 15, 29
nb_filters = 32
nb_pool = 2
nb_conv = 5
# if load_from is not None, we always load from model files for keras and xgb
# if load_from is None, we re-train keras/xgb and then if save_to is not None, we save the models
save_to = "model.h5"
load_from = "model.h5"
# load_from = None

# xmva_data, x_data, y_data = utils.load_data(inputdir="outputs/",prefix="flip_",nfiles=35)
xmva_data, x_data, y_data = utils.load_data(inputdir="outputs/",prefix="flip_",nevents=10000)
# xmva_data, x_data, y_data = utils.load_data(inputdir="outputs/",prefix="flip_",nevents=100000)
print("Loaded data")

truth, extra = y_data[:,0], y_data[:,range(1,y_data.shape[1])]
pt = extra[:,0]
eta = extra[:,1]
mva = extra[:,2]

# matchtype == 2 is our signal, so screw everything else
truth[truth != 2] = 0
truth[truth == 2] = 1


binedges, sfs = utils.get_ptweight_info(pt,eta,truth,do_print=False)
weights = utils.get_ptweight(pt,eta,truth,binedges,sfs)
print("Got the pt/eta weights")


x_train, x_test, xmva_train, xmva_test, y_train, y_test, extra_train, extra_test, weights_train, weights_test = utils.train_test_split(x_data, xmva_data, truth, extra, weights, test_size=0.3, random_state=43)
sieie_test = xmva_test[:,0]
print("Did splitting")

xshape_train = xmva_train[:,range(6)]
xshape_test = xmva_test[:,range(6)]
print(xshape_train)
print(y_train)
dtrain = xgb.DMatrix( xshape_train, label=y_train, weight=np.abs(weights_train))
dtest = xgb.DMatrix( xshape_test, label=y_test, weight=np.abs(weights_test))
ihalf = len(xshape_test)
evallist  = [(dtrain,'train'), (dtest,'eval')]
num_round = 300
param = {}
param['objective'] = 'binary:logistic'
param['max_depth'] = 3
param['silent'] = 1
param['nthread'] = 12
param['eval_metric'] = "auc"
if load_from:
    bst = xgb.Booster(model_file=load_from.replace("h5","xgb"))
else:
    bst = xgb.train( param.items(), dtrain, num_round, evallist, early_stopping_rounds=50 )
    featurenames = ["ele_oldsigmaietaieta_", "ele_oldsigmaiphiiphi_", "ele_oldcircularity_", "ele_oldr9_", "ele_scletawidth_", "ele_sclphiwidth_"]
    print("XGB feature rankings:")
    for rank,name in sorted(zip(map(bst.get_fscore().get, bst.feature_names),featurenames), reverse=True):
        print("  {:4s} {:15s}".format(str(rank),str(name)))

    if save_to:
        bst.save_model(save_to.replace("h5","xgb"))
print("Predicting with xgb")
xgb_y_pred = bst.predict(dtest)


print("Counts:")
print("  train ntot",len(x_train))
print("  train nsig",(y_train==1).sum())
print("  train nbkg",(y_train==0).sum())
print("  test ntot",len(x_test))
print("  test nsig",(y_test==1).sum())
print("  test nbkg",(y_test==0).sum())

xmva_train = xmva_train.reshape(xmva_train.shape[0], len(xmva_train[0]))
xmva_test = xmva_test.reshape(xmva_test.shape[0], len(xmva_test[0]))

if use_tf:
    x_train = x_train.reshape(x_train.shape[0], img_rows, img_cols, 1)
    x_test = x_test.reshape(x_test.shape[0], img_rows, img_cols, 1)
    input_shape = (img_rows, img_cols, 1)
else:
    x_train = x_train.reshape(x_train.shape[0], 1, img_rows, img_cols)
    x_test = x_test.reshape(x_test.shape[0], 1, img_rows, img_cols)
    input_shape = (1, img_rows, img_cols)


x_train = x_train.astype('float32')
x_test = x_test.astype('float32')

# convert class vectors to binary class matrices
y_train = np_utils.to_categorical(y_train, num_classes)
y_test = np_utils.to_categorical(y_test, num_classes)

if load_from:
    print("Loading CNN model")
    model = load_model(load_from)
else:

    model = Sequential()
    model.add(Convolution2D(nb_filters, nb_conv, nb_conv,
                            border_mode='valid',
                            input_shape = input_shape))
    model.add(Activation('relu'))
    model.add(Convolution2D(nb_filters*2, nb_conv, nb_conv))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(nb_pool, nb_pool)))
    model.add(Dropout(0.2))
    model.add(Flatten())
    model.add(Dense(1024))
    model.add(Activation('relu'))
    model.add(Dense(256))
    model.add(Activation('relu'))
    model.add(Dropout(0.3))
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

    if save_to:
        model.save(save_to)

print("Predicting")
y_pred = model.predict(x_test)
# score = model.evaluate(x_test, y_test, verbose=0)
# print('Test loss:', score[0])
# print('Test accuracy:', score[1])

# test:         y_test, y_pred, pt, mva, weights sigmaietaieta
print("Dumping extra data")
todump = np.c_[ y_test[:,1], y_pred[:,1], extra_test[:,0], extra_test[:,2], weights_test, sieie_test, xgb_y_pred ]
np.array(todump, dtype=np.float32).dump("todump.npa")

print("AUC total (CNN):",roc_auc_score(y_test[:,1],y_pred[:,1]))
print("AUC total (XGB):",roc_auc_score(y_test[:,1],xgb_y_pred))
print("AUC total (SIEIE):",roc_auc_score(y_test[:,1],1.-sieie_test))
print("AUC total (CMS):",roc_auc_score(y_test[:,1],extra_test[:,2]))
