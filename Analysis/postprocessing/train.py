
'''Trains a simple convnet on the MNIST dataset.

Gets to 99.25% test accuracy after 12 epochs
(there is still a lot of margin for parameter tuning).
16 seconds per epoch on a GRID K520 GPU.
'''

from __future__ import print_function
import numpy as np
np.random.seed(1337)  # for reproducibility

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
epochs = 5

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

# matchtype == 2 is our signal, so screw everything else
y_data[:,0][y_data[:,0] != 2] = 0
y_data[:,0][y_data[:,0] == 2] = 1

print(len(y_data))
print(sum(y_data[:,0]==1))
print(sum(y_data[:,0]==0))

def train_test_split(*args,**kwargs):
    test_size = kwargs.get("test_size", 0.5)
    for arg in args:
        n_total = arg.shape[0]
        n_train = int(test_size*n_total)
        train = arg[:n_train]
        test = arg[n_train-n_total:]
        yield train
        yield test

x_train, x_test, y_train, y_test, extra_train, extra_test = train_test_split(x_data, y_data[:,0], y_data[:,range(1,y_data.shape[1])], test_size=0.7, random_state=43)

x_train = x_train.reshape(x_train.shape[0], 1, img_rows, img_cols)
x_test = x_test.reshape(x_test.shape[0], 1, img_rows, img_cols)

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
          validation_data=(x_test, y_test))
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
