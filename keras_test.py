from keras.models import Sequential
from keras.layers.convolutional import Conv3D
from keras.layers import Conv3D, MaxPooling3D
from keras.layers.normalization import BatchNormalization
import numpy as np
import pylab as plt



seq = Sequential()
seq.add(Conv3D(filters=1, kernel_size=(120,120,120),
               activation='sigmoid',
               padding='same', data_format='channels_last'))


seq.add(MaxPooling3D(pool_size=(3,3,3),strides=(2,2,2))

seq.add(Dense(activation='relu')



seq.add(BatchNormalization())
