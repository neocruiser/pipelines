#!/usr/bin/env python3



## embedding with neural network
## keras starter NN with embeddings
## lightgbm fixing unbalanced data
## RNN with keras ridge SGDR


import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')
import os
os.environ['OMP_NUM_THREADS'] = '4'
import gc


## testing dir
path = '../input/'


train_df = pd.read_csv(path+"train.csv",l
                       dtype=dtypes,
                       skiprows = range(1, 131886954),
                       usecols=['ip',
                                'app',
                                'device',
                                'os',
                                'channel',
                                'click_time',
                                'is_attributed'])
print('load test....')
test_df = pd.read_csv(path+"test.csv",
                      dtype=dtypes,
                      usecols=['ip',
                               'app',
                               'device',
                               'os',
                               'channel',
                               'click_time',
                               'click_id'])


len_train = len(train_df)
train_df=train_df.append(test_df)
del test_df; gc.collect()


from keras.layers import Input, Embedding, Dense, Flatten, Dropout, concatenate
from keras.layers import BatchNormalization, SpatialDropout1D
from keras.callbacks import Callback
from keras.models import Model
from keras.optimizers import Adam
emb_n = 50
dense_n = 1000

batch_size = 20000
epochs = 2
exp_decay = lambda init, fin, steps: (init/fin)**(1/(steps-1)) - 1
steps = int(len(train_df) / batch_size) * epochs
lr_init, lr_fin = 0.001, 0.0001
lr_decay = exp_decay(lr_init, lr_fin, steps)
optimizer_adam = Adam(lr=0.001, decay=lr_decay)
model.compile(loss='binary_crossentropy',optimizer=optimizer_adam,metrics=['accuracy'])

model.summary()

model.fit(train_df, y_train, batch_size=batch_size, epochs=2, shuffle=True, verbose=2)
del train_df, y_train; gc.collect()
model.save_weights('dl_support.h5')

print("predicting....")
sub['is_attributed'] = model.predict(test_df, batch_size=batch_size, verbose=2)
del test_df; gc.collect()
print("writing....")
sub.to_csv('dl_support.csv',index=False)
