#!/usr/bin/env python3



##  two-bidirectional-GRU-layer neural net

import time
import sys, os, re, csv, codecs, numpy as np, pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score


from keras.preprocessing.text import Tokenizer
from keras.preprocessing.sequence import pad_sequences
from keras import initializers, regularizers, constraints, optimizers, layers, callbacks
from keras import backend as K
from keras.engine import InputSpec, Layer
from keras.layers import Dense, Input, LSTM, Embedding, Dropout, Activation, Conv1D, GRU, CuDNNGRU, CuDNNLSTM
from keras.layers import Bidirectional, GlobalMaxPool1D, MaxPooling1D, Add, Flatten
from keras.layers import GlobalAveragePooling1D, GlobalMaxPooling1D, concatenate, SpatialDropout1D
from keras.layers import GRU, BatchNormalization, Conv1D, MaxPooling1D
from keras.models import Model, load_model
from keras.callbacks import Callback
from keras.callbacks import EarlyStopping, ModelCheckpoint, LearningRateScheduler
from keras.optimizers import Adam, RMSprop



import logging



start_time = time.time()
np.random.seed(5079)



class RocAucEvaluation(Callback):
    def __init__(self, validation_data=(), interval=1):
        super(Callback, self).__init__()
        self.interval = interval
        self.X_val, self.y_val = validation_data



    def on_epoch_end(self, epoch, logs={}):
        if epoch % self.interval == 0:
            y_pred = self.model.predict(self.X_val, verbose=0)
            score = roc_auc_score(self.y_val, y_pred)
            print("\n ROC-AUC -- EPOCH: {:d} - SCORE: {:.6f}".format(epoch+1, score))

            

            

training = pd.read_csv("train.csv")
testing = pd.read_csv("test.csv")



embed_size = 300
max_features = 200000
max_len = 200



list_classes = [1:6]
y = training[list_classes].values
train["CT"].fillna("NA")
test["CT"].fillna("NA")
X_train, X_valid, Y_train, Y_valid = train_test_split(training, y, test_size = 0.1)



raw_text_train = X_train["CT"].str.lower()
raw_text_valid = X_valid["CT"].str.lower()
raw_text_test = testing["CT"].str.lower()



tk = Tokenizer(num_words = max_features, lower = True)
tk.fit_on_texts(raw_text_train)
X_train["comment_seq"] = tk.texts_to_sequences(raw_text_train)
X_valid["comment_seq"] = tk.texts_to_sequences(raw_text_valid)
testing["comment_seq"] = tk.texts_to_sequences(raw_text_test)



X_train = pad_sequences(X_train.comment_seq, maxlen = max_len)
X_valid = pad_sequences(X_valid.comment_seq, maxlen = max_len)
testing = pad_sequences(test.comment_seq, maxlen = max_len)



def get_coefs(word,*arr): return word, np.asarray(arr, dtype='float32')
embedding_index = dict(get_coefs(*o.strip().split(" ")) for o in open(embedding_path))



word_index = tk.word_index
nb_words = min(max_features, len(word_index))
embedding_matrix = np.zeros((nb_words, embed_size))
for word, i in word_index.items():
    if i >= max_features: continue
    embedding_vector = embedding_index.get(word)
    if embedding_vector is not None: embedding_matrix[i] = embedding_vector

    

    


def build_model(lr = 0.0, lr_d = 0.0, units = 0, dr = 0.0):
    inp = Input(shape = (max_len,))
    x = Embedding(max_features, embed_size, weights = [embedding_matrix], trainable = False)(inp)
    x1 = SpatialDropout1D(dr)(x)

    x = Bidirectional(CuDNNGRU(units, return_sequences = True))(x1)
    x = Conv1D(64, kernel_size = 2, padding = "valid", kernel_initializer = "he_uniform")(x)
    y = Bidirectional(CuDNNLSTM(units, return_sequences = True))(x1)
    y = Conv1D(64, kernel_size = 2, padding = "valid", kernel_initializer = "he_uniform")(y)

    
    avg_pool1 = GlobalAveragePooling1D()(x)
    max_pool1 = GlobalMaxPooling1D()(x)
    avg_pool2 = GlobalAveragePooling1D()(y)
    max_pool2 = GlobalMaxPooling1D()(y)

    

    x = concatenate([avg_pool1, max_pool1, avg_pool2, max_pool2])
    x = Dense(6, activation = "sigmoid")(x)


    model = Model(inputs = inp, outputs = x)
    model.compile(loss = "binary_crossentropy", optimizer = Adam(lr = lr, decay = lr_d), metrics = ["accuracy"])
    history = model.fit(X_train, Y_train, batch_size = 128, epochs = 3, validation_data = (X_valid, Y_valid), 
                        verbose = 1, callbacks = [ra_val, check_point, early_stop])
    model = load_model(file_path)

    return model

    

    
pred = 0
n_seeds = 10

for i in range(n_seeds):
    model = build_model(lr = 1e-3, lr_d = 0, units = 128, dr = 0.2)
    pred += model.predict(testing, batch_size = 1024, verbose = 1)/n_seeds



print("[{}] Completed!".format(time.time() - start_time))

