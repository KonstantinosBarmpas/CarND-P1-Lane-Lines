##-----------IMPORT LIBRARIES--------------##
import csv
import os
import cv2
import numpy as np
import sklearn
from keras.models import Model
from keras.models import Sequential
from keras.layers import Flatten, Dense, Lambda, Dropout
from keras.layers import Conv2D
from keras.layers.pooling import MaxPooling2D
from keras.layers import Cropping2D
from keras import optimizers
from sklearn.model_selection import train_test_split
from scipy import ndimage
import matplotlib.pyplot as plt
import random

#Load Images
samples = []
csv_file = '/Users/user/CarND-Behavioral-Cloning-P3-master/data_run5/driving_log.csv'
with open(csv_file, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        samples.append(row)

#Generator funtion
def generator(samples, batch_size):
    num_samples = len(samples)
    correction = 0.22
    while 1:
        for offset in range(0, num_samples, batch_size):
            batch_samples = samples[offset:offset+batch_size]
            images = []
            measurements = []
            for batch_sample in batch_samples:
                for i in range(3):
                    current_path = '/Users/user/CarND-Behavioral-Cloning-P3-master/data_run5/IMG/'+batch_sample[i].split('/')[-1]
                    image = ndimage.imread(current_path)
                    measurement = float(batch_sample[3])
                    if (i == 1):
                        measurement = measurement + correction
                    if (i == 2):
                        measurement = measurement - correction
                    images.append(image)
                    measurements.append(measurement)
            X_train = np.array(images)
            y_train = np.array(measurements)
            yield sklearn.utils.shuffle(X_train, y_train)

def NVDIA_network():
    #Dropout parameter
    dropout=0.4
    #Normalize inputs
    model=Sequential()
    model.add(Lambda(lambda x:x /255.0 - 0.5, input_shape=(160,320,3)))
    #Cropping the original image to keep only the important part
    model.add(Cropping2D(cropping=((70,25),(0,0))))
    #Layer 1 convolution
    model.add(Conv2D(24, (5, 5), activation="relu", strides=(2, 2)))
    #Layer 2 convolution
    model.add(Conv2D(36, (5, 5), activation="relu", strides=(2, 2)))
    #Dropout
    model.add(Dropout(dropout))
    #Layer 3 convolution
    model.add(Conv2D(48, (5, 5), activation="relu", strides=(2, 2)))
    #Layer 4 convolution
    model.add(Conv2D(64, (3, 3), activation="relu", data_format="channels_first"))
    #Layer 5 convolution
    model.add(Conv2D(64, (3, 3), activation="relu"))
    #Dropout
    model.add(Dropout(dropout))
    #Flatten before connected layers
    model.add(Flatten())
    #Fully Connected Layer 1
    model.add(Dense(100, activation="relu"))
    #Dropout
    model.add(Dropout(dropout))
    #Fully Connected Layer 2
    model.add(Dense(50, activation="relu"))
    #Dropout
    model.add(Dropout(dropout))
    #Fully Connected Layer 3
    model.add(Dense(10, activation="relu"))
    #Dropout
    model.add(Dense(1))
    return model;


train_samples, validation_samples = train_test_split(samples, test_size=0.2)


batch_size = 64
#Compile and train the model using the generator function
train_generator = generator(train_samples, batch_size=batch_size)
validation_generator = generator(validation_samples, batch_size=batch_size)
model=NVDIA_network()
adam = optimizers.Adam(lr=0.0001)
model.compile(loss='mse', optimizer=adam)
my_history = model.fit_generator(train_generator, samples_per_epoch= int((len(train_samples) * 3) / batch_size),
                                     validation_data=validation_generator,
                                     nb_val_samples=int((len(validation_samples) * 3) / batch_size),
                                     nb_epoch=2,
                                     verbose=1)
model.save('model.h5')

### print the keys contained in the history object
print(my_history.history.keys())

### Create the loss plot
plt.plot(my_history.history['loss'])
plt.plot(my_history.history['val_loss'])
plt.title('model mean squared error loss')
plt.ylabel('mean squared error loss')
plt.xlabel('epoch')
plt.legend(['training set', 'validation set'], loc='upper right')
plt.show()


