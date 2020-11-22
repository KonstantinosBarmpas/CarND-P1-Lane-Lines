# **Behavioral Cloning** 

## Writeup - Project Report

---

**Behavioral Cloning Project**

The goals / steps of this project are the following:
* Use the simulator to collect data of good driving behavior
* Build, a convolution neural network in Keras that predicts steering angles from images
* Train and validate the model with a training and validation set
* Test that the model successfully drives around track one without leaving the road
* Summarize the results with a written report


[//]: # (Image References)

[image1]: ./examples/nvdia.png "Model Visualization"
[image2]: ./examples/center.png "Center"
[image3]: ./examples/recovery1.png "Recovery Image"
[image4]: ./examples/recovery1.png "Recovery Image"
[image5]: ./examples/loss.png "Loss"

---

Watch the final video here:

https://www.barmpas.com/self-driving-cars

### Files Submitted & Code Quality

#### 1. Submission includes all required files and can be used to run the simulator in autonomous mode

My project includes the following files:
* model.py containing the script to create and train the model
* drive.py for driving the car in autonomous mode
* model.h5 containing a trained convolution neural network 
* writeup_report.md or writeup_report.pdf summarizing the results

#### 2. Submission includes functional code
Using the Udacity provided simulator and my drive.py file, the car can be driven autonomously around the track by executing 
```sh
python drive.py model.h5
```

#### 3. Submission code is usable and readable

The model.py file contains the code for training and saving the convolution neural network. The file shows the pipeline I used for training and validating the model, and it contains comments to explain how the code works.

### Model Architecture and Training Strategy

#### 1. An appropriate model architecture has been employed

My model was basically inspired by NVDIA's network model presented in the classroom.  It consists of 5 convolutional layers. The first three have filters with dimensions 5x5 while the last two have 3x3. All convolutional layers depths are between 24-64. (model.py lines 49-85). The convolutional layers are followed by 3 fully connected layers with outputs of 100,50,10 respectively.

The model includes RELU layers to introduce nonlinearity (code line 58), and the data is normalized in the model using a Keras lambda layer (code line 54).

#### 2. Attempts to reduce overfitting in the model

The model contains dropout layers in order to reduce overfitting (model.py lines 51).

The model was trained and validated on different data sets to ensure that the model was not overfitting (code line 88). The model was tested by running it through the simulator and ensuring that the vehicle could stay on the track.

#### 3. Model parameter tuning

The model used an adam optimizer with learning rate of 0.0001 (model.py lines 97).

#### 4. Appropriate training data

Appropriate data were chosen from different angles of the route and also by collecting data doing the route clockwise.

For details about how I created the training data, see the next section. 

### Model Architecture and Training Strategy

#### 1. Solution Design Approach

My overal solution was to design a model based on NVDIA's proposed model. I also added dropouts to avoid overfitting.

First I collected the data and trained the model. After the model was tested I tuned the correction and dropout parameters as well as collecting more data. The model was retrained and the same procedure was repeated until I reached the final model.

#### 2. Final Model Architecture

My model was basically inspired by NVDIA's network model presented in the classroom.  It consists of 5 convolutional layers. The first three have filters with dimensions 5x5 while the last two have 3x3. All convolutional layers depths are between 24-64. (model.py lines 49-85). The convolutional layers are followed by 3 fully connected layers with outputs of 100,50,10 respectively.

The model includes RELU layers to introduce nonlinearity (code line 58), and the data is normalized in the model using a Keras lambda layer (code line 54).

![alt text][image1]

#### 3. Creation of the Training Set & Training Process

To capture good driving behavior, I first recorded two laps on track one using center lane driving. Here is an example image of center lane driving:

![alt text][image2]

I then recorded the vehicle recovering from the left side and right sides of the road back to center so that the vehicle would learn to to recover.

![alt text][image3]
![alt text][image4]

Then I repeated this process on track two in order to get more data points.

After the collection process, I had X number of data points. I then preprocessed this data by normalizing them and peform cropping (model.py line 56)

Finally I split the data into training and validation set and trained the network. Visualization of the training process mean square loss:

![alt text][image5]
