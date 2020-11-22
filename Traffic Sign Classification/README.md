# **Traffic Sign Recognition** 

## Writeup - Project Report

---

**Build a Traffic Sign Recognition Project**

The goals / steps of this project are the following:
* Load the data set (see below for links to the project data set)
* Explore, summarize and visualize the data set
* Design, train and test a model architecture
* Use the model to make predictions on new images
* Analyze the softmax probabilities of the new images
* Summarize the results with a written report


[//]: # (Image References)

[image1]: ./examples/vis.png "Visualization"
[image2]: ./examples/standard.png "Standarizing"
[image3]: ./examples/top5.png "Top5"

[image4]: ./my_web_examples_signs/1.jpg "Traffic Sign 1"
[image5]: ./my_web_examples_signs/9.jpg  "Traffic Sign 2"
[image6]: ./my_web_examples_signs/27.jpg  "Traffic Sign 3"
[image7]: ./my_web_examples_signs/39.jpg  "Traffic Sign 4"
[image8]: ./my_web_examples_signs/18.jpg "Traffic Sign 5"


#### 1. Basic summary of the data set.

I used the pandas library to calculate summary statistics of the traffic
signs data set:

* The size of training set is 34799
* The size of the validation set is 4410
* The size of test set is 12630
* The shape of a traffic sign image is 32x32x3
* The number of unique classes/labels in the data set is 43

#### 2. Include an exploratory visualization of the dataset.

Here is an exploratory visualization of the Training data set. It is a bar chart showing how the data are distributed in the 43 unique classes.

![alt text][image1]

#### 2. Preprocess of the image data

The image data should be normalized so that the data has mean zero and equal variance. For image data, (pixel - 128)/ 128 is suggested as a quick way to approximately normalize the data and this is what i followed.

Also in the examples from the web I had to standarize them first in order to be 32x32x3 before feeding them into the neural network.

![alt text][image2]


#### 3. Final model architecture

My final model consisted of the following layers:

| Layer         		|     Description	        					| 
|:---------------------:|:---------------------------------------------:| 
| Input         		| 32x32x3 RGB image   							| 
| Convolution 3x3     	| 1x1 stride, valid padding, outputs 28x28x6 	|
| RELU					|												|
| Max pooling	      	| 1x1 stride,  outputs 14x14x6 				|
| Convolution 3x3	    | 1x1 stride, valid padding, outputs 10x10x6  									|
| RELU                    |                                                |
| Max pooling              | 1x1 stride,  outputs 5x5x6                 |
| DROPOUT                    | rate=0.8                                              |
| Fully connected		| output = 84       									|
| RELU                    |                                                |
| DROPOUT                    | rate=0.8                                              |
| Fully connected        | output = 43                                           |


My model was basically inspired by  LeNet network model presented in the classroom.  It consists of 2 convolutional layers with filters with dimensions 3x3. The convolutional layers are followed by 2 fully connected layers with outputs of 84 and 43 respectively.

The model includes RELU layers to introduce nonlinearity. The model contains dropout layers between convolutional and fully connected layers in order to reduce overfitting.
 
Through testing and checking the accuracy of the validation test I managed to tune and set the hyperparameters, dropout rate and number of epochs to the final value.

#### 4. Model Training and results

To train the model, I used an the LeNet architecture that was taught in the lessons but in order to increase the accuracy and avoid overfitting I included additions dropouts layers.

The model used:
1) an adam optimizer with learning rate of 0.001.
2) epoch 30 which after testing proved to be effiecient
3) Batch size of 64 which after testing proved to be effiecient
4) mu = 0 unchnaged from LeNet
5) sigma = 0.1 modified to increase accuracy

My final model results were:
* training set accuracy of 0.985
* validation set accuracy of 0.934
* test set accuracy of 0.917

#### 5. My Web-Images

Here are five German traffic signs that I found on the web:

![alt text][image4] ![alt text][image5] ![alt text][image6] 
![alt text][image7] ![alt text][image8]

Here are the results of the prediction:

| Image			        |     Prediction	        					| 
|:---------------------:|:---------------------------------------------:| 
| No passing     		| No passing									|
| Keep left    			| Keep left									|
| Pedestrians					|  Right-of-way at the next intersection										|
| General caution	      		| Turn left ahead					 				|
| Speed limit (30km/h)			| Slippery Road      							|

The model did surpringly well with 60% accuracy compared to 90% in test accuracy.

Since most images need to be resized infromation will be lost. Another factor that may affect accuarcy of the new images is of course the resolution which changes when we resize the pictures. All of these were confirmed in the above 5 examples as well as other images I tried for further testing in my spare time.

#### 6.  Top5 Softmax

![alt text][image3]

