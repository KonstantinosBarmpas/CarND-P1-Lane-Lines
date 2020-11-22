# ***Vehicles Detection**

## Writeup Template - Project Report


**Vehicle Detection Project**

Watch the video at:

https://www.barmpas.com/self-driving-cars

The goals / steps of this project are the following:

* Perform a Histogram of Oriented Gradients (HOG) feature extraction on a labeled training set of images and train a classifier Linear SVM classifier
* Optionally, you can also apply a color transform and append binned color features, as well as histograms of color, to your HOG feature vector.
* Note: for those first two steps don't forget to normalize your features and randomize a selection for training and testing.
* Implement a sliding-window technique and use your trained classifier to search for vehicles in images.
* Run your pipeline on a video stream (start with the test_video.mp4 and later implement on full project_video.mp4) and create a heat map of recurring detections frame by frame to reject outliers and follow detected vehicles.
* Estimate a bounding box for vehicles detected.

[//]: # (Image References)
[image1]: ./examples/1.png
[image2]: ./examples/2.png
[image3]: ./examples/3.png
[image4]: ./examples/4.png

### Histogram of Oriented Gradients (HOG)

#### 1. Explain how (and identify where in your code) you extracted HOG features from the training images.

The code for this step is contained in the

I started by reading in all the `vehicle` and `non-vehicle` images.  Here is an example of one of each of the `vehicle` and `non-vehicle` classes:

![alt text][image1]

I then explored different color spaces and different `skimage.hog()` parameters (`orientations`, `pixels_per_cell`, and `cells_per_block`).  I grabbed random images from each of the two classes and displayed them to get a feel for what the `skimage.hog()` output looks like:

![alt text][image2]

#### 2. Final choice of HOG parameters.

I tried various combinations of parameters and and the final values of the parameters are :

color_space = 'RGB' # RGB, HSV, LUV, HLS, YUV, YCrCb
orient = 5  # HOG orientations
pix_per_cell = 9 # HOG pixels per cell
cell_per_block = 2 # HOG cells per block
hog_channel = 'ALL' # 0, 1, 2,"ALL"
spatial_size = (16, 16) # Spatial binning dimensions
hist_bins = 16    # Number of histogram bins
spatial_feat = True # Spatial features True or False
hist_feat = True # Histogram features True or False
hog_feat = True # HOG features True or False

The main criterion behind this choice was that this combination has a high accuracy, eliminates the error to identify cars and can be used perfectly with windowing and heatmap techniques to give me the final beuatiful result.

#### 3. Describe how (and identify where in your code) you trained a classifier using your selected HOG features (and color features if you used them).

I trained a linear SVM using the exctract_features function. This combines three methods introduced in the classroom. First the HOG feutures described above and two other functions color_hist and bin_spatial to give additional feutures to the classifier.

Then I created a StandardScaler() to fit the data and normalized them. After that I splitted the data to a training set and a test set.

Finally, I used a SVC classifier and trained it in classfier() function using the training set and tested it using the test set.

### Sliding Window Search

#### 1. Describe how (and identify where in your code) you implemented a sliding window search.  How did you decide what scales to search and how much to overlap windows?

y_start_stop = [350, 656] # Min and max in y to search in the sliding window
scale = 2

These were the parameters I used for the sliding window to draw boxes. Everytime I found one window I added the heat to the heatmap. Finally, I used that heatmap and threshold  of 2 to avoid false positives in the apply_threshold and draw_labeled_bboxes functions.

I decided to search random window positions at random scales all over the image and came up with this (ok just kidding I didn't actually ;):

![alt text][image3]

#### 2. Test images

Here are the test images:

![alt text][image4]
---

### Combine the pipelines

I also combined the two pipelines from project 4 and project 5 to have a final Computer Vision project

### Discussion

#### 1. Briefly discuss any problems / issues you faced in your implementation of this project.  Where will your pipeline likely fail?  What could you do to make it more robust?

I believe that the pipeline works pretty good. There are cases of false positives but after investigation I noticed that happens in the dark side of the road. Possible ways to fix it would be more data, better choice of classifier, better choice of classifier's parameters or better choice of features to exctract. Thats a common computer vision problem and an area where deep learning and neural networks have been proved to amazingly smart.

