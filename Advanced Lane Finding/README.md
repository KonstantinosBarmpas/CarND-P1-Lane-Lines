# ***Advanced Lane Finding Project**

## Writeup Template - Project Report

---

Watch the video:

https://www.barmpas.com/self-driving-cars

**Advanced Lane Finding Project**

The goals / steps of this project are the following:

* Compute the camera calibration matrix and distortion coefficients given a set of chessboard images.
* Apply a distortion correction to raw images.
* Use color transforms, gradients, etc., to create a thresholded binary image.
* Apply a perspective transform to rectify binary image ("birds-eye view").
* Detect lane pixels and fit to find the lane boundary.
* Determine the curvature of the lane and vehicle position with respect to center.
* Warp the detected lane boundaries back onto the original image.
* Output visual display of the lane boundaries and numerical estimation of lane curvature and vehicle position.

[//]: # (Image References)

[image1]: ./examples/undistort_output.png "Undistorted"
[image2]: ./examples/undi_image.png "Road Transformed"
[image3]: ./examples/gradient.png "Binary Example"
[image4]: ./examples/trans.jpg "Warp Example"
[image5]: ./examples/color_fit_lines.jpg "Fit Visual"
[image6]: ./examples/draw.png "Output"


### Camera Calibration

#### 1. Calibration image.

The code for this step is contained in the first code cell of the IPython notebook located in "P4.ipynb" (calibrate_camera() function).

I start by preparing "object points", which will be the (x, y, z) coordinates of the chessboard corners in the world. Here I am assuming the chessboard is fixed on the (x, y) plane at z=0, such that the object points are the same for each calibration image.  Thus, `objp` is just a replicated array of coordinates, and `objpoints` will be appended with a copy of it every time I successfully detect all chessboard corners in a test image.  `imgpoints` will be appended with the (x, y) pixel position of each of the corners in the image plane with each successful chessboard detection.  

I then used the output `objpoints` and `imgpoints` to compute the camera calibration and distortion coefficients using the `cv2.calibrateCamera()` function.  I applied this distortion correction to the test image using the `cv2.undistort()` function and obtained this result: 

![alt text][image1]

### Pipeline (single images)

#### 1. Provide an example of a distortion-corrected image.

To demonstrate this step, I will describe how I apply the distortion correction to one of the test images like this one:
![alt text][image2]

#### 2. Describe how (and identify where in your code) you used color transforms, gradients or other methods to create a thresholded binary image.  Provide an example of a binary image result.

I used a combination of gradient thresholds to generate a binary image (in the functions: abs_sobel, mag_thresh and dir_threshold).  Here's an example of my output for this step:

![alt text][image3]

I also included a color threshold to isolate S-channel. Here's is an example:

![alt text][image7]

Finally, I combined both of them into one single file:

![alt text][image8]


#### 3. Describe how (and identify where in your code) you performed a perspective transform and provide an example of a transformed image.

The code for my perspective transform includes a function called "perspective_transform| and I chose the hardcode the source and destination points in the following manner:

```python
src = np.float32(
    [[(img_size[0] / 2) - 55, img_size[1] / 2 + 100],
    [((img_size[0] / 6) - 10), img_size[1]],
    [(img_size[0] * 5 / 6) + 60, img_size[1]],
    [(img_size[0] / 2 + 55), img_size[1] / 2 + 100]])
dst = np.float32(
    [[(img_size[0] / 4), 0],
    [(img_size[0] / 4), img_size[1]],
    [(img_size[0] * 3 / 4), img_size[1]],
    [(img_size[0] * 3 / 4), 0]])
```

This resulted in the following source and destination points:

| Source        | Destination   | 
|:-------------:|:-------------:| 
| 585, 460      | 320, 0        | 
| 203, 720      | 320, 720      |
| 1127, 720     | 960, 720      |
| 695, 460      | 960, 0        |


![alt text][image4]

#### 4. Describe how (and identify where in your code) you identified lane-line pixels and fit their positions with a polynomial?

I used windowing and sliding windowing to find lines in the picture.

Then I fit my lane lines with a 2nd order polynomial kinda like this:

![alt text][image5]

#### 5. Describe how (and identify where in your code) you calculated the radius of curvature of the lane and the position of the vehicle with respect to center.

I did this in the functions:
1) curvature(leftx, rightx, lefty, righty, left_fitx, right_fitx, left_fit, right_fit, ploty)
2) position_to_center(left_fitx, right_fitx,result)

The first method was recommended during the lesson and it basically uses a fit second order polynomial which takes the gradient.

![alt text][image5]

The latter I used the approximate width of the road and the fact that the center of the image is the center of the camera and substracted.

#### 6. Provide an example image of your result plotted back down onto the road such that the lane area is identified clearly.

I used the following functions to draw the lines and write on the image respectively:
1) draw_the_lines(binary_warped,left_fitx, ploty,Minv,image,undist)
2) write_on_image(img,curve,distance)

![alt text][image6]

---

### Discussion

#### 1. Briefly discuss any problems / issues you faced in your implementation of this project.  Where will your pipeline likely fail?  What could you do to make it more robust?

I implemented all the techniques proposed in the simple and the advanced lines classes. My model although works pretty well I believe it fails in areas with big curves. One way to better improve it would be to set better thresholds. Another approach can be the deep learning techniques weve learnt in the previous class.
