# **Finding Lane Lines on the Road** 

## Writeup Template - Project Reflection

### This is my writeup file as part of my project submission

---

**Finding Lane Lines on the Road**

The goals / steps of this project are the following:
* Make a pipeline that finds lane lines on the road
* Reflect on your work in a written report


[//]: # (Image References)

[image1]: ./examples/grayscale.jpg "Grayscale"

[image2]: ./tests_images/output.jpg "Final Output"

---

### Reflection

### 1. Describe your pipeline. As part of the description, explain how you modified the draw_lines() function.

My pipeline consisted of 5 steps. First, I converted the images to grayscale, then I did the Guassian filtering using kernel size of 5.

Next I did Kanny filtering to find edges in the image and masked the image.

Last I did Hough Transform in the image using threshold values based on measurements and experimentation. The way i implemented the draw_lines function for the transform is the following:

1. Find the slope and based on the sign decide is its a left or right point
2. Using the polyfit function I find the average ideal polynomial to the coordinates (the definition of the function although was not covered during the lesson, based on internet reserach is advisable to use in that kind of tasks. It gives back the Ax+B coordinates mapping to y=mx+b. Defining yman and ymin we can use the polyfit coordinates to find the ideal xmax and xmin).
3. Correct the coordinates for the picture
4. Draw the lines

Here is an example image:

![image1]
![image2]

### 2. Identify potential shortcomings with your current pipeline


One potential shortcoming would be that the threshold levels are based on measurements for the video examples. Making the algorithm not universal.


### 3. Suggest possible improvements to your pipeline

Better selection in the threshold levels is a possible improvement. Also the concept of intersection in the Hough plane which is something taught but not used in my implementation.
