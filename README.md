# BTBR
All scripts and data for our BTBR study.
The original image files can be requested. Here we provided the images already imported into MATLAB.
For the images generated for Figure 1 use DetailedMOrphImFinal.m. Just run each section. 
The image enhancing algorithm is basically a histogram equalizer that runs over tiles over the images. The images are also thresholded. 
THere is a final pass of histogram equalization and a Gaussian convoluation. The tiling is done to enhance images locally and to avoid edge articats.

For the data presented in Figures 2 and 3 run FinalCleanAnalysis_ForPublication.m
This file is self contained and can be easily run. 

All inquiries to Prof. Fidel Santamaria fidel.santamaria@utsa.edu
