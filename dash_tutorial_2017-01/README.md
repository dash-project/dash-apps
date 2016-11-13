# DASH tutorial

This is a set of assignments with solutions for an upcoming DASH tutorial. It uses a very large astronomical image that is loaded into a block-distributed 2D DASH matrix.

## Input image

I recommend the 31813x19425 928 MB image 'heic1620a.tif' from http://www.spacetelescope.org/images/heic1620a/. Any other astronomical image in TIFF format will do if it it stored in a 'stripped' manner.

## Requirements:

The astro examples need:
* libtiff5-dev
* libsdl1.2-dev (not libsdl2) This is not required but optionally, build with -DWITHSDL in this case.

## Assignments

The assignments are given codes with some DASH data container allocations or operations on them missing. Participants need to fill in few lines of code in order to:

### 1. A hello world intro

### 2. Load image, optionally show a part of the image
* Declare the 2D block-distributed DASH matrix.
* Copy the images line by line from the libtiff buffer to the distributed matrix using a DASH algorithm.

### 3. Compute a brightness histogram of the pixel values
* First, compute independent local histograms in each local block in parallel.
* Then, sum the local histograms using a DASH algorithm to produce the global result.
* Finally, looking at the brightness histogram define the limit between what is considered as bright pixels and dark pixels.

### 4. Choose a marker color that is not appearing in the image yet
* Do a quick parallel scan over all pixels to confirm, that the marker color is not found.

### 5. Count the number of bright objects in the image.
* For every pixel, run a given 'checkobject' function. It will report, if it was a bright pixel and, if so, mark all adjacent pixels with the marker color, so that they won't be reported again.
* Count in parallel for every block and combine the result into the global result.
* (For simplixity, we don't correct for bright objects at the border between distributed block, which may be countet multiple time. Thie correction is left as a further exercise ...)
