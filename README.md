### Introduction

This is an implementation of the dominant vanishing point detection method described in the paper:

[*Detecting Dominant Vanishing Points In Natural Scenes with Application to Composition-Sensitive Image Retrieval.*](https://faculty.ist.psu.edu/zzhou/projects/vpdetection/) Zihan Zhou, Farshid Farhat, and James Z. Wang. IEEE Transactions on Multimedia (TMM), 2017. 

### Usage
* run main.m
* The code is tested on Mac OS X El Capitan. If you want to use it on other platform, you need to compile the .mex files using the source code included.
* The current package includes the [original Berkeley contour detection algorithm](https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html), which was used in the paper. However, this algorithm is quite slow (>1 min for an image). For fast computation of the ultrametic contour map (UCM), you can replace this module with newer and faster algorithms such as [this one](https://github.com/jponttuset/mcg).
* The size of the input images is assumed to be 500 pixels by default. It is recommended to resize the images to 500 pixels before running the code.

### References
If you use this software you have to reference ALL of these papers:

1. *Detecting Dominant Vanishing Points In Natural Scenes with Application to Composition-Sensitive Image Retrieval.* Zihan Zhou, Farshid Farhat, and James Z. Wang. IEEE Transactions on Multimedia, 2017.

2. *Contour Detection and Hierarchical Image Segmentation*. 
P. Arbelaez, M. Maire, C. Fowlkes and J. Malik.
IEEE Transactions on Pattern Recognition and Machine Intelligence, Vol. 33, No. 5, pp. 898-916, May 2011.

3. *Non-iterative Approach for Fast and Accurate Vanishing Point Detection*, Jean-Philippe Tardif. ICCV, 2009.


### Copyright and License

Copyright 2017 Zihan Zhou

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

Thanks
