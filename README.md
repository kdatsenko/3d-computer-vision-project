# 3d-computer-vision-project

3-D computer vision project developed during an NSERC Undergraduate Student Research Internship at University of Toronto.  

Designed and implemented algorithms for segmenting (partitioning) RGBD images into regions/surfaces. Generated multi-scale region/surface segmentation proposals from which viewpoint invariant features can be extracted, match to database of object models. (C++, PCL).  

Contains a region-growing/merging procedure which does an initial growth of regions using PCL libraries, and then computes a neighbour map of regions and creates a priority-queue of region pairs (r1, r2) based on merge cost. Implemented an application that preserves the history of merges, and contains a visualizer that can be used to go forward and backward through the merging procedure. 

Presentation Slides: https://drive.google.com/file/d/0B2Alcaq6E4I3dXlyWVJtTmRVSU0/
