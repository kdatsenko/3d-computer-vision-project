/*
 * helpers.h
 *
 *  Created on: 10 Jul 2015
 *      Author: osboxes
 */

#ifndef HELPERS_H_
#define HELPERS_H_

//=================================
// included dependencies
#include <pcl/point_types.h>
#include <iostream>
#include <pcl/io/pcd_io.h>

/* Compute Hausdorff (maximum possible) distance between two Point Clouds */
float compute (pcl::PointCloud<pcl::PointXYZ> &cloud_a, pcl::PointCloud<pcl::PointXYZ> &cloud_b);





#endif /* HELPERS_H_ */





