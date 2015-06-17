/* 
 * File:   point_types.h
 * Author: Katie
 *
 * Created on June 12, 2015, 03:30 PM
 */

#ifndef POINT_TYPES_H
#define POINT_TYPES_H

//=================================
// included dependencies
#include <pcl/point_types.h>
#include <iostream>
#include <string>
#include <pcl/io/pcd_io.h>
#include "boost/program_options.hpp"
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/visualization/pcl_visualizer.h>

//Macros
#define VISUALIZATION_SIMPLE 0
#define VISUALIZATION_SHAPE 1
#define PLANE 0
#define CYLINDER 1
#define SPHERE 2

// Using declarations
using namespace std;
using boost::program_options::options_description;

//Point Types
typedef pcl::PointXYZ PointT;

// Globals
extern bool VERBOSE;
extern bool VISUAL;
extern int TYPE_VISUALIZATION;
extern bool OUTPUT;
//Change to segment particular models
extern bool search_plane;
extern bool search_sphere;
extern bool search_cyl;

//Constants
const int SHAPES= 3;
extern const string shape_names[3]; //Maps from Macro NUM to shape name (as text).


//=======================================================================================
// Segmentation helper functions

/**
 * Checks whether any of the models defined in PCL can still be parametrized on the dataset.
 * If not, the segmentation process will end.
 **/
bool search_valid();

/* Flips the values of the flag that indicate a green light in regards to searching
 * for a particular model on the point cloud data. Useful for initializing flags
 * based on user preferences at start. */
bool search_flip();

/**
 * Update a cloud given a list of inliers. Dumps inliers to segment cloud.
 * APPENDS segment to LIST of clouds.
 */
void updateCloud(pcl::PointCloud<PointT>::Ptr cloud, vector<int> inliers,
		vector<pcl::PointCloud<PointT>::Ptr>* clouds, pcl::PointCloud<PointT>::Ptr updated,
		pcl::PointCloud<pcl::Normal>::Ptr working_normals);

/**
 * Overloaded function. Update a cloud given a list of inliers. Dumps inliers to segment cloud.
 */
void updateCloud(pcl::PointCloud<PointT>::Ptr cloud, vector<int> inliers,
		pcl::PointCloud<PointT>::Ptr * segment_cloud, pcl::PointCloud<PointT>::Ptr updated,
		pcl::PointCloud<pcl::Normal>::Ptr working_normals);

// Compute the best ransac fit for a single cylinder
vector<int> ransacCylinder(pcl::PointCloud<PointT>::Ptr updated, pcl::PointCloud<pcl::Normal>::Ptr cloud_normals,
		pcl::ModelCoefficients::Ptr coefficients,
		float threshold=0.02, float max_radius=0.1, float normal_weight=0.1, int max_iterations=5000);

// Compute the best ransac fit for a single plane
vector<int> ransacPlane(pcl::PointCloud<PointT>::Ptr updated, pcl::PointCloud<pcl::Normal>::Ptr cloud_normals,
		pcl::ModelCoefficients::Ptr coefficients,
		float threshold=0.03, float normal_weight=0.1, int max_iterations=500);

// Compute the best ransac fit for a single sphere
vector<int> ransacSphere(pcl::PointCloud<PointT>::Ptr updated, pcl::PointCloud<pcl::Normal>::Ptr cloud_normals,
		pcl::ModelCoefficients::Ptr coefficients,
		float threshold=0.01, float normal_weight=0.1, int max_iterations=100);

//=======================================================================================
// Visualization functions
boost::shared_ptr<pcl::visualization::PCLVisualizer> rgbVis (vector<pcl::PointCloud<PointT>::Ptr> clouds);

boost::shared_ptr<pcl::visualization::PCLVisualizer> viewportsVis (
		vector<pcl::PointCloud<PointT>::Ptr> clouds,
		std::vector<std::vector<float> > models);

// Visualize the point clouds in different colors
void visualizeCloud(vector<pcl::PointCloud<PointT>::Ptr> clouds, std::vector<std::vector<float> > models);


#endif	/* POINT_TYPES_H */

