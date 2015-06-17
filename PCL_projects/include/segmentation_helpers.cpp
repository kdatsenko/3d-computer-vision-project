#include <pcl/ModelCoefficients.h>
#include <pcl/point_types.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/features/normal_3d.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/sample_consensus/ransac.h>
#include <pcl/sample_consensus/sac_model_sphere.h>
#include <point_types.h>







/*
#include <iostream>
#include <pcl/sample_consensus/sac_model_plane.h>
#include <pcl/sample_consensus/sac_model_sphere.h>

*/

/* Maps the names of the shapes to their numerical identifier (as the index). */
const string shape_names[3] = {"plane", "cylinder", "sphere"};

/* Checks whether one of the corresponding models can still be parametrized on the
 * remaining data points. */
bool search_valid(){
	return (search_plane | search_sphere | search_cyl);
}

/* Flips the values of the flag that indicate a green light in regards to searching
 * for a particular model on the point cloud data. */
bool search_flip(){
  search_plane = !search_plane;
  search_cyl = !search_cyl;
  search_sphere = !search_sphere;
}

/**
 * Update a cloud given a list of inliers. Dumps inliers to segment cloud. APPENDS
 * segment to LIST of clouds.
 */
void updateCloud(pcl::PointCloud<PointT>::Ptr cloud, vector<int> inliers,
		vector<pcl::PointCloud<PointT>::Ptr>* clouds, pcl::PointCloud<PointT>::Ptr updated,
		pcl::PointCloud<pcl::Normal>::Ptr working_normals){

	// copies all inliers of the model computed to another PointCloud
	pcl::PointCloud<PointT>::Ptr final (new pcl::PointCloud<PointT>);
	pcl::copyPointCloud<PointT>(*cloud, inliers, *final);
	(*clouds).push_back(final);

	// Remove inliers from original cloud
	pcl::PointIndices::Ptr fInliers (new pcl::PointIndices);
	fInliers->indices= inliers;
	pcl::ExtractIndices<PointT> extract;
	extract.setInputCloud(cloud);
	extract.setIndices(fInliers);
	extract.setNegative(true); // Removes part_of_cloud from full cloud  and keep the rest
	extract.filter(*updated);

	//Remove inliers from the original normals
	pcl::ExtractIndices<pcl::Normal> extractNormals;
	extractNormals.setInputCloud(working_normals);
	extractNormals.setIndices(fInliers);
	extractNormals.setNegative(true); // Removes part_of_cloud from full cloud  and keep the rest
	extractNormals.filter(*working_normals);
	assert(updated->size() == working_normals->size() );

}

/**
 * ---------------------------OVERLOADED---------------------------------
 * Update a cloud given a list of inliers. Dumps inliers to segment cloud.
 */
void updateCloud(pcl::PointCloud<PointT>::Ptr cloud, vector<int> inliers,
		pcl::PointCloud<PointT>::Ptr * segment_cloud, pcl::PointCloud<PointT>::Ptr updated,
		pcl::PointCloud<pcl::Normal>::Ptr working_normals){

	// copies all inliers of the model computed to another PointCloud
	pcl::PointCloud<PointT>::Ptr final (new pcl::PointCloud<PointT>);
	pcl::copyPointCloud<PointT>(*cloud, inliers, *final);
	*segment_cloud = final;

	// Remove inliers from original cloud
	pcl::PointIndices::Ptr fInliers (new pcl::PointIndices);
	fInliers->indices= inliers;
	pcl::ExtractIndices<PointT> extract;
	extract.setInputCloud(cloud);
	extract.setIndices(fInliers);
	extract.setNegative(true); // Removes part_of_cloud from full cloud  and keep the rest
	extract.filter(*updated);

	//Remove inliers from the original normals
	pcl::ExtractIndices<pcl::Normal> extractNormals;
	extractNormals.setInputCloud(working_normals);
	extractNormals.setIndices(fInliers);
	extractNormals.setNegative(true); // Removes part_of_cloud from full cloud  and keep the rest
	extractNormals.filter(*working_normals);
	assert(updated->size() == working_normals->size() );
}


// RANSAC Methods ---------------------------------------------------------------------------------------------

/**
 * Functions for the segmentation of various Models using RANSAC
 **/

// Compute the best ransac fit for a single cylinder
vector<int> ransacCylinder(pcl::PointCloud<PointT>::Ptr updated, pcl::PointCloud<pcl::Normal>::Ptr cloud_normals,
		pcl::ModelCoefficients::Ptr coefficients,
		float threshold /*=0.02*/, float max_radius /*=0.1*/, float normal_weight /*=0.1*/, int max_iterations /*=5000*/){
	pcl::SACSegmentationFromNormals<PointT, pcl::Normal> seg;
	pcl::PointIndices::Ptr inliers_cyl (new pcl::PointIndices);
	//pcl::ModelCoefficients::Ptr coefficients_cylinder (new pcl::ModelCoefficients);

	// Create the segmentation object for cylinder segmentation and set all the parameters
	seg.setOptimizeCoefficients (true);
	seg.setModelType (pcl::SACMODEL_CYLINDER);
	seg.setMethodType (pcl::SAC_RANSAC);
	seg.setNormalDistanceWeight (normal_weight);
	seg.setMaxIterations (5000);
	//seg->setDistanceThreshold (0.02);
	seg.setRadiusLimits (0, max_radius);
	seg.setDistanceThreshold (threshold);//THRESHOLD*CYL_THRESH_MOD);
	seg.setInputCloud (updated);
	seg.setInputNormals (cloud_normals);

	// Obtain the cylinder inliers and coefficients
	seg.segment (*inliers_cyl, *coefficients);
	if (inliers_cyl->indices.size() == 0){
		if (VERBOSE) PCL_ERROR ("Could not estimate a cylindrical model for the given dataset.\n");
	}
	return inliers_cyl->indices;
}

// Compute the best ransac fit for a single plane
vector<int> ransacPlane(pcl::PointCloud<PointT>::Ptr updated, pcl::PointCloud<pcl::Normal>::Ptr cloud_normals,
		pcl::ModelCoefficients::Ptr coefficients,
		float threshold /*=0.03*/, float normal_weight /*=0.1*/, int max_iterations /*=500*/){

	pcl::SACSegmentationFromNormals<PointT, pcl::Normal> seg;
	pcl::PointIndices::Ptr inliers (new pcl::PointIndices);
	// Create the segmentation object for the planar model and set all the parameters
	seg.setOptimizeCoefficients(true);
	seg.setModelType (pcl::SACMODEL_NORMAL_PLANE);
	seg.setNormalDistanceWeight (normal_weight);
	seg.setMethodType (pcl::SAC_RANSAC);
	seg.setMaxIterations (max_iterations);
	//maximum distance points may lie from the model, to be an inlier
	seg.setDistanceThreshold (threshold); //0.03
	seg.setInputCloud (updated);
	seg.setInputNormals (cloud_normals);

	// Obtain the plane inliers and coefficients
	seg.segment (*inliers, *coefficients);
	if (inliers->indices.size() == 0){
		if (VERBOSE) PCL_ERROR ("Could not estimate a planar model for the given dataset.\n");
	}
	return inliers->indices;
}

// Compute the best ransac fit for a single sphere
vector<int> ransacSphere(pcl::PointCloud<PointT>::Ptr updated, pcl::PointCloud<pcl::Normal>::Ptr cloud_normals,
		pcl::ModelCoefficients::Ptr coefficients,
		float threshold /*=0.01*/, float normal_weight /*=0.1*/, int max_iterations /*=100*/){

	/*pcl::SACSegmentationFromNormals<PointT, pcl::Normal> seg;
	pcl::PointIndices::Ptr inliers (new pcl::PointIndices);

	seg.setOptimizeCoefficients(true);
	seg.setModelType (pcl::SACMODEL_NORMAL_SPHERE);
	seg.setNormalDistanceWeight (normal_weight);
	seg.setMethodType (pcl::SAC_RANSAC);
	seg.setMaxIterations (max_iterations);
	//maximum distance points may lie from the model, to be an inlier
	seg.setDistanceThreshold (threshold); //0.03
	seg.setInputCloud (updated);
	seg.setInputNormals (cloud_normals);
	// Obtain the plane inliers and coefficients
	seg.segment (*inliers, *coefficients);
	if (inliers->indices.size() == 0){
		if (VERBOSE) PCL_ERROR ("Could not estimate a spherical model for the given dataset.");
	}
	*/

	vector<int> inliers;
	pcl::SampleConsensusModelSphere<PointT>::Ptr model_sph (new pcl::SampleConsensusModelSphere<PointT> (updated));
	pcl::RandomSampleConsensus<PointT> ransac (model_sph);
	ransac.setDistanceThreshold (threshold);
	ransac.setMaxIterations(max_iterations);
	ransac.computeModel();
	ransac.getInliers(inliers);

	Eigen::VectorXf coeff;
	ransac.getModelCoefficients(coeff);
	for (int i = 0; i < coeff.size(); i++)
	  coefficients->values.push_back(coeff[i]);

	if (inliers.size() == 0){
	  if (VERBOSE) PCL_ERROR ("Could not estimate a spherical model for the given dataset.\n");
	}
	return inliers;
}

// End RANSAC Methods ---------------------------------------------------------------------------------------------






