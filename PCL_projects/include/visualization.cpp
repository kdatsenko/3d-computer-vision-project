#include <iostream>
#include <pcl/ModelCoefficients.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/common/common.h>
#include <pcl/common/transforms.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <point_types.h>

/*
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/filters/extract_indices.h>
*/


//-----------------------------------------------------------------------------------------------------------------

// Cloud Visualization Methods ---------------------------------------------------------------------------------------------
// Open 3D viewer to visualize the segmented pointcloud
boost::shared_ptr<pcl::visualization::PCLVisualizer> rgbVis (vector<pcl::PointCloud<PointT>::Ptr> clouds) {
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
	viewer->setBackgroundColor (0, 0, 0);
	for (int i=0; i<(int)clouds.size(); i++){
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_xyzrgb (new pcl::PointCloud<pcl::PointXYZRGB>);
		copyPointCloud(*clouds[i], *cloud_xyzrgb);
		char label[50];
		sprintf(label, "cloud %d%s", i, "\0");
		pcl::visualization::PointCloudColorHandlerRandom<pcl::PointXYZRGB> handler (cloud_xyzrgb);
		viewer->addPointCloud (cloud_xyzrgb, handler, label);
		viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, label);
	}
	viewer->addCoordinateSystem	(1.0);
	viewer->initCameraParameters ();

	return (viewer);
}

boost::shared_ptr<pcl::visualization::PCLVisualizer> viewportsVis (
		vector<pcl::PointCloud<PointT>::Ptr> clouds,
		std::vector<std::vector<float> > models)
				{
	// --------------------------------------------------------
	// -----Open 3D viewer and add point cloud and normals-----
	// --------------------------------------------------------
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
	viewer->initCameraParameters ();

	/* Initialize the two viewports */
	int v1(0);
	viewer->createViewPort(0.0, 0.0, 0.5, 1.0, v1);
	viewer->setBackgroundColor (0, 0, 0, v1);
	viewer->addText("Radius: 0.01", 10, 10, "v1 text", v1);
	int v2(0);
	viewer->createViewPort(0.5, 0.0, 1.0, 1.0, v2);
	viewer->setBackgroundColor (0.3, 0.3, 0.3, v2);
	viewer->addText("Radius: 0.1", 10, 10, "v2 text", v2);

	for (int i=0; i<(int)clouds.size(); i++){
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_xyzrgb (new pcl::PointCloud<pcl::PointXYZRGB>);
		copyPointCloud(*clouds[i], *cloud_xyzrgb);
		char label1[50];
		char label2[50];
		sprintf(label1, "cloud v1 %d%s", i, "\0");
		sprintf(label2, "cloud v2 %d%s", i, "\0");

		double c[3] = {(rand() % (255 + 1)), (rand() % (255 + 1)), (rand() % (255 + 1))};
		pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZRGB> handler(cloud_xyzrgb, c[0], c[1], c[2]);
		viewer->addPointCloud<pcl::PointXYZRGB> (cloud_xyzrgb, handler, label1, v1);

		//point cloud is monochrome in second viewport
		pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZRGB> single_color(cloud_xyzrgb, 0, 255, 0);
		viewer->addPointCloud<pcl::PointXYZRGB> (cloud_xyzrgb, single_color, label2, v2);
		//borrow colours from matching point clouds to render shapes

		pcl::ModelCoefficients coeffs;
		for (int j=1; j < models[i].size(); j++)
			coeffs.values.push_back (models[i][j]);

		Eigen::Vector4f centroid;
		pcl::compute3DCentroid(*cloud_xyzrgb, centroid);

		//cout << "centroid: " << centroid[0] << " " << centroid[1] << " " << centroid[2] << "\n";
		//centroid.
		char shape_label[50];
		switch ((int)models[i][0]) {
		case PLANE:
			sprintf(shape_label, "plane_%d%s", i, "\0");
			viewer->addPlane(coeffs, (double)centroid[0], (double)centroid[1], (double)centroid[2], shape_label, v2); //
			break;
		case CYLINDER:
		  sprintf(shape_label, "  cylinder_%d%s", i, "\0");
		  viewer->addCylinder(coeffs, shape_label, v2);
		  break;
		case SPHERE:
		  sprintf(shape_label, "sphere_%d%s", i, "\0");
		  viewer->addSphere(coeffs, shape_label, v2);
		  break;
		default:
		  cout << "value of x unknown";
		}
		//A value of 1.0 is equivalent to 255, a value of 0.0 to 0.
		double r = c[0] / 255.0;
		double g = c[1] / 255.0;
		double b = c[2] / 255.0;

		//cout << "colour: " << c[0] << " " << c[1] << " " << c[2] << "\n";
		viewer->setShapeRenderingProperties (pcl::visualization::PCL_VISUALIZER_COLOR, r, g, b /*R,G,B*/, shape_label, v2);
		viewer->setShapeRenderingProperties (pcl::visualization::PCL_VISUALIZER_OPACITY, 0.9, shape_label, v2);
		viewer->setShapeRenderingProperties (pcl::visualization::PCL_VISUALIZER_REPRESENTATION,
				pcl::visualization::PCL_VISUALIZER_REPRESENTATION_SURFACE, shape_label, v2);
		viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, label1);
		viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 0.1, label2);
	}

	viewer->addCoordinateSystem (1.0);
	return (viewer);
}

// Visualize the point clouds in different colors
void visualizeCloud(vector<pcl::PointCloud<PointT>::Ptr> clouds, std::vector<std::vector<float> > models){
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer;
	switch (TYPE_VISUALIZATION) {
		case VISUALIZATION_SIMPLE:
			viewer= rgbVis(clouds);
			break;
		case VISUALIZATION_SHAPE:
			viewer= viewportsVis(clouds, models);
			break;
		default:
			viewer= rgbVis(clouds);
	}

	while (!viewer->wasStopped ())
	{
		viewer->spinOnce (100);
		boost::this_thread::sleep (boost::posix_time::microseconds (100000));
	}
}
