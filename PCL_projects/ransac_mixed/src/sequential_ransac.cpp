/*
 * cylinder_plane.cpp
 *
 *  Created on: 25 May 2015
 *      Author: osboxes
 */

#include <iostream>
#include <pcl/ModelCoefficients.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/features/normal_3d.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/filters/extract_indices.h>
#include <point_types.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include "boost/program_options.hpp"


//I/O
bool VERBOSE = false;
bool VISUAL = false;
int TYPE_VISUALIZATION = 0;
bool OUTPUT = false;

//Change to segment particular models
//(Plane SACSegmentation object incorporates normals)
bool search_plane = false;
bool search_sphere = false;
bool search_cyl = false;

//Parameters for SACSegmentationFromNormals object
float cyl_thresh;
int max_cyl_iterations;
float max_cyl_radius;
float cyl_normal_weight;
float plane_thresh;
float plane_normal_weight;
int max_plane_iterations;
float sphere_thresh;
float sphere_normal_weight;
int max_sphere_iterations;

//Globals
int MAX_OUTLIERS;
const float PERCENT_OUTLIERS= 0.1; //MODIFY

int
main (int argc, char** argv)
{

	ofstream output_file;
	string output_file_name;
	string output_directory;
	double threshold;
	try {

		int visual_arg;
		namespace opt = boost::program_options;
		opt::options_description desc("Options");
		desc.add_options()
								("help,h", "Print usage options.")
								("verbose,v", "Verbose mode, print out progress during computation.")
								("visual,s", opt::value<int>(&visual_arg), "ARG (0=RGB colour | 1=Viewport shape) Visual mode, display a colored segmentation of the cloud.")
								("output,o", opt::value<string>(&output_file_name), "<output_file> Save the model parameters to this file.")
								("outdir,d", opt::value<string>(&output_directory)->default_value("./"), "<./directory> Save all output to this RELATIVE (to location of execution file) dir.")
								("noplane", opt::bool_switch(&search_plane), "Exclude planes.")
								("nocylinder", opt::bool_switch(&search_cyl), "Exclude cylinders.")
								("nosphere", opt::bool_switch(&search_sphere), "Exclude spheres.")
								("cyl_thresh,a", opt::value<float>(&cyl_thresh)->default_value(0.02), "Cylinder Distance threshold - units are relative to data pcd file, usually metres.")
								("max_cyl_iterations,b", opt::value<int>(&max_cyl_iterations)->default_value(5000), "Maximum number of iterations to try to find a Cylinder model")
								("max_cyl_radius,c", opt::value<float>(&max_cyl_radius)->default_value(0.1), "Cylinder")
								("cyl_normal_weight,d", opt::value<float>(&cyl_normal_weight)->default_value(0.1), "Set the relative weight (between 0 and 1) to give to the angular distance (0 to pi/2) between point normals and the Cylinder normal.")
								("plane_thresh,l", opt::value<float>(&plane_thresh)->default_value(0.03), "Plane Distance threshold - units are relative to data pcd file, usually metres")
								("plane_normal_weight,m", opt::value<float>(&plane_normal_weight)->default_value(0.1), "Set the relative weight (between 0 and 1) to give to the angular distance (0 to pi/2) between point normals and the Plane normal.")
								("max_plane_iterations,n", opt::value<int>(&max_plane_iterations)->default_value(500), "Maximum number of iterations to try to find a Planar model")
								("sphere_thresh,x", opt::value<float>(&sphere_thresh)->default_value(0.01), "Sphere Distance threshold - units are relative to data pcd file, usually metres")
								("sphere_normal_weight,y", opt::value<float>(&sphere_normal_weight)->default_value(0.1), "Set the relative weight (between 0 and 1) to give to the angular distance (0 to pi/2) between point normals and the Sphere normal.")
								("max_sphere_iterations,z", opt::value<int>(&max_sphere_iterations)->default_value(100), "Maximum number of iterations to try to find a Spherical model");
		opt::variables_map vm;

		try
		{	// Parse arguments and, if necessary, print usage
		  namespace cls = boost::program_options::command_line_style;
			opt::store(opt::command_line_parser(argc, argv) // can throw
						.options(desc)
						.style(
						    opt::command_line_style::unix_style) //allow shorthand opt with > 1 char
								.run(), vm);
			if (argc < 2){
				std::cerr << desc << std::endl;
				exit(0);
			}

			/** --help option */
			if ( vm.count("help")  ){
				cout << "\n Sequential Multi-Object Ransac Segmentation" << endl
						<< "\nUsage: cylinder_plane <input_cloud.pcd> [-hv] [-s (0|1)] [-o <output_file>] [-t NUM]" << endl
						<< "Purpose: Uses RANSAC algorithm to segment a PCL pointcloud (.pcd) file into "
						"multiple objects. Segments shapes one-by-one in a sequence."
						" Outputs the segmented .pcd plane cloud files." << endl
						<< desc << endl;
				exit(0);
			}
			opt::notify(vm); // throws on error, so do after help in case

			// if there are any problems
		} catch(boost::program_options::required_option& e)
		{
			std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
			std::cerr << desc << std::endl;
			return (1);
		}
		catch(boost::program_options::error& e)
		{
			std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
			std::cerr << desc << std::endl;
			return (1);
		}

		if (vm.count("visual")){ // Check visual flag
			VISUAL= true;
			TYPE_VISUALIZATION = (visual_arg == 0) ? VISUALIZATION_SIMPLE : VISUALIZATION_SHAPE;
		} if (vm.count("verbose")){ // Check verbose flag
			VERBOSE= true;
		} if (vm.count("output")){ // Check output flag
			string p = output_directory + "/" + output_file_name;
			const char * path = p.c_str();
      boost::filesystem::path dir(output_directory);

      if(!boost::filesystem::exists(dir)){
        if (boost::filesystem::create_directory(dir))
          std::cout << "Directory " << output_directory << " was created." << std::endl;
      }

			output_file.open (path, ios::out | ios::trunc);
			if (!output_file.is_open()) {
				output_directory = "./";
				std::cerr << "Unable to open/create file. "
						"Will proceed to save cloud segments in cwd.\n";
			} else {
				OUTPUT= true;
			}
		}
	} catch(std::exception& e) {
		std::cerr << "Unhandled Exception reached the top of main: "
				<< e.what() << ", application will now exit" << std::endl;
		return(2);

	}
	search_flip(); //initialize search flags
	char * inputFile = argv[1];
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>),
			cloud_filtered (new pcl::PointCloud<pcl::PointXYZ>);
	//All objects needed for segmentation
	pcl::PointCloud<pcl::Normal>::Ptr cloud_normals (new pcl::PointCloud<pcl::Normal>);
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ> ());

	if (pcl::io::loadPCDFile<PointT> (inputFile, *cloud) == -1) //* load the file
	{
		PCL_ERROR ("Couldn't read .pcd file\n");
		return (-1);
	}

	if (VERBOSE) cout<<"Loaded point cloud...\n";
	//used to identify segmented point cloud files
	char cloud_name[strlen(inputFile) - 4 + 1]; //remove .pcd
	strncpy(cloud_name, inputFile, strlen(inputFile) - 4);

	// Create the filtering object & filter point cloud
	pcl::StatisticalOutlierRemoval<pcl::PointXYZ> sor;
	sor.setInputCloud (cloud);
	sor.setMeanK (50); //number of neighbours
	sor.setStddevMulThresh (5.0);
	sor.filter (*cloud_filtered);

	// Estimate point normals
	ne.setSearchMethod (tree);
	ne.setInputCloud (cloud_filtered);
	ne.setKSearch (50);
	ne.compute (*cloud_normals);

	vector<int> inliers;
	vector<int> outliers;
	vector<pcl::PointCloud<PointT>::Ptr> clouds;
	std::vector<std::vector<float> > models;

	if (OUTPUT){
	  output_file << inputFile << " Sequential RANSAC Models Extracted\n";
	  output_file << "File Format:\n";
	  output_file << "Plane model format:\n";
	  output_file << "\t<segment.pcd> <x> data points. Model # <x> coefficients <ax + by + cz + d = 0>\n";
	  output_file << "Cylinder model format:\n";
	  output_file << "\t<segment.pcd> <x> data points. Model # <x> coefficients Axis: point 1 <x, y, z> point 2 <x, y, z> Radius <R>\n";
	  output_file << "Sphere model format:\n";
	  output_file << "\t<segment.pcd> <x> data points. Model # <x> coefficients Center <x, y, z> Radius <R>\n";
	  output_file << "Segmentation:\n";
	}

	if (VERBOSE) cout<<"Fitting model...\n";

	pcl::PCDWriter writer;
	int j = -1;//shape number (initially incremented to 0 (PLANE) in the loop)
	int i = 1, nr_points = (int) cloud_filtered->points.size ();
	MAX_OUTLIERS= nr_points * PERCENT_OUTLIERS;
	// While 30% of the original cloud is still there
	while (cloud_filtered->points.size () > MAX_OUTLIERS & search_valid())
	{

	  j++; //change the shape
		if (j==SHAPES) j = 0; //reinitialize to the first shape
		//Skip the shapes that were eliminated due to ransac's segmentation failure
		if (j==PLANE){
			if (!search_plane) continue;
		} else if (j==CYLINDER){
			if (!search_cyl) continue;
		} else if (j==SPHERE){
			if (!search_sphere) continue;
		}

		// Segment the largest component from the remaining cloud, with fixed sequence of shapes
		pcl::ModelCoefficients::Ptr coefficients (new pcl::ModelCoefficients);
		string shape_name;

		//Segment shape
		switch(j) {
		case PLANE:
			inliers= ransacPlane(cloud_filtered, cloud_normals, coefficients,
					plane_thresh, plane_normal_weight, max_plane_iterations);
			if (inliers.size() == 0){	 search_plane = false;	}
			shape_name="plane";
			break;
		case CYLINDER:
			inliers= ransacCylinder(cloud_filtered, cloud_normals, coefficients,
					cyl_thresh, max_cyl_radius, cyl_normal_weight, max_cyl_iterations);
			if (inliers.size() == 0){	 search_cyl = false;	}
			shape_name="cylinder";
			break;
		case SPHERE:
			inliers= ransacSphere(cloud_filtered, cloud_normals, coefficients,
					sphere_thresh, sphere_normal_weight, max_sphere_iterations);
			if (inliers.size() == 0){	 search_sphere = false;	}
			shape_name="sphere";
			break;
		}

		//Skip model & cloud update if segmentation unsuccessful
		//Another one bites the dust!
		if (inliers.size() == 0) continue;

		vector<float> coeffs;
		coeffs.push_back(j); //SHAPE identifier
		for (int c = 0; c < coefficients->values.size(); c++)
			coeffs.push_back(coefficients->values[c]);
		models.push_back(coeffs);

		if (VERBOSE) cout << "PointCloud representing the " << shape_name << " component: " << inliers.size() << " data points." << endl;
		if (VERBOSE) cout << "Model #" << i << " " << shape_name <<" coefficients: " << *coefficients << std::endl;
		if (OUTPUT){
			output_file << cloud_name << "_" << shape_name << "_" << i << ".pcd "
					<< inliers.size() << " data points. "
					<< "Model #" << i << " coefficients: "
					<< *coefficients << "\n";
		}
		pcl::PointCloud<PointT>::Ptr cloud_temp;

		if (VISUAL) { //Extract inliers, update lists of cloud segments (&clouds)
			updateCloud(cloud_filtered, inliers, &clouds, cloud_filtered, cloud_normals);
			cloud_temp = clouds.back();
		} else { //Extract inliers, dump them to diff cloud, save ref in cloud_temp
			updateCloud(cloud_filtered, inliers, &cloud_temp, cloud_filtered, cloud_normals);
		}

		//Write segments to seperate files
		std::stringstream ss;
		ss << output_directory << "/" << cloud_name << "_" << shape_name << "_" << i << ".pcd";
		writer.write<PointT> (ss.str (), *cloud_temp, false);
		i++; //another model extracted
		if (VERBOSE) cout << "Fitting next model...\n";
	}

	if (OUTPUT) {
		output_file << "Number of segments found: " << i - 1 << "\n";
		if (VERBOSE) cout << "Written model parameters to file \"" << output_file_name << "\"\n";
		output_file.close();
	}
	if (VERBOSE) cout << "Number of segments found: " << i - 1 << "\n";

	// Display segments visually
	if (VISUAL) visualizeCloud(clouds, models);

	return (0);
}
