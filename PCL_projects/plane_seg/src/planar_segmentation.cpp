#include <iostream>
#include <fstream>
#include <pcl/ModelCoefficients.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include "boost/program_options.hpp"
#include <pcl/common/common.h>
#include <pcl/common/transforms.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/filters/extract_indices.h>
#include <point_types.h>

bool VERBOSE = false;
bool VISUAL = false;
int TYPE_VISUALIZATION = 0;
bool OUTPUT = false;

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
								("threshold,t", opt::value<double>(&threshold)->default_value(0.01), "Set error threshold to model");

		opt::variables_map vm;
		try
		{	// Parse arguments and, if necessary, print usage
			opt::store(opt::parse_command_line(argc, argv, desc), vm); // can throw

			if (argc < 2){
				std::cerr << desc << std::endl;
				exit(0);
			}

			/** --help option */
			if ( vm.count("help")  ){
				cout << "\nMulti-Plane Ransac Segmentation" << endl
						<< "\nUsage: plane_segmentation <input_cloud.pcd> [OPTIONS] [-o <output_file>] [-d <./relative/path/to/dir>]" << endl
						<< "Purpose: Uses RANSAC algorithm to segment a PCL pointcloud (.pcd) file into "
						"multiple planes. Outputs the segmented .pcd plane cloud files." << endl
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
		}
		if (vm.count("outdir")){ //create directory if necessary
		  boost::filesystem::path dir(output_directory);
		  if(!boost::filesystem::exists(dir)){
		    if (boost::filesystem::create_directory(dir)){
		      cout << "Directory " << output_directory << " was created." << std::endl;
		    } else {
		      std::cerr << "Directory invalid. Will proceed to save output to cwd." << std::endl;
		      output_directory = "./"; //reset to cwd.
		    }
		  }
		} if (vm.count("output")){ // Check output flag
		  string p = output_directory + "/" + output_file_name;
		  const char * path = p.c_str();
		  output_file.open (path, ios::out | ios::trunc);
		  if (!output_file.is_open()) {
		    output_directory = "./";
		    std::cerr << "Unable to open/create file. "
		        "Will skip output file, but proceed to save cloud segments in cwd.\n";
		  } else {
		    OUTPUT= true;
		  }
		}
	} catch(std::exception& e) {
		std::cerr << "Unhandled Exception reached the top of main: "
				<< e.what() << ", application will now exit" << std::endl;
		return 2;

	}

	char* inputFile = argv[1];
	vector<pcl::PointCloud<PointT>::Ptr> clouds;
	std::vector<std::vector<float> > models;

	pcl::PointCloud<PointT>::Ptr cloud(new pcl::PointCloud<PointT>),
			cloud_f (new pcl::PointCloud<PointT>);

	if (pcl::io::loadPCDFile<PointT> (inputFile, *cloud) == -1) //* load the file
	{
		PCL_ERROR ("Couldn't read .pcd file\n");
		return (-1);
	}

	if (VERBOSE) cout<<"Loaded point cloud...\n";
	//used to identify segmented point cloud files
	char cloud_name[strlen(inputFile) - 4 + 1]; //remove .pcd
	strncpy(cloud_name, inputFile, strlen(inputFile) - 4);

	pcl::ModelCoefficients::Ptr coefficients (new pcl::ModelCoefficients);
	pcl::PointIndices::Ptr inliers (new pcl::PointIndices);
	// Create the segmentation object
	pcl::SACSegmentation<PointT> seg;
	// Optional
	seg.setOptimizeCoefficients (true);
	// Mandatory
	seg.setModelType (pcl::SACMODEL_PLANE);
	seg.setMethodType (pcl::SAC_RANSAC);
	seg.setDistanceThreshold (threshold);

	// Create the filtering object
	pcl::ExtractIndices<PointT> extract;
	pcl::PCDWriter writer;

	if (OUTPUT){
		output_file << inputFile << " Planar Models Extracted\n";
		output_file << "File Format:\n";
		output_file << "<segment.pcd> <x> data points. Model # <x> coefficients <ax + by + cz + d = 0>\n";
		output_file << "Segmentation:\n";
	}

	if (VERBOSE) cout<<"Fitting model...\n";
	int i = 1, nr_points = (int) cloud->points.size ();
	// While 30% of the original cloud is still there
	while (cloud->points.size () > 0.05 * nr_points) //noise is estimated to be 5%
	{
		// Segment the largest planar component from the remaining cloud
		seg.setInputCloud (cloud);
		seg.segment (*inliers, *coefficients);

		if (inliers->indices.size () == 0)
		{
			if (VERBOSE) std::cerr << "Could not estimate a planar model for the given dataset." << std::endl;
			break;
		}
		if (VERBOSE) {
			cout << "Model #" << i << " coefficients: " << coefficients->values[0] << " "
					<< coefficients->values[1] << " "
					<< coefficients->values[2] << " "
					<< coefficients->values[3] << std::endl;
		}

		// Extract the inliers
		// copies all inliers of the model computed to another PointCloud
		pcl::PointCloud<PointT>::Ptr cloud_temp (new pcl::PointCloud<PointT>);
		extract.setInputCloud (cloud);
		extract.setIndices (inliers);
		extract.setNegative (false);
		extract.filter (*cloud_temp);
		if (VERBOSE) cout << "PointCloud representing the planar component: " << cloud_temp->width * cloud_temp->height << " data points." << std::endl;

		std::stringstream ss;
		ss << output_directory << "/" << cloud_name << "_plane_" << i << ".pcd";
		writer.write<PointT> (ss.str (), *cloud_temp, false);

		if (VISUAL){
			//saves the segment cloud
			clouds.push_back(cloud_temp);
			//saves the model coefficients
			vector<float> plane_coeffs;
			plane_coeffs.push_back(PLANE); //identify shape
			for (int j = 0; j < coefficients->values.size(); j++)
				plane_coeffs.push_back(coefficients->values[j]);
			models.push_back(plane_coeffs);
		}

		if (OUTPUT){
			output_file << cloud_name << "_plane_" << i << ".pcd "
					<< cloud_temp->width * cloud_temp->height << " data points. "
					<< "Model #" << i << " coefficients: "
					<< coefficients->values[0] << " "
					<< coefficients->values[1] << " "
					<< coefficients->values[2] << " "
					<< coefficients->values[3] << "\n";
		}

		// Create the filtering object
		extract.setNegative (true);
		extract.filter (*cloud_f);
		cloud.swap (cloud_f);
		i++;
	}
	if (OUTPUT) {
		output_file << "Number of Planar segments found: " << i - 1 << "\n";
		if (VERBOSE) cout << "Written model parameters to file \"" << output_file_name << "\"\n";
		output_file.close();
	}
	if (VERBOSE) cout << "Number of Planar segments found: " << i - 1 << "\n";

	// Display segments visually
	if (VISUAL) visualizeCloud(clouds, models);

	return (0);
}
