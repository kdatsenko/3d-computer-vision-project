//Title: mixed_objects.cpp
//Purpose: To implement a generalized multi-shape ransac algorithm for segmenting PCL point clouds
//Results: PCL ransac for cylinder is poor -

// Include the necessary standard libraries
#include <iostream>
#include <cstdlib>
#include <queue>
#include <cstring>
#include <cmath>
#include <string>
// Include necessary ros/pcl libraries
//#include <pcl/console/parse.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/sample_consensus/ransac.h>
#include <pcl/sample_consensus/sac_model_plane.h>
#include <pcl/sample_consensus/sac_model_sphere.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/extract_indices.h>
#include <boost/thread/thread.hpp>
#include <pcl/ModelCoefficients.h>
#include <pcl/features/normal_3d.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/filters/passthrough.h>
#include <pcl/filters/statistical_outlier_removal.h>

// Include our own header files
#include <point_types.h>

// Using declarations
using std::cout;
using std::vector;
using std::queue;
using std::flush;
using std::string;

// Set constants
int CLUSTERS= 0;
int MAX_OUTLIERS;
const float PERCENT_OUTLIERS= .1; //MODIFY

//I/O
bool VERBOSE = false;
bool VISUAL = false;
int TYPE_VISUALIZATION = 0;
bool OUTPUT = false;

//Segment these models
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

/*
const float PERCENT_OUTLIERS= .025;
const float THRESHOLD= 0.0025; //0.0025
const float CYL_THRESH_MOD= 5;
const int MAX_CYL_ITERATIONS= 10000;

const int SHAPES= 3;
bool TREE= true;
bool NEIGHBORS= true;
bool SEARCH_PLANE = true;
bool SEARCH_CYL = true;
bool SEARCH_SPHERE = false;
 */


// A tree structure to keep track of the shape sequences ----------------------------------------------------------
struct Tree {
	bool isRoot;
	int shape; //0 - plane; 1 - cylinder; 2 - sphere;
	vector<int>* outliers;
	int score;
	struct Tree* parent;
	vector<struct Tree*>* children;
};

typedef struct Tree treeNode;

//Returns a set of outliers given a cloud and a set of inliers
vector<int>* getOutliers(pcl::PointCloud<PointT>::Ptr cloud, vector<int> inliers){
	vector<int>* outliers= new vector<int>;
	int index=0;
	for (int i=0; i<(int)(cloud->width * cloud->height); i++){
		if (i==inliers[index]){
			index++;
		} else {
			outliers->push_back(i);
		}
	}
	return outliers;
}

//Print a shape given its integer encoding
void printShape(int encoded){
	if (encoded==PLANE) {
		cout << "plane, ";
	} else if (encoded==CYLINDER){
		cout << "cylinder, ";
	} else if (encoded==SPHERE){
		cout << "sphere, ";
	}
}

// Main method
int main(int argc, char** argv) {

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
                ("cyl_thresh,a", opt::value<float>(&cyl_thresh)->default_value(0.02, "0.02"), "Cylinder Distance threshold - units are relative to data pcd file, usually metres.")
                ("max_cyl_iterations,b", opt::value<int>(&max_cyl_iterations)->default_value(5000), "Maximum number of iterations to try to find a Cylinder model")
                ("max_cyl_radius,c", opt::value<float>(&max_cyl_radius)->default_value(0.1, "0.1"), "Cylinder")
                ("cyl_normal_weight,f", opt::value<float>(&cyl_normal_weight)->default_value(0.1, "0.1"), "Set the relative weight (between 0 and 1) to give to the angular distance (0 to pi/2) between point normals and the Cylinder normal.")
                ("plane_thresh,l", opt::value<float>(&plane_thresh)->default_value(0.03, "0.03"), "Plane Distance threshold - units are relative to data pcd file, usually metres")
                ("plane_normal_weight,m", opt::value<float>(&plane_normal_weight)->default_value(0.1, "0.1"), "Set the relative weight (between 0 and 1) to give to the angular distance (0 to pi/2) between point normals and the Plane normal.")
                ("max_plane_iterations,n", opt::value<int>(&max_plane_iterations)->default_value(500), "Maximum number of iterations to try to find a Planar model")
                ("sphere_thresh,x", opt::value<float>(&sphere_thresh)->default_value(0.01, "0.01"), "Sphere Distance threshold - units are relative to data pcd file, usually metres")
                ("sphere_normal_weight,y", opt::value<float>(&sphere_normal_weight)->default_value(0.1, "0.1"), "Set the relative weight (between 0 and 1) to give to the angular distance (0 to pi/2) between point normals and the Sphere normal.")
                ("max_sphere_iterations,z", opt::value<int>(&max_sphere_iterations)->default_value(100), "Maximum number of iterations to try to find a Spherical model");
    opt::variables_map vm;
		try
		{	// Parse arguments and, if necessary, print usage
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
				cout << "\nGeneralized Minimal-Model Ransac Segmentation" << endl
						<< "\nUsage: mixed_objects_tree <input_cloud.pcd> [OPTIONS]" << endl
						<< "Purpose: Uses RANSAC algorithm to segment a PCL pointcloud (.pcd) file into "
						"multiple objects. Supported shapes are planes, cylinders, and spheres. "
						"Computes and stores all combinations of the sequence in which models are segmented, and uses BFS to find the minimal set of "
						"models that encompass ~95% of the points. Outputs the segmented .pcd plane cloud files." << endl
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
		return(2);

	}

	// Initialize queue of nodes containing segment outliers to process.
	// Useful for BFS of shortest branch of segments.
	queue<treeNode*> toCompute;
	search_flip(); //initialize search flags
	// initialize PointClouds
	pcl::PointCloud<PointT>::Ptr cloud_unfil (new pcl::PointCloud<PointT>);
	pcl::PointCloud<PointT>::Ptr updated (new pcl::PointCloud<PointT>);
	pcl::PointCloud<PointT>::Ptr cloud (new pcl::PointCloud<PointT>);

	char * inputFile = argv[1];

	if (pcl::io::loadPCDFile<PointT> (inputFile, *cloud_unfil) == -1) //* load the file
	{
		PCL_ERROR ("Couldn't read .pcd file\n");
		return (-1);
	}

	if (VERBOSE) cout<<"Loaded point cloud...\n";
	//used to identify segmented point cloud files
	char cloud_name[strlen(inputFile) - 4 + 1]; //remove .pcd
	strncpy(cloud_name, inputFile, strlen(inputFile) - 4);

	if (VERBOSE) cout << "PointCloud before filtering has: " << cloud_unfil->size() << " data points.\n";

	// Filter point cloud for spurious NaNs
	// Create the filtering object
	pcl::StatisticalOutlierRemoval<pcl::PointXYZ> sor;
	sor.setInputCloud (cloud_unfil);
	sor.setMeanK (50); //number of neighbours
	sor.setStddevMulThresh (5.0);
	sor.filter (*cloud);

	/* PASSTHROUGH filter (more severe - downsample of point greater)
	pcl::PassThrough<PointT> pass;
	// Build a passthrough filter to remove spurious NaNs
	pass.setInputCloud (cloud_unfil);
	pass.setFilterFieldName ("z");
	pass.setFilterLimits (0, 1.5);
	pass.filter (*cloud);
	 */

	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ> ());
	pcl::PointCloud<pcl::Normal>::Ptr cloud_normals (new pcl::PointCloud<pcl::Normal>);
	//Compute the normals
	pcl::NormalEstimation<PointT, pcl::Normal> ne;
	ne.setSearchMethod(tree);
	ne.setInputCloud(cloud);
	ne.setKSearch(50);
	ne.compute(*cloud_normals);
	if (VERBOSE) cout << "Computed normals...\n";

	vector<int> inliers; //List of model inlier point indexes
	pcl::ModelCoefficients::Ptr coeffs (new pcl::ModelCoefficients);
	pcl::PointCloud<PointT>::Ptr cloud_temp;

	//REMOVE Ground Plane
	inliers = ransacPlane(cloud, cloud_normals, coeffs); //0.025
	//outlier_points = (getOutliers(cloud, inliers));
	if (VERBOSE) cout << "Removed ground planes...\n";
	updateCloud(cloud, inliers, &cloud_temp, cloud, cloud_normals);

//updateCloud(cloud_filtered, inliers, &cloud_temp, cloud_filtered, cloud_normals);

	/* View points left over
	pcl::PCDWriter writer;
	std::stringstream ss;
	ss << ".../filtered_cloud.pcd";
	writer.write<PointT> (ss.str (), *cloud, false);*/

	if (VERBOSE) cout << "PointCloud after filtering has: " << cloud->size() << " data points.\n";

	//Set the maximum outliers allowed before ending segmentation process
	MAX_OUTLIERS= (int)(cloud->size() * PERCENT_OUTLIERS);
	if (VERBOSE) cout<<"Loaded point cloud...\n";

	// Instantiate the tree root node
	treeNode* root= new treeNode;
	root->children= new vector<treeNode*>;
	vector<int>* cloudIndices= new vector<int>;
	for (int i=0; i<(int)(cloud->width * cloud->height); i++){ cloudIndices->push_back(i); }
	root->outliers= cloudIndices;
	root->score= root->outliers->size();
	root->isRoot= true;
	if (VERBOSE) cout << "Instantiated parse tree...\n";


	vector<int> outliers; //List of model outlier point indexes
	vector<pcl::PointCloud<PointT>::Ptr> clouds;
	std::vector<std::vector<float> > models;
	pcl::ModelCoefficients::Ptr coefficients (new pcl::ModelCoefficients);
	pcl::copyPointCloud<PointT>(*cloud, *cloudIndices, *updated);
	treeNode* result=0;
	bool stillRunning= true;
	if (VERBOSE) cout << "Built working cloud...\n";

	treeNode* workingNode; //Current node being segmented
	//The normals of the points of current working node
	pcl::PointCloud<pcl::Normal>::Ptr working_normals (new pcl::PointCloud<pcl::Normal>);

	pcl::PCDWriter writer;

	//Add the root treeNode to the queue
	toCompute.push(root);
	if (VERBOSE) cout<<"Fitting model";
	// Iterate over ransac models
	while (toCompute.size() > 0 && stillRunning && search_valid()){
		if (VERBOSE) cout << "." << flush;
		// Pop off the head of the queue
		workingNode= toCompute.front();
		toCompute.pop();

		// Build cloud of points remaining to be fit
		pcl::copyPointCloud<PointT>(*cloud, *(workingNode->outliers), *updated);
		pcl::copyPointCloud<pcl::Normal>(*cloud_normals, *(workingNode->outliers), *working_normals);

		//workingNode->outliers->clear();
		//delete workingNode->outliers;
		//cout << "Outliers: " << *(workingNode->outliers) << "\n";

		// Compute children of this node, and add them to the queue
		for (int i=0; i<SHAPES; i++){
			//Skip the shapes that were eliminated due to ransac's segmentation failure
			if (i==PLANE){
				if (!search_plane) continue;
			} else if (i==CYLINDER){
				if (!search_cyl) continue;
			} else if (i==SPHERE){
				if (!search_sphere) continue;
			}
			//Build a child node with a new shape
			treeNode* childNode= new treeNode;
			childNode->isRoot= false;
			childNode->shape= i;
			childNode->parent= workingNode;
			workingNode->children->push_back(childNode);
			childNode->children= new vector<treeNode*>;


			//Fit a ransac shape and compute remaining outliers
			if (childNode->shape==PLANE){
				inliers= ransacPlane(updated, working_normals, coefficients,
						plane_thresh, plane_normal_weight, max_plane_iterations);
				if (inliers.size() == 0)
					search_plane = false;

			} else if (childNode->shape==CYLINDER) {
				inliers= ransacCylinder(updated, working_normals, coefficients,
						cyl_thresh, max_cyl_radius, cyl_normal_weight, max_cyl_iterations);
				if (inliers.size() == 0)
					search_cyl = false;

			} else if (childNode->shape==SPHERE) {
				inliers= ransacSphere(updated, working_normals, coefficients,
						sphere_thresh, sphere_normal_weight, max_sphere_iterations);
				if (inliers.size() == 0)
					search_sphere = false;
			}

			//If this is the end of segmentation, the end node of the shortest path is the
			//one we just popped off of the queue
			if (!search_valid()){
				result = workingNode;
				stillRunning = false;
				break;
			} else if (inliers.size() == 0){
			    continue;
			} else { //Assign the outliers to node representing this segment
				childNode->outliers= getOutliers(updated, inliers);
				childNode->score= childNode->outliers->size();
			}
			// Push child onto queue, check score
			if (childNode->score < MAX_OUTLIERS){
				stillRunning= false;
				result= childNode;
				break;
			} else {
				toCompute.push(childNode);
				if (VERBOSE) cout<<"push\n";
			}

		} //End of shape iteration

	} //End of main loop

	// Cascade back to find the path of shapes
	vector<int> reversePath;
	while (!result->isRoot){
		reversePath.push_back(result->shape);
		result= result->parent;
	}

	if (OUTPUT){
		output_file << inputFile << " Mixed RANSAC Models Extracted\n";
    output_file << "File Format:\n";
    output_file << "Plane model format:\n";
    output_file << "\t<segment.pcd> <x> data points. Model # <x> coefficients <ax + by + cz + d = 0>\n";
    output_file << "Cylinder model format:\n";
    output_file << "\t<segment.pcd> <x> data points. Model # <x> coefficients Axis: point 1 <x, y, z> point 2 <x, y, z> Radius <R>\n";
    output_file << "Sphere model format:\n";
    output_file << "\t<segment.pcd> <x> data points. Model # <x> coefficients Center <x, y, z> Radius <R>\n";
    output_file << "Segmentation:\n";
	}
	//Compute series of fit clouds
	pcl::copyPointCloud<PointT>(*cloud, *cloudIndices, *updated);
	pcl::copyPointCloud<pcl::Normal>(*cloud_normals, *cloudIndices, *working_normals);
	cout << "Shapes Fit: ";
	for (int i=reversePath.size()-1; i>=0; i--){
	  pcl::ModelCoefficients::Ptr coefficients_f (new pcl::ModelCoefficients);
		printShape(reversePath[i]);
		if (reversePath[i]==1){
			inliers= ransacPlane(updated, working_normals, coefficients_f,
					plane_thresh, plane_normal_weight, max_plane_iterations);
		} else if (reversePath[i]==0){
			inliers= ransacCylinder(updated, working_normals, coefficients_f,
					cyl_thresh, max_cyl_radius, cyl_normal_weight, max_cyl_iterations);
		} else if (reversePath[i]==2){
			inliers= ransacSphere(updated, working_normals, coefficients_f,
					sphere_thresh, sphere_normal_weight, max_sphere_iterations);
		}
		updateCloud(updated, inliers, &clouds, updated, working_normals);
		CLUSTERS++;

		//Optional
		if (VISUAL) {
			vector<float> coeffs;
			coeffs.push_back(reversePath[i]); //SHAPE identifier
			for (int j = 0; j < coefficients_f->values.size(); j++)
				coeffs.push_back(coefficients_f->values[j]);
			models.push_back(coeffs);
		}
		if (OUTPUT){
			output_file << cloud_name << "_" << shape_names[reversePath[i]] << "_" << CLUSTERS << ".pcd "
					<< inliers.size() << " data points. "
					<< "Model #" << i << " coefficients: "
					<< *coefficients_f << "\n";
		}

	}

	cout << "\n";

	// Print the number of iterations
	cout << "Number of clusters found: " << CLUSTERS << "\n";
	if (OUTPUT) {
		output_file << "Number of segments found: " << CLUSTERS << "\n";
		if (VERBOSE) cout << "Written model parameters to file \"" << output_file_name << "\"\n";
		output_file.close();
	}

	if (VERBOSE) cout << "Building segmented clouds..." << endl;
	// Build cloud segments

	int shape_index = reversePath.size() - 1;
	for (int i=0; i<(int)clouds.size(); i++){
		std::stringstream ss;
		ss << output_directory << "/" << cloud_name << "_"
				<< shape_names[reversePath[shape_index]] << "_" << i << ".pcd";
		writer.write<PointT> (ss.str (), *clouds[i], false);
		shape_index--;
	}

	//Save points that were identified as not belonging to any model, for comparision.
	pcl::PointCloud<PointT>::Ptr cloud_f (new pcl::PointCloud<PointT>);
	std::stringstream ss2;
	vector<int> * outlier_points = (getOutliers(updated, inliers));
	pcl::copyPointCloud<PointT>(*updated, *(outlier_points), *cloud_f);
	ss2 << output_directory << "/" << cloud_name << "_" << "left_over.pcd";
	writer.write<PointT> (ss2.str (), *cloud_f, false);

	// Display segments visually
	if (VISUAL) visualizeCloud(clouds, models);

	return 0;
}
