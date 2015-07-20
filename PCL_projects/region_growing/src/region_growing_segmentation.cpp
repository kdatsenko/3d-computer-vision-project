#include <iostream>
#include <vector>
#include <pcl/common/geometry.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/search/search.h>
#include <pcl/search/kdtree.h>
#include <pcl/features/normal_3d.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/filters/passthrough.h>
#include <pcl/segmentation/region_growing.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <string>
// Include our own header file's
#include <point_types.h>
#include <vtkLODActor.h>
#include <boost/heap/fibonacci_heap.hpp>
#include <helpers.h>

//I/O
bool VERBOSE = true;
bool VISUAL = true;
int TYPE_VISUALIZATION = 0;
bool OUTPUT = false;

using std::vector;

struct Tree {
  bool isBase;
  int score;
  int id;
  int cluster_index; //if isBase=true
  //record if it exists, in case for split/merge
  std::list<pcl::PointCloud<PointT>::Ptr>::iterator cloud_ref;
  //pcl::PointCloud<PointT>::Ptr * cloud_ref;
  //int cloud_index; Deprecated
  struct Tree* parent; //0 or 1
  vector<struct Tree*>* children; //0 or 2

};
typedef struct Tree treeNode;

struct Merge_Heap {
  int N1;
  int N2;
  int i_edge;
  float cost;
  Merge_Heap(int N1, int N2, int i_edge, int cost)
    : N1(N1) // initializes first node
    , N2(N2) // initializes second node
    , i_edge(i_edge) // initializes index of edge repr in EDGE table
    , cost(cost)
  {}
};

class PCL_EXPORTS RegionGrowing : public pcl::RegionGrowing<PointT, pcl::Normal>{
  public:
    using pcl::RegionGrowing<PointT, pcl::Normal>::point_neighbours_;
    using pcl::RegionGrowing<PointT, pcl::Normal>::point_labels_;

    void extract (std::vector <pcl::PointIndices>& clusters)
    {
      clusters_.clear ();
      clusters.clear ();
      point_neighbours_.clear ();
      point_labels_.clear ();
      num_pts_in_segment_.clear ();
      number_of_segments_ = 0;

      bool segmentation_is_possible = initCompute ();
      if ( !segmentation_is_possible )
      {
        deinitCompute ();
        return;
      }

      segmentation_is_possible = prepareForSegmentation ();
      if ( !segmentation_is_possible )
      {
        deinitCompute ();
        return;
      }

      findPointNeighbours ();
      applySmoothRegionGrowingAlgorithm ();
      assembleRegions ();

      clusters.resize (number_of_segments_);
      std::vector<pcl::PointIndices>::iterator cluster_iter_input = clusters.begin ();
      for (std::vector<pcl::PointIndices>::const_iterator cluster_iter = clusters_.begin (); cluster_iter != clusters_.end (); cluster_iter++)
      {
        if ((cluster_iter->indices.size () >= min_pts_per_cluster_) &&
            (cluster_iter->indices.size () <= max_pts_per_cluster_))
        {
          *cluster_iter_input = *cluster_iter;
          cluster_iter_input++;
        }
      }


      clusters_ = std::vector<pcl::PointIndices> (clusters.begin (), cluster_iter_input);
      deinitCompute ();
    }

    void assembleRegions ()
    {
      int number_of_segments = static_cast<int> (num_pts_in_segment_.size ());
      int number_of_points = static_cast<int> (input_->points.size ());

      pcl::PointIndices segment;
      clusters_.resize (number_of_segments, segment);

      vector<int> segment_indices;
      segment_indices.resize(number_of_segments, 0);
      int count_valid_seg = 0;

      for (int i_seg = 0; i_seg < number_of_segments; i_seg++)
      {
        clusters_[i_seg].indices.resize ( num_pts_in_segment_[i_seg], 0);

        if ((num_pts_in_segment_[i_seg] >= min_pts_per_cluster_) &&
            (num_pts_in_segment_[i_seg] <= max_pts_per_cluster_)){
          segment_indices[i_seg] = count_valid_seg;
          count_valid_seg++;
        } else {
          segment_indices[i_seg] = -1;
        }

      }

      std::vector<int> counter;
      counter.resize (number_of_segments, 0);

      for (int i_point = 0; i_point < number_of_points; i_point++) //for each point
      {
        //if num_pts_in_segment_[i_seg] is approp,
        int segment_index = point_labels_[i_point];
        if (segment_index != -1)
        {
          int point_index = counter[segment_index];
          clusters_[segment_index].indices[point_index] = i_point;
          counter[segment_index] = point_index + 1;
        }
        point_labels_[i_point] = segment_indices[segment_index];

      }

      number_of_segments_ = count_valid_seg;
    }



};


typedef struct Merge_Heap n_pair;

/* Min-Heap compare function. */
struct compare_node
{
    bool operator()(const n_pair& n1, const n_pair& n2) const
    {
        return n1.cost > n2.cost;
    }
};




//handle_type handle_t;


// clusters from base level region growing
std::vector <pcl::PointIndices> clusters;

//Heirarchy of merges, used to reconstruct point cloud segmentation at any level.
vector< pair<treeNode*, treeNode*> > merge_order; //vector => should be able to expand

//next merge operation
vector< pair<treeNode*, treeNode*> >::iterator merge_op; // = merge_order.begin();
//Rules of use:
//cout << *it << end; // output value with dereference operator

//list of current level clouds
std::list<pcl::PointCloud<PointT>::Ptr> cloud_list;

//Unorganized initial input point cloud
pcl::PointCloud<PointT>::Ptr cloud;

int NUM_EDGES;

//EDGE TABLE doubly-linked pointer indices
const int NODE_ONE_NEXT = 2;
const int NODE_ONE_PREC = 3;
const int NODE_TWO_NEXT = 4;
const int NODE_TWO_PREC = 5;


typedef boost::heap::fibonacci_heap<n_pair, boost::heap::compare<compare_node> > fibonacci_heap;
typedef typename boost::heap::fibonacci_heap<n_pair, boost::heap::compare<compare_node> >::handle_type handle_t;
/* Cost function type for conciseness. */
typedef float (*Cost_Function)(pcl::PointCloud<PointT>&, pcl::PointCloud<PointT>&);

/* Change cost computation function HERE.
 * Used in wrapper cost function calcCost(). */
Cost_Function compute_cost = compute;

/*PROTOTYPES*/

int getNodeNextColumn(int edge, int node, int edge_map[][6]);

void mergeEdge (//HEAP,
                   int i_edge,
                   handle_t heap_handler_list[],
                   int edge_map[][6],
                   int node_list[],
                   fibonacci_heap& heap);

void removeNode (int node, int edge, int edge_map[][6],
                 int node_list[],
                 int link_fwd= -1, int link_bck= -1);

void findSegmentNeighbours(int edge_map[][6],
                           int node_list[],
                           //vector<vector<int>>& segment_neighbours,
                           RegionGrowing& reg,
                           fibonacci_heap& heap,
                           handle_t heap_handler_list[],
                           treeNode** tree_nodes,
                           float max_dist);

void findRegionsKNN (int index,
                     RegionGrowing& reg,
                     int edge_map[][6],
                     int& edge_count,
                     int **edge_hash_tab,
                     int node_list[],
                     //std::vector<int>& nghbrs,
                     fibonacci_heap& heap,
                     handle_t heap_handler_list[],
                     treeNode** tree_nodes,
                     float max_dist);
vector<int> getCloudIndices(treeNode * cloud);
string constructCloudLabel (int id);
void getPointCloudRenderingProperties_Color (
    boost::shared_ptr<pcl::visualization::PCLVisualizer>& viewer,
    string label, double& r, double& g, double &b);
void keyboardEventOccurred (const pcl::visualization::KeyboardEvent &event,
                            void* viewer_void);
int getNumEdges(RegionGrowing& reg, float max_dist);
float calcCost(int N1, int N2, treeNode** tree_nodes);
float getDistance(PointT& p1, PointT& p2);


//Return a set of outliers given a cloud and a set of inliers
vector<int> unionIndices(vector<int>& U_inliers, vector<int>& inliers){
	vector<int> result;
	int index=0;
	int size = inliers.size();
	int i = 0;
	while (i < U_inliers.size() && index < size){
		if (U_inliers[i] < inliers[index]){
			result.push_back(U_inliers[i]);
			i++;
		} else {
			result.push_back(inliers[index]);
			index++;
		}
	}
	if (U_inliers.size() != i){
		for (; index != size; index++)
			result.push_back(inliers[index]);
	} else {
		for (; i < U_inliers.size(); i++)
			result.push_back(U_inliers[i]);
	}
	return result;
}


//Return a set of outliers given a cloud and a set of inliers
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

void keyboardEventOccurred (const pcl::visualization::KeyboardEvent &event,
                            void* viewer_void);

int
main (int argc, char** argv)
{


  /*float sm_threshold = 0.037;
    float c_threshold = 1.0;
    float max_dist = 5.0;*/

  ofstream output_file;
  string output_file_name;
  string output_directory;
  float sm_threshold;
  float c_threshold;
  float max_dist;

  try {

    namespace opt = boost::program_options;

    opt::options_description desc("Options");
    desc.add_options()
                  ("help,h", "Print usage options.")
                ("verbose,v", "Verbose mode, print out progress during computation.")
                ("visual,s", "(Type=RGB colour) Visual mode, display a colored segmentation of the cloud.")
                ("output,o", opt::value<string>(&output_file_name), "<output_file> Save the model parameters to this file.")
                ("outdir,d", opt::value<string>(&output_directory)->default_value("./"), "<./directory> Save all output to this RELATIVE (to location of execution file) dir.")
                ("smoothness threshold,t", opt::value<float>(&sm_threshold)->default_value(2.0 / 180.0 * M_PI), "Set smoothness threshold to region")
                ("curvature threshold,c", opt::value<float>(&c_threshold)->default_value(1.0), "Set curvature threshold to region")
                ("max distance,m", opt::value<float>(&max_dist)->default_value(5), "Set maximum distance between regions considered to be neighbours.");


    opt::variables_map vm;
    try
    { // Parse arguments and, if necessary, print usage
      opt::store(opt::parse_command_line(argc, argv, desc), vm); // can throw

      if (argc < 2){
        std::cerr << desc << std::endl;
        exit(0);
      }

      /** --help option */
     if ( vm.count("help")  ){
        cout << "\nMulti-Scale Bottom-up Region Growing Segmentation" << endl
            << "\nUsage: region_growing_segmentation <input_cloud.pcd> [OPTIONS] [-o <output_file>] [-d <./relative/path/to/dir>]" << endl
            << "Purpose: Uses Region-growing to segment a PCL pointcloud (.pcd) file into "
            "multiple segments. Outputs the segmented .pcd plane cloud files." << endl
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
      TYPE_VISUALIZATION = VISUALIZATION_CUSTOM;
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

  //Initialize new point cloud, save pointer.
  cloud = pcl::PointCloud<pcl::PointXYZ>::Ptr (new pcl::PointCloud<pcl::PointXYZ>);
  cloud.get()->points.clear();
  pcl::PointCloud<PointT>::Ptr cloud_unfiltered (new pcl::PointCloud<PointT>);

  if (pcl::io::loadPCDFile<PointT> (inputFile, *cloud_unfiltered) == -1) //* load the file
  {
    PCL_ERROR ("Couldn't read .pcd file\n");
    return (-1);
  }

  if (VERBOSE) cout<<"Loaded point cloud...\n";
  //used to identify segmented point cloud files
 // char cloud_name[strlen(inputFile) - 4 + 1]; //remove .pcd
 // strncpy(cloud_name, inputFile, strlen(inputFile) - 4);

  pcl::search::Search<PointT>::Ptr tree = boost::shared_ptr<pcl::search::Search<PointT> > (new pcl::search::KdTree<PointT>);
  pcl::PointCloud <pcl::Normal>::Ptr cloud_normals (new pcl::PointCloud <pcl::Normal>);
  pcl::PCDWriter writer;

  if (VERBOSE) cout << "PointCloud before filtering has: " << cloud_unfiltered->points.size ()  << " data points." << std::endl;

  // Filter point cloud for spurious NaNs
  // Create the filtering object
  pcl::StatisticalOutlierRemoval<pcl::PointXYZ> sor;
  sor.setInputCloud (cloud_unfiltered);
  sor.setMeanK (50); //number of neighbours
  sor.setStddevMulThresh (5.0);
  sor.filter (*cloud);

  //PASSTHROUGH filter (more severe - downsample of point greater)
  /*pcl::PassThrough<PointT> pass;
  // Build a passthrough filter to remove spurious NaNs
  pass.setInputCloud (cloud_unfiltered);
  pass.setFilterFieldName ("z");
  pass.setFilterLimits (0.0, 5.0);
  pass.filter (*cloud);*/


  pcl::NormalEstimation<PointT, pcl::Normal> normal_estimator;
  normal_estimator.setSearchMethod (tree);
  normal_estimator.setInputCloud (cloud);
  normal_estimator.setKSearch (50);
  normal_estimator.compute (*cloud_normals);



  if (VERBOSE) cout << "PointCloud after filtering has: " << cloud->points.size ()  << " data points." << std::endl;

  RegionGrowing reg;
  reg.setMinClusterSize (5);
  reg.setMaxClusterSize (100000);
  reg.setSearchMethod (tree);
  reg.setNumberOfNeighbours (20);
  reg.setInputCloud (cloud);
  //reg.setIndices (indices);
  reg.setInputNormals (cloud_normals);
  reg.setSmoothnessThreshold (sm_threshold);
  reg.setCurvatureThreshold (c_threshold);
  reg.extract (clusters);

  if (VERBOSE) cout <<"Initial region growth is finished.\n";
  if (VERBOSE) cout << "Number of clusters is equal to " << clusters.size () << std::endl;

  int x = 0;
  for (int i = 0; i < clusters.size(); i++){
    x += clusters[i].indices.size();
  }
  std::cout << "Clustered points size: " << x << std::endl;

  // Create the filtering object
  pcl::ExtractIndices<PointT> extract;
  vector<pcl::PointCloud<PointT>::Ptr> clouds;

  int cloud_id = clusters.size(); //initialize identifier for cloud label beyond base-level regions.

  /* For each node number (N1, N2, etc) derived from the corresponding indices of the base-level region
   * stored in the point indices array 'clusters', tree_nodes stores a pointer to its representational node.  */
  // treeNode* tree_nodes[clusters.size()]; //Size = number of nodes (n1, n2, n3, etc)

  treeNode** tree_nodes = new treeNode*[clusters.size()]; /* clusters.size() ^ 2 array providing the row index of the edge for (NodeX, NodeY) */

  //Cloud creation
  if (VERBOSE) cout<<"Creating the Visualization cloud list.\n";
  for (int i = 0; i < clusters.size(); i++){
    // Extract the inliers
    // copies all inliers of the model computed to another PointCloud
    pcl::PointCloud<PointT>::Ptr cloud_temp (new pcl::PointCloud<PointT>);
    pcl::ExtractIndices<PointT> extract;
    extract.setInputCloud (cloud);

    pcl::PointIndicesPtr cluster_ptr(new pcl::PointIndices);
    cluster_ptr->indices = clusters[i].indices;

    extract.setIndices (cluster_ptr);
    extract.setNegative (false);
    extract.filter (*cloud_temp);
    clouds.push_back(cloud_temp);
    cloud_list.push_back(cloud_temp);

    //Create node for this cluster:
    treeNode* node = new treeNode();
    node->isBase = true;
    node->id = i;
    node->cluster_index = i;
    //&cloud_list.back() <-- for when you need to save ref
    node->cloud_ref = --cloud_list.end(); //Save reference to corres. PointCloud in cloud_list.
    //cloud_list.
    pcl::PointCloud<PointT>::Ptr cl = *(node->cloud_ref);
    //no children, no parent yet!
    tree_nodes[i] = node;
  }


  NUM_EDGES = getNumEdges(reg, max_dist); //Find total number of neighbouring edges between segments.
  cout << "NUM_EDGES: " << NUM_EDGES << "\n";
  fibonacci_heap heap; //Min-heap for merge procedure organized by decreasing cost of merge btw segments.
  handle_t heap_handler_list[NUM_EDGES]; //For each edge index in Edge Map, links to its place in the heap.

  /* The graph of segments and their neighbours is represented by node_list and edge_map as follows:*/

  /* For each node (separate cloud), provides a link to row in the Edge Map representing
   * the head of the doubly-linked list of the pairs of connected nodes. */
   int node_list[clusters.size()];

  /* Contains various linked lists of arcs (arc = neighbours to a node)*/
   int edge_map[NUM_EDGES][6];

  //Create neighbour map
  if (VERBOSE) cout<<"Creating the linked-list map of neighbours.\n";
  findSegmentNeighbours(edge_map, node_list, reg, heap, heap_handler_list, tree_nodes, max_dist);
  if (VERBOSE) cout<<"Populated the Min-Heap with initial neighbour pairs.\n";

/*
 * //Contents of heap_handler list and node pairs in edge map
  for (int i = 0; i< NUM_EDGES; i++){
    //cout << "heap_handler_list[i_edge]: " << &heap_handler_list[i] << "\n";
    handle_t handle = heap_handler_list[i];
    n_pair linked_edge = *handle;
    cout << " linked_edge: N1 __ N2: " <<  linked_edge.N1 << " " << linked_edge.N2 << "\n";

  }
  for (int i = 0; i< NUM_EDGES; i++){
  cout << "edge_map: N1 __ N2: " <<  edge_map[i][0] << " " << edge_map[i][1] << "\n";
  }*/

  if (VERBOSE) cout<<"Starting the merge algorithm...\n";
  while (!heap.empty()){
    n_pair pnode = heap.top(); //Get min-cost node pair
    heap.pop(); //Erase pnode

    treeNode* merged_node = new treeNode();
    merged_node->isBase = false;
    merged_node->id = cloud_id;

    merged_node->children = new vector<treeNode*>;
    treeNode * n1 = tree_nodes[pnode.N1];
    treeNode * n2 = tree_nodes[pnode.N2];
    merged_node->children->push_back(n1);
    merged_node->children->push_back(n2);
    n1->parent = merged_node;
    n2->parent = merged_node;
    //node->cloud_ref: don't need this since cloud_ref only relevant for visualization
    tree_nodes[pnode.N1] = merged_node; //update ref

    //Take record of new merge operation
    std::pair <treeNode*, treeNode*> m_op;
    m_op = std::make_pair (n1,n2);
    merge_order.push_back(m_op);

    /*
     * //Edge map entries for N1, N2
    cout << "Nodes to merge: " << pnode.N1 << " " << pnode.N2 << "\n";
    if (pnode.N1 != 0 | pnode.N2 != 1){
    int k = node_list[pnode.N1];
    while (k != -1){
      cout << "N1's negihs: " << edge_map[k][0] << " " << edge_map[k][1] << "\n";
      k = edge_map[k][getNodeNextColumn(k, pnode.N1, edge_map)];
    }
    int j = node_list[pnode.N2];
    while (j != -1){
      cout << "N2's negihs: " << edge_map[j][0] << " " << edge_map[j][1] << "\n";
      j = edge_map[j][getNodeNextColumn(j, pnode.N2, edge_map)];
    }
    } else {
      cout << "Hey Oh no! \n";
    }*/

    //Merge edge into one node (N1), remove all obsolete refs to N2.
    mergeEdge (pnode.i_edge, heap_handler_list, edge_map, node_list, heap);

    /* Go through arc of N1, recalculate cost of each adjacent edge and increase/decrease
     its position in the heap respectively. */
    int edge = node_list[pnode.N1]; //head of N1 arc
    int next;
    while (edge != -1){
      //All refs to pnode.i_edge already removed thanks to mergeEdge (above).
      next = getNodeNextColumn(edge, pnode.N1, edge_map);
      handle_t lhandle = heap_handler_list[edge];
      n_pair linked_edge = *lhandle; //Get heap node

      //Update any N2 reference to N1 (N1 now refers to combined cloud)
      if (linked_edge.N1 == pnode.N2){
        linked_edge.N1 = pnode.N1;
        heap.update(heap_handler_list[edge], linked_edge);
      } else if (linked_edge.N2 == pnode.N2) {
        linked_edge.N2 = pnode.N1;
        heap.update(heap_handler_list[edge], linked_edge);
      }

      //Update cost
      int old_cost = linked_edge.cost;
      linked_edge.cost = calcCost(linked_edge.N1, linked_edge.N2, tree_nodes);

      //Update heap position
      if (old_cost < linked_edge.cost){
        heap.increase(lhandle); //O(1)
      } else if (old_cost > linked_edge.cost){
        heap.decrease(lhandle); //O(log(n))
      }

      edge = edge_map[edge][next];
    }
    cloud_id++;
  }


  //Now prepare the merge_op iterator!
  merge_op = merge_order.begin();

  //Do Cost_Function!
  //Continue with visualisation code here (should be two short steps!!)

  if (VERBOSE) cout<<"Finished merge procedure.\n";
  if (VERBOSE) cout<<"Launching visualization for initial regions.\n";

  //If (VISUAL)
  visualizeCloud(clouds, interactionCustomizationVis, keyboardEventOccurred);

  /*for(int i = 0; i < clusters.size(); i++)
        delete tree_nodes[i];
     */


  //delete tree_nodes[69];

  //delete[] tree_nodes;
  //return (0);
  return 0;


}

// Temporary
//or treeNode table, plus N1 and N2 index
float calcCost(int N1, int N2, treeNode** tree_nodes){

  pcl::PointCloud<PointT>::Ptr cloud_n1 (new pcl::PointCloud<PointT>);
  pcl::PointCloud<PointT>::Ptr cloud_n2 (new pcl::PointCloud<PointT>);
  pcl::ExtractIndices<PointT> extract;

  extract.setInputCloud (cloud);
  //Type conversion necessary for matching ExtractIndices::setIndices prototype
  pcl::IndicesPtr indices_a = boost::make_shared<std::vector<int> >(getCloudIndices(tree_nodes[N1]));
  pcl::IndicesPtr indices_b = boost::make_shared<std::vector<int> >(getCloudIndices(tree_nodes[N2]));

  extract.setIndices (indices_a);
  extract.setNegative (false);
  extract.filter (*cloud_n1);

  extract.setIndices (indices_b);
  extract.setNegative (false);
  extract.filter (*cloud_n2);

  return compute_cost (*cloud_n1, *cloud_n2);
}


/* Merges edge (N1---N2) into one node. N1 swallows N2; removes all obsolete refs to N2.
 * 1. For all edges containing a node (N3) linked originally to BOTH N1 & N2,
 * deletes redundant edge with connection to N2.
 * 2. Change all edges inherited from N2 to link to N1 instead.
 * 3. Remove all references to this edge scheduled for merging.
 */
void mergeEdge (int i_edge,                   //edge to merge
                handle_t heap_handler_list[], //Edge indices to heap nodes
                int edge_map[][6],            //Linked lists of arcs
                int node_list[],              //Pointers to linked-list Head of
                                              //connected edges for each node
                fibonacci_heap& heap){        //Heap of Merge pairs

  //1. Intersecting Neighbours of N1, N2 (nodes of i_edge)
  //Find all shared nodes between the two
  int N1 = edge_map[i_edge][0]; //first node
  int N2 = edge_map[i_edge][1]; //second node
  int i_n1 = node_list[N1]; //head of edge list for node_1

  int n2_end = -1; //Last edge linked to N2 (if exists) that is disjoint from N1.
  while (i_n1 != -1) {//For loop: unclear which 'next' to use?
    int n1_nghbr; //Node != N1
    int next; //Iterator for N1

    int x = edge_map[i_n1][1];
    n1_nghbr = (x == N1) ? edge_map[i_n1][0] : x; //get adjacent node
    //get value from n1_nghbr's Next column based on its position (0 or 1).
    next = getNodeNextColumn(i_n1, N1, edge_map);

    if (n1_nghbr != N2) { //Skip merged edge for now.

      int i_n2 = node_list[N2]; //Head of arc list for N2
      n2_end = -1;
      while (i_n2 != -1){

        int x = edge_map[i_n2][1];
        int n2_nghbr = (x == N2) ? edge_map[i_n2][0] : x; //Node adjacent to N2 in i_n2 edge
        int next = getNodeNextColumn(i_n2, N2, edge_map); //Iterator of N2

        if (n2_nghbr != N1) { //Skip merged edge.
          if (n2_nghbr == n1_nghbr){ //Shared node between merged edge endpoints identified.
            //intersecting edge! Remove all refs.
            //cout << "Pair to remove: " << n2_nghbr << " " << N2 << "\n";
            heap.erase(heap_handler_list[i_n2]); //Remove Node pair
            //Cut out the middle man: i_n2 edge must be removed.
            removeNode(n2_nghbr, i_n2, edge_map, node_list);
            removeNode(N2, i_n2, edge_map, node_list);
          } else {
            n2_end = i_n2; // There is at least one either unprocessed OR **non-intersecting node**
          }
        }
        i_n2 = edge_map[i_n2][next];
      }//END while N2
    }
    i_n1 = edge_map[i_n1][next];
  }//END while N1

  removeNode(N2, i_edge, edge_map, node_list);
  if (node_list[N2] != -1){ //There is an edge connected to N2 that was disjoint from N1. Add ref for N1.
    /* N1---N2 Merge!
     * Obliterate merged edge and all refs to N2 (above)
     * 1. Remove N1--N2 from N2's side right off the bat, renew head if needed
     * 2. Scour all the "N2" edges using node list head, change N2 column to N1 instead!
     * 3. Include disjoint N2-chain in N1 ref's. Get preceeding and next edges from N1's
     * side of N1--N2, change their links (prec->) to point to head, (<-next) point to node end.
     */
    int i_n2 = node_list[N2];
    while (i_n2 != -1){ //Update all remaining N2 refs to N1.
      n2_end = i_n2;
      cout << "Pair to update!: " << edge_map[i_n2][0] << " " << edge_map[i_n2][1] << "\n";
      int next = edge_map[i_n2][getNodeNextColumn(i_n2, N2, edge_map)];
      edge_map[i_n2][(edge_map[i_n2][1] == N2) ? 1 : 0] = N1;
      i_n2 = next;
    }
    int n2_head = node_list[N2]; //Head of left over N2 chain

    /* Use additional parameter option to doubly remove i_edge ref in N1 arc,
     * and insert left over N2 arc in its place. */
    //cout << "Start of insert chain: " << edge_map[n2_head][0] << " " << edge_map[n2_head][1] << "\n";
    //cout << "End of insert chain: " << edge_map[n2_end][0] << " " << edge_map[n2_end][1] << "\n";
    removeNode(N1, i_edge, edge_map, node_list, n2_head, n2_end);

  } else {
    removeNode(N1, i_edge, edge_map, node_list);
  }

}

/**
 * Based on Edge-Map Row layout below, identify the NEXT column of the arc
 * containing the node specified in the parameter list.
 * Layout = [N1, N2, NODE_ONE_NEXT, NODE_ONE_PREC, NODE_TWO_NEXT, NODE_TWO_PREC]
 * */
int getNodeNextColumn(int edge, int node, int edge_map[][6]){
  return (edge_map[edge][1] == node) ? NODE_TWO_NEXT : NODE_ONE_NEXT;

}

/**
 * Delete edge from linked-list of arc containing parameter node.
 * OPTIONAL: specify before (link_fwd) and after (link_bck) links
 * to replace those previously referencing edge to-be-removed.
 */
void removeNode (int node, int edge,
                 int edge_map[][6], int node_list[],
                 /* Following 2 params: Used to insert a new edge chain in place of removed edge. */
                 int link_fwd /*DEFAULT: -1*/, int link_bck /*DEFAULT: -1*/){

  int node_next = getNodeNextColumn(edge, node, edge_map);
  int rm_preceeding = edge_map[edge][node_next + 1];
  int rm_next = edge_map[edge][node_next];

  if (link_fwd == -1){
    link_fwd = rm_next;
    link_bck = rm_preceeding;
  } else {
    edge_map[link_fwd][getNodeNextColumn(link_fwd, node, edge_map) + 1] = rm_preceeding;
    edge_map[link_bck][getNodeNextColumn(link_bck, node, edge_map)] = rm_next;
  }

  if (rm_preceeding == -1){ //This is neighbour's head node
    node_list[node] = link_fwd; //change pointer to new head
  } else {
    //Modify Prec, link to next-next node
    edge_map[rm_preceeding][getNodeNextColumn(rm_preceeding, node, edge_map)] = link_fwd;
  }
  if (rm_next != -1)
    edge_map[rm_next][getNodeNextColumn(rm_next, node, edge_map) + 1] = link_bck; //link to prev-prev node
}


/* Calls findRegionsKNN (K-nearest-neighbours) for each segment and saves the results for later use */
void findSegmentNeighbours(int edge_map[][6],            //Linked lists of arcs
                           int node_list[],              //Pointers to linked-list Head of
                                                         //connected edges for each node
                           RegionGrowing& reg,           //Region Growing object
                           fibonacci_heap& heap,         //Heap of Merge pairs
                           handle_t heap_handler_list[], //Edges to heap nodes
                           treeNode** tree_nodes,       //Nodes to tree nodes
                           float max_dist){              //Maximum distance between neighbours
                          //vector<vector<int>>& segment_neighbours,

  int edge_count = 0;
  //int edge_hash_tab[NUM_EDGES][NUM_EDGES] = {0}; //0 is like -1: get() as i - 1, set() as i + 1
  int** edge_hash_tab; /* clusters.size() ^ 2 array providing the row index of the edge for (NodeX, NodeY) */
  edge_hash_tab = new int*[clusters.size()]; //initialize
  for (int i = 0; i < clusters.size(); i++){
    edge_hash_tab[i] = new int[clusters.size()] (); //Brackets enables us to init all values in row to 0.
  }

  //segment_neighbours.resize(clusters.size());
  /*For each segments in clusters array, find and map it to its neighbours. Create a node
   * in the heap for each unique edge.*/
  for (int i_seg = 0; i_seg < clusters.size(); i_seg++){
    //std::vector<int> nghbrs;
    findRegionsKNN(i_seg, reg, edge_map, edge_count, edge_hash_tab,
                   node_list, heap, heap_handler_list, tree_nodes, max_dist);
    //segment_neighbours[i_seg].swap(nghbrs);
  }
  //cout << "Edge count: " << edge_count < "\n";
  // cleanup of edge_hash_table
  for(int i = 0; i < clusters.size(); i++)
    delete[] edge_hash_tab[i];
  delete[] edge_hash_tab;

}

/* Finds the K nearest neighbours of the given segment */
void findRegionsKNN (int index, //segment
                     RegionGrowing& reg,           //Region Growing object
                     int edge_map[][6],            //Linked lists of arcs
                     int& edge_count,              //the number of edges mapped so far
                     int **edge_hash_tab,          //Node pairs to edge indices
                     int node_list[],              //Pointers to linked-list Head of
                                                   //connected edges for each node
                     fibonacci_heap& heap,         //Heap of Merge pairs
                     handle_t heap_handler_list[], //Edges to heap nodes
                     treeNode** tree_nodes,       //Nodes to tree nodes
                     float max_dist){              //Maximum distance between neighbours
                     //std::vector<int>& nghbrs,
  int number_of_points = clusters[index].indices.size(); //num points inside cloud segment
  std::vector<float> distances;
  //float max_dist = std::numeric_limits<float>::max ();
  distances.resize (clusters.size (), max_dist);

  /* FIND ALL neighbours of point, record their distance */
  for (int i_point = 0; i_point < number_of_points; i_point++){
    int point_index = clusters[index].indices[i_point];
    int number_of_neighbours = static_cast<int>(reg.point_neighbours_[point_index].size ());
    for (int i_nghbr = 0; i_nghbr < number_of_neighbours; i_nghbr++){
      // find segment
      int nghbr_point_index = reg.point_neighbours_[point_index][i_nghbr];
      int ngbr_segment_index = -1;
      ngbr_segment_index = reg.point_labels_[nghbr_point_index];

      if ( ngbr_segment_index != index && ngbr_segment_index != -1) //alien (outside) segment
      {
        // try to push it to the queue
        float distance = getDistance(cloud->points[point_index], cloud->points[nghbr_point_index]);
        if (distances[ngbr_segment_index] > distance) //Get closest distance
          distances[ngbr_segment_index] = distance;
      }
    }
  }


  /* KEEP ALL neighbours within MAX_DIST, BUILD doubly-linked list for neighbours of this segment,
   * save the head of neighbour edge list associated with this segment node. */

  int last_accessed_ind = -1; //PRECEEDING
  int lai_next = 0; //Either NODE_ONE_NEXT or NODE_TWO_NEXT of the "last-accessed" edge above.
  int i_edge; //Edge map index of current Edge
  bool head = true; //tracks whether current edge is the head of any node arcs.
  for (int i_segment = 0; i_segment < clusters.size(); i_segment++){
    if (distances[i_segment] < max_dist | (edge_hash_tab[i_segment][index] != 0)) {

      if (edge_hash_tab[index][i_segment] == 0){ //Edge entry is empty; NEW edge
        i_edge = edge_count;
        //+1 since hashtable is indexed starting from 1 (0= -1)
        edge_hash_tab[index][i_segment] = i_edge + 1; //+1 since hashtable is indexed starting from 1 (0 <=> -1)
        edge_hash_tab[i_segment][index] = i_edge + 1;

        edge_map[i_edge][0] = index; //THIS segment
        edge_map[i_edge][1] = i_segment; //Segment neighbour
        edge_map[i_edge][NODE_ONE_PREC] = last_accessed_ind;
        edge_map[i_edge][NODE_ONE_NEXT] = -1;

        if (last_accessed_ind > -1) //If there was a preceeding segment
          edge_map[last_accessed_ind][lai_next] = i_edge;

        last_accessed_ind = i_edge;
        lai_next = NODE_ONE_NEXT;

        //Calculate cost of merge, insert new Node Pair into the Heap.
        n_pair nnode(index, i_segment, i_edge,
                     calcCost(index, i_segment, tree_nodes));
        heap_handler_list[i_edge] = heap.push(nnode); //Save at edge index the node handler (pointer to node).
        edge_count++;
      } else {
        i_edge = edge_hash_tab[index][i_segment] - 1; //Indexed starting from 1. Invalid index = 0.

        edge_map[i_edge][NODE_TWO_PREC] = last_accessed_ind;
        edge_map[i_edge][NODE_TWO_NEXT] = -1;

        if (last_accessed_ind > -1) //If there was a preceeding segment
          edge_map[last_accessed_ind][lai_next] = i_edge;

        last_accessed_ind = i_edge;
        lai_next = NODE_TWO_NEXT;

      }
      if (head){ //Save the head of the linked edge list for future access
        node_list[index] = i_edge;
        head = false;
      }

      //OLD code
      //nghbrs.push_back(i_segment);
    }

  }


}



/* Return total number of edges in segment-neighbours graph. */
int getNumEdges(RegionGrowing& reg, float max_dist){
  int size = 0;
  int edges[clusters.size ()][clusters.size ()] = { 0 };
  for (int c = 0; c < clusters.size(); c++){
    //For each cluster, find neighbouring segments
    for (int i = 0; i < clusters[c].indices.size(); i++){
      int point_index = clusters[c].indices[i];
      for (int j = 0; j < reg.point_neighbours_[point_index].size(); j++){
        int i_nghbr = reg.point_neighbours_[point_index][j];
        int i_seg = reg.point_labels_[i_nghbr];
        if (i_seg == -1 | i_seg == c){
          continue; //Skip non-clustered points
        }
        if (edges[c][i_seg] == 0){
          /* NOTE: KNN nearest neighbours search finds "closest" neighbours for point A
           * at a certain limit N. There is no guarantee of parallel computation
           * such as X is found to be a neighbour of Y, but repeating the same
           * computation for Y, X may not be among the neighbours chosen. */
          float distance = getDistance(cloud->points[point_index], cloud->points[i_nghbr]);
          if (distance < max_dist){
            size++;
            if (i_seg < c){
              //Caught an inconsistent neighbour pairing, add to adjacent point's nghbrs list.
              reg.point_neighbours_[i_nghbr].push_back(point_index);
            }
            edges[c][i_seg] = 1;
            edges[i_seg][c] = 1;

          }

        }

      }
    }
  }
  //if (VERBOSE) cout << "Size: " << size << "\n";
  //assert(size%2 == 0); //Each valid edge should have been double counted.
  return size;
}


float getDistance(PointT& p1, PointT& p2){
  Eigen::Vector3f diff = p1.getVector3fMap() - p2.getVector3fMap();
  return diff.norm ();

}





/** Does a depth first search through the cloud Tree and
 *  gathers all indices of the merged cloud.
 **/
vector<int> getCloudIndices(treeNode * cloud){
  //std::vector <pcl::PointIndices> clusters;
  if (cloud->isBase){
    return clusters[cloud->cluster_index].indices;
  }
  vector<treeNode*>& children = *(cloud->children); // vector is not copied here
  vector<int> i1 = getCloudIndices(children[0]);
  vector<int> i2 = getCloudIndices(children[1]);
  return unionIndices(i1, i2);

}

/* The standard for cloud labelling in the visualization viewer. */
string constructCloudLabel (int id){
  char label[50];
  sprintf(label, "cloud %d%s", id, "\0");
  //cloud v1 %d%s
  string str(label);
  return str;
}

/* Store color in r, g, b variables of point cloud with corresponding label.
 * Created since functionality is missing from PCL libraries. */
void getPointCloudRenderingProperties_Color (
    boost::shared_ptr<pcl::visualization::PCLVisualizer>& viewer,
    string label, double& r, double& g, double &b){
  //cout << "getPointCloudRenderingProperties_Color:  " << &(viewer->getCloudActorMap()->find(label)->second.actor) <<"\n";
  vtkLODActor* actor = vtkLODActor::SafeDownCast( viewer->getCloudActorMap()->find(label)->second.actor);
  actor->GetProperty()->GetColor(r, g, b);
  //cout << "Color: " << std::setprecision (4) << r << " " << g << " " << b << "\n";
}


void keyboardEventOccurred (const pcl::visualization::KeyboardEvent &event,
                            void* viewer_void)
{
  boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer = *static_cast<boost::shared_ptr<pcl::visualization::PCLVisualizer> *> (viewer_void);

  //SPLIT OPERATION
  if (event.getKeySym () == "a" && event.keyDown ())
  {
    std::cout << "A (Split) was pressed" << std::endl;

    if (merge_op == merge_order.begin())
      return; //exit if splitting is at base (beginning)
    merge_op--;
    pair<treeNode*, treeNode*> mp = *merge_op; //pair of split clouds
    pcl::PointCloud<PointT>::Ptr cloud_a (new pcl::PointCloud<PointT>);
    pcl::PointCloud<PointT>::Ptr cloud_b (new pcl::PointCloud<PointT>);

    pcl::ExtractIndices<PointT> extract;
    extract.setInputCloud (cloud);
    //Type conversion necessary for matching ExtractIndices::setIndices prototype
    pcl::IndicesPtr indices_a = boost::make_shared<std::vector<int> >(getCloudIndices(mp.first));
    pcl::IndicesPtr indices_b = boost::make_shared<std::vector<int> >(getCloudIndices(mp.second));

    extract.setIndices (indices_a);
    extract.setNegative (false);
    extract.filter (*cloud_a);

    extract.setIndices (indices_b);
    extract.setNegative (false);
    extract.filter (*cloud_b);

    cloud_list.erase(mp.first->parent->cloud_ref);
    cloud_list.push_back(cloud_a);
    mp.first->cloud_ref = --cloud_list.end(); //last element
    cloud_list.push_back(cloud_b);
    mp.second->cloud_ref = --cloud_list.end(); //last element

    string old_label = constructCloudLabel(mp.first->parent->id);
    /*cout << "Split, first cloud: " << (*(mp.first->cloud_ref))->points.size() << "\n";
    cout << "Split, second cloud: " << (*(mp.second->cloud_ref))->points.size() << "\n";
    cout << "Split IDs: " << mp.first->id << " " << mp.second->id << "\n";*/

    double old_r, old_g, old_b;
    getPointCloudRenderingProperties_Color(viewer, old_label, old_r, old_g, old_b);
    viewer->removePointCloud(old_label);

    string label_a = constructCloudLabel(mp.first->id);
    string label_b = constructCloudLabel(mp.second->id);
    viewer->addPointCloud(cloud_a, label_a);
    //reuse color for consistency
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR,
                                             old_r, old_g, old_b /*R,G,B*/, label_a);
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, label_a);

    viewer->addPointCloud(cloud_b, label_b);
    //Random colour for one of the clouds from split
    //A value of 1.0 is equivalent to 255, a value of 0.0 to 0.
    double r = (rand() % 256) / 255.0;
    double g = (rand() % 256) / 255.0;
    double b = (rand() % 256) / 255.0;
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR,
                                                 r, g, b /*R,G,B*/, label_b);
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, label_b);

  }
  //MERGE EVENT
  else if (event.getKeySym () == "d" && event.keyDown ()){
    std::cout << "D (Merge) was pressed" << std::endl;
    if (merge_op == merge_order.end())
      return; //exit if merge is at maximum

    pair<treeNode*, treeNode*> mp = *merge_op;

    pcl::PointCloud<PointT>::Ptr m_cloud (new pcl::PointCloud<PointT>);

    *m_cloud = **(mp.first->cloud_ref);
    *m_cloud += **(mp.second->cloud_ref); //concatenate clouds

    /* cout << "Merge, first cloud: " << (*(mp.first->cloud_ref))->points.size() << "\n";
    cout << "Merge, second cloud: " << (*(mp.second->cloud_ref))->points.size() << "\n";
    cout << "Merge, IDS: " << mp.first->id << " " << mp.second->id << "\n";
    cout << "Merge, parent: " << mp.first->parent->id << "\n";
    cout << "Merged cloud size: " << " " << m_cloud->points.size() <<"\n"; */
    //hopefully, the above worked

    cloud_list.erase(mp.first->cloud_ref);
    cloud_list.erase(mp.second->cloud_ref);
    cloud_list.push_back(m_cloud);

    //Parent tree node cloud creation ref
    mp.first->parent->cloud_ref = --cloud_list.end(); //last element

    merge_op++; //increment to next merge op

    string label_a = constructCloudLabel(mp.first->id);
    string label_b = constructCloudLabel(mp.second->id);
    double old_r, old_g, old_b;
    //reuse a color from the one of the two clouds
    getPointCloudRenderingProperties_Color(viewer, label_a, old_r, old_g, old_b);
    viewer->removePointCloud(label_a);
    viewer->removePointCloud(label_b);
    string new_label = constructCloudLabel(mp.first->parent->id);
    viewer->addPointCloud(m_cloud, new_label);
    //reuse color for consistency
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR,
                                            old_r, old_g, old_b /*R,G,B*/, new_label);
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 8, new_label);
  }

  //SLOWS everything down:
  //viewer->spinOnce (100); //update the screen

}

