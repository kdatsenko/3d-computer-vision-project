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
bool VERBOSE = false;
bool VISUAL = false;
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
const int FIRST_NODE_NEXT = 2;

//EDGE TABLE doubly-linked pointer indices
const int NODE_ONE_NEXT = 2;
const int NODE_ONE_PREC = 4;
const int NODE_TWO_NEXT = 3;
const int NODE_TWO_PREC = 5;


typedef boost::heap::fibonacci_heap<n_pair, boost::heap::compare<compare_node> > fibonacci_heap;
typedef typename boost::heap::fibonacci_heap<n_pair, boost::heap::compare<compare_node> >::handle_type handle_t;


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
                           treeNode* tree_nodes[],
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
                     treeNode* tree_nodes[],
                     float max_dist);
vector<int>* getCloudIndices(treeNode * cloud);
string constructCloudLabel (int id);
void getPointCloudRenderingProperties_Color (
    boost::shared_ptr<pcl::visualization::PCLVisualizer>& viewer,
    string label, double& r, double& g, double &b);
void keyboardEventOccurred (const pcl::visualization::KeyboardEvent &event,
                            void* viewer_void);
int getNumEdges(RegionGrowing& reg, float max_dist);
float Cost_Function(int N1, int N2, treeNode* tree_nodes[]);


//Return a set of outliers given a cloud and a set of inliers
vector<int>* unionIndices(vector<int>& U_inliers, vector<int>& inliers){
	vector<int> * result= new vector<int>; //heap
	int index=0;
	int size = inliers.size();
	int i;
	for (i=0; i < U_inliers.size() | (index == size); i++){
		if (U_inliers[i] < inliers[index]){
			result->push_back(U_inliers[i]);
		} else {
			result->push_back(inliers[index]);
			index++;
		}
	}
	if (U_inliers.size() != i){
		for (; index != size; index++)
			result->push_back(inliers[index]);
	} else {
		for (; i != U_inliers.size(); i++)
			result->push_back(U_inliers[i]);
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
  std::vector<std::vector<float> > models;

  //Initialize new point cloud, save pointer.
  cloud = pcl::PointCloud<pcl::PointXYZ>::Ptr (new pcl::PointCloud<pcl::PointXYZ>);
  cloud.get()->points.clear();
  pcl::PointCloud<PointT>::Ptr cloud_filtered (new pcl::PointCloud<PointT>);

  if (pcl::io::loadPCDFile<PointT> (inputFile, *cloud) == -1) //* load the file
  {
    PCL_ERROR ("Couldn't read .pcd file\n");
    return (-1);
  }

  if (VERBOSE) cout<<"Loaded point cloud...\n";
  //used to identify segmented point cloud files
  char cloud_name[strlen(inputFile) - 4 + 1]; //remove .pcd
  strncpy(cloud_name, inputFile, strlen(inputFile) - 4);

  pcl::search::Search<PointT>::Ptr tree = boost::shared_ptr<pcl::search::Search<PointT> > (new pcl::search::KdTree<PointT>);
  pcl::PointCloud <pcl::Normal>::Ptr cloud_normals (new pcl::PointCloud <pcl::Normal>);
  pcl::PCDWriter writer;

  std::cout << "PointCloud before filtering has: " << cloud->points.size ()  << " data points." << std::endl;

  // Filter point cloud for spurious NaNs
  // Create the filtering object
  pcl::StatisticalOutlierRemoval<pcl::PointXYZ> sor;
  sor.setInputCloud (cloud);
  sor.setMeanK (50); //number of neighbours
  sor.setStddevMulThresh (5.0);
  sor.filter (*cloud_filtered);

  pcl::NormalEstimation<PointT, pcl::Normal> normal_estimator;
  normal_estimator.setSearchMethod (tree);
  normal_estimator.setInputCloud (cloud_filtered);
  normal_estimator.setKSearch (50);
  normal_estimator.compute (*cloud_normals);



  std::cout << "PointCloud after filtering has: " << cloud_filtered->points.size ()  << " data points." << std::endl;

  RegionGrowing reg;
  reg.setMinClusterSize (50);
  reg.setMaxClusterSize (100000);
  reg.setSearchMethod (tree);
  reg.setNumberOfNeighbours (30);
  reg.setInputCloud (cloud_filtered);
  //reg.setIndices (indices);
  reg.setInputNormals (cloud_normals);
  reg.setSmoothnessThreshold (2.0 / 180.0 * M_PI);
  reg.setCurvatureThreshold (1.0);

  reg.extract (clusters);

  int x = 0;
  std::cout << "Number of clusters is equal to " << clusters.size () << std::endl;
  for (int i = 0; i < clusters.size(); i++){
    x += clusters[i].indices.size();


    /*pcl::PointCloud<PointT>::Ptr cloud_cluster (new pcl::PointCloud<PointT>);
	  int counter = 0;
	  while (counter < clusters[i].indices.size ()){
	    cloud_cluster->points.push_back (cloud_filtered->points[clusters[i].indices[counter]]);
	    counter++;
	  }
	  cloud_cluster->width = cloud_cluster->points.size ();
	  cloud_cluster->height = 1;
	  cloud_cluster->is_dense = true;
	  std::cout << "PointCloud representing the Cluster: " << cloud_cluster->points.size () << " data points." << std::endl;
	  std::stringstream ss;
	  ss << "clusters/" << cluster_folder << "/cloud_cluster_" << i << ".pcd";
	  writer.write<PointT> (ss.str (), *cloud_cluster, false); //* */
  }

  std::cout << "Size: " << x << std::endl;

  // Create the filtering object
  pcl::ExtractIndices<PointT> extract;
  vector<pcl::PointCloud<PointT>::Ptr> clouds;

  int cloud_id = clusters.size(); //initialize identifier for cloud label
  treeNode* tree_nodes[clusters.size()]; //NUMber of nodes (n1, n2, n3, etc)

  //Cloud creation
  for (int i = 0; i < clusters.size(); i++){
    // Extract the inliers
    // copies all inliers of the model computed to another PointCloud
    pcl::PointCloud<PointT>::Ptr cloud_temp (new pcl::PointCloud<PointT>);
    extract.setInputCloud (cloud_filtered);

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
    node->cloud_ref;
    node->cloud_ref = cloud_list.end();

    //no children, no parent yet!
    tree_nodes[i] = node;
  }

  //std::make_pair(-1, -1)

  // put merge code HERE
  float max_dist = 5;
  int NUM_EDGES = getNumEdges(reg, max_dist);

  fibonacci_heap heap;
  handle_t heap_handler_list[NUM_EDGES];


  int node_list[clusters.size()]; //Links to rows in the Edge Map

  int edge_map[NUM_EDGES][6];

  //Create neighbour map
  findSegmentNeighbours(edge_map, node_list, reg, heap, heap_handler_list, tree_nodes, max_dist);

  while (!heap.empty()){
    n_pair pnode = heap.top();

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
    //node->cloud_ref = &cloud_list.back(); don't need this since cloud_ref hasn't really been touched yet.
    tree_nodes[pnode.N1] = merged_node;

    //Take record of new merge operation
    vector< pair<treeNode*, treeNode*> > merge_order;
    std::pair <treeNode*, treeNode*> m_op;
    m_op = std::make_pair (n1,n2);
    merge_order.push_back(m_op);

    //Merge edge into one node, remove all obsolete refs to N2.
    mergeEdge (pnode.i_edge, heap_handler_list, edge_map, node_list, heap);

    //Go through all links, recalc cost! & increcr/decrease
    int edge = node_list[pnode.N1];
    int next;
    while (edge != -1){
      //won't ever come across i_edge, don't worry
      next = getNodeNextColumn(edge, pnode.N1, edge_map);
      handle_t lhandle = heap_handler_list[edge];
      n_pair linked_edge = *lhandle;

      //Remove Update any N2 reference to N1 (N1 now refers to combined cloud)
      if (linked_edge.N1 == pnode.N2)
        linked_edge.N1 = pnode.N1;
      else if (linked_edge.N2 == pnode.N2)
        linked_edge.N2 = pnode.N1;

      //Update cost
      int old_cost = linked_edge.cost;
      linked_edge.cost = Cost_Function(linked_edge.N1, linked_edge.N2, tree_nodes);

      if (old_cost < linked_edge.cost)
        heap.increase(lhandle);
      else if (old_cost > linked_edge.cost)
        heap.decrease(lhandle);

      edge = edge_map[edge][next];
    }

    //don't forget to destroy this pnode!
    //pnode = NULL;
    heap.pop(); //Erase pnode

    cloud_id++;
  }

  //Now prepare the merge_op iterator!
  merge_op = merge_order.begin();

  //Do Cost_Function!
  //Continue with visualization code here (shoud be two short steps!!)

  //std::vector <pcl::PointIndices> clusters; Global now







  return (0);
}

//treenode1, treeNode2,
//or treeNode table, plus N1 and N2 index
float Cost_Function(int N1, int N2, treeNode* tree_nodes[]){

  pcl::PointCloud<PointT>::Ptr cloud_n1 (new pcl::PointCloud<PointT>);
  pcl::PointCloud<PointT>::Ptr cloud_n2 (new pcl::PointCloud<PointT>);
  pcl::ExtractIndices<PointT> extract;
  extract.setInputCloud (cloud);
  //Type conversion necessary for matching ExtractIndices::setIndices prototype
  pcl::IndicesPtr indices_a = boost::make_shared<std::vector<int> >(*getCloudIndices(tree_nodes[N1]));
  pcl::IndicesPtr indices_b = boost::make_shared<std::vector<int> >(*getCloudIndices(tree_nodes[N2]));

  extract.setIndices (indices_a);
  extract.setNegative (false);
  extract.filter (*cloud_n1);

  extract.setIndices (indices_b);
  extract.setNegative (false);
  extract.filter (*cloud_n2);

  return compute (*cloud_n1, *cloud_n2);
}



int getNumEdges(RegionGrowing& reg, float max_dist){
  int size = 0;
  for (int i = 0; i < reg.point_neighbours_.size(); i++){
    int segments[clusters.size ()] = {0};
    for (int j = 0; j < reg.point_neighbours_[i].size(); j++){
      int i_seg = reg.point_labels_[reg.point_neighbours_[i][j]];
      if (segments[i_seg] != 0){
        float distance = pcl::geometry::distance(cloud->points[i], cloud->points[reg.point_neighbours_[i][j]]);
        if (distance < max_dist){
          size++;
          segments[i_seg] = 1;
        }
      }
    }
  }
  assert(size%2 == 0); //Each valid edge should have been double counted.
  return size / 2;
}


//class PCL_EXPORTS RegionGrowing : public pcl::PCLBase<PointT>

void mergeEdge (//HEAP,
                  int i_edge,
                  handle_t heap_handler_list[],
                  int edge_map[][6],
                  int node_list[],
                  fibonacci_heap& heap){
  //1. Intersecting Neighnours between N1, N2 of i_edge
  //Find all shared nodes between the two
  int N1 = edge_map[i_edge][0]; //first node
  int N2 = edge_map[i_edge][1]; //second node
  int i_n1 = node_list[N1]; //head of edge list for node_1

  int n2_end = -1;
  while (i_n1 != -1) {//For: which next to use?
    int n1_nghbr;
    int next; //Iterator for N1

    int x = edge_map[i_n1][1];
    n1_nghbr = (x == N1) ? edge_map[i_n1][0] : x;
    next = getNodeNextColumn(i_n1, N1, edge_map);

    if (n1_nghbr != N2) { //Skip merged edge for now.

      int i_n2 = node_list[N2];
      n2_end = -1;
      while (i_n2 != -1){

        int x = edge_map[i_n2][1];
        int n2_nghbr = (x == N2) ? edge_map[i_n2][0] : x;
        next = getNodeNextColumn(i_n2, N2, edge_map); //Iterator of N2

        if (n2_nghbr != N1) { //Otherwise, skip this edge.
          if (n2_nghbr == n1_nghbr){ //Shared node between merged edge endpoints
            //intersecting edge!

            //Insert Remove NODE Nerd code [HERE]
            heap.erase(heap_handler_list[i_n2]);
            //heap.er
            ///(heap_handler_list[i_n2]);
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

  if (n2_end != -1){ //There is an edge connected to N2 that was disjoint from N1. Add ref for N1.

    //N1---------N2 Merge!
    //Obliterate merged edge and all refs to N2
    /*
     * 1. Remove N1--N2 from N2's side right off the bat, renew head if needed
     * 2. Scour all the "N2" edges using node list head, change N2 column to N1 instead!
     * 3. Include disjoint N2-chain in N1 ref's. Get preceeding and next edges from N1's
     * side of N1--N2, change their links (prec) to point to head, (2) points to node end.
     */
    removeNode(N2, i_edge, edge_map, node_list);

    int i_n2 = node_list[N2];
    while (i_n2 != -1){
      int next = edge_map[i_n2][getNodeNextColumn(i_n2, N2, edge_map)];
      edge_map[i_n2][(edge_map[i_n2][1] == N2) ? 1 : 0] = N1;
      i_n2 = next;
    }
  }

  int n2_head = node_list[N2];
  removeNode(N1, i_edge, edge_map, node_list, n2_head, n2_end);

}


int getNodeNextColumn(int edge, int node, int edge_map[][6]){
  return (edge_map[edge][1] == node) ? NODE_TWO_NEXT : NODE_ONE_NEXT;

}

void removeNode (int node, int edge,
                 int edge_map[][6], int node_list[],
                 int link_fwd, int link_bck /* For inserting a new edge chain in place of one to-be-removed. */){

  int node_next = getNodeNextColumn(edge, node, edge_map);
  int rm_preceeding = edge_map[edge][node_next + 1];
  int rm_next = edge_map[edge][node_next];

  if (link_fwd == -1){
    link_fwd = rm_next;
    link_bck = rm_preceeding;
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
void findSegmentNeighbours(int edge_map[][6],
                           int node_list[],
                           //vector<vector<int>>& segment_neighbours,
                           RegionGrowing& reg,
                           fibonacci_heap& heap,
                           handle_t heap_handler_list[],
                           treeNode* tree_nodes[],
                           float max_dist){

  int edge_count = 0;
  //int edge_hash_tab[NUM_EDGES][NUM_EDGES] = {0}; //0 is like -1: get() as i - 1, set() as i + 1
  int** edge_hash_tab;
  edge_hash_tab = new int*[NUM_EDGES];
  for (int i = 0; i < NUM_EDGES; i++){
    edge_hash_tab[i] = new int[NUM_EDGES]();
  }

  /*double** buf2d;
  buf2d = new double*[number_of_parts];
  for(i = 0; i < number_of_parts; ++i)
    buf2d[i] = new double[12]*/

  int num_seg = clusters.size();
  //segment_neighbours.resize(num_seg);
  for (int i_seg = 0; i_seg < num_seg; i_seg++){
    //std::vector<int> nghbrs;
    //hello(heap_handler_list, heap);
    findRegionsKNN(i_seg, reg, edge_map, edge_count, edge_hash_tab,
                   node_list, heap, heap_handler_list, tree_nodes, max_dist);
    //segment_neighbours[i_seg].swap(nghbrs);
  }
  // cleanup of edge_hash_table
  for(int i = 0; i < NUM_EDGES; ++i)
    delete[] edge_hash_tab[i];
  delete[] edge_hash_tab;

}

/* Finds the K nearest neighbours of the given segment */
void findRegionsKNN (int index,
                     RegionGrowing& reg,
                     int edge_map[][6],
                     int& edge_count,
                     int **edge_hash_tab,
                     int node_list[],
                     //std::vector<int>& nghbrs,
                     fibonacci_heap& heap,
                     handle_t heap_handler_list[],
                     treeNode* tree_nodes[],
                     float max_dist){
  int number_of_points = clusters[index].indices.size();
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

      if ( ngbr_segment_index != index ) //alien (outside) segment
      {
        // try to push it to the queue
        float distance = pcl::geometry::distance(cloud->points[point_index], cloud->points[nghbr_point_index]);
        if (distances[ngbr_segment_index] > distance)
          distances[ngbr_segment_index] = distance;
      }
    }
  }

  /* KEEP ALL neighbours within MAX_DIST, BUILD doubly-linked list for neighbours of this segment,
   * save the head of neighbour edge list associated with this segment node. */

  int last_accessed_ind = -1; //PRECEEDING
  int lai_next = 0; //Either NODE_ONE_NEXT or NODE_TWO_NEXT of the LA edge above.
  int i_edge; //Edge map index of current Edge
  bool head = true;
  for (int i_segment = 0; i_segment < clusters.size(); i_segment++){
    if (distances[i_segment] < max_dist) {

      if (edge_hash_tab[index][i_segment] == 0){ //Edge entry is empty; NEW edge
        i_edge = edge_count;
        //+1 since hashtable is indexed starting from 1 (0= -1)
        edge_hash_tab[index][i_segment] = i_edge + 1; //+1 since hashtable is indexed starting from 1 (0 <=> -1)

        edge_map[i_edge][0] = index; //THIS segment
        edge_map[i_edge][1] = i_segment; //Segment neighbour
        edge_map[i_edge][NODE_ONE_PREC] = last_accessed_ind;
        edge_map[i_edge][NODE_ONE_NEXT] = -1;

        if (last_accessed_ind > -1) //If there was a preceeding segment
          edge_map[last_accessed_ind][lai_next] = i_edge;

        last_accessed_ind = i_edge;
        lai_next = NODE_ONE_NEXT;

        //Insert new Node Pair into the Heap.
        n_pair nnode(index, i_segment, i_edge,
                     Cost_Function(index, i_segment, tree_nodes));
        heap_handler_list[i_edge] = heap.push(nnode); //Save at edge index the node handler (pointer to node).
        //calc Cost for this node?? [HERE]

        edge_count++;
      } else {
        i_edge = edge_hash_tab[index][i_segment] - 1;

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


/** Does a depth first search through the cloud Tree and
 * gathers all indices of the merged cloud.
 **/
vector<int>* getCloudIndices(treeNode * cloud){
  //std::vector <pcl::PointIndices> clusters;
  if (cloud->isBase){
    return &clusters[cloud->cluster_index].indices;
  }
  vector<treeNode*>& children = *(cloud->children); // vector is not copied here
  vector<int>* i1 = getCloudIndices(children[0]);
  vector<int>* i2 = getCloudIndices(children[1]);
  return unionIndices(*i1, *i2);

}

/* The standard for cloud label */
string constructCloudLabel (int id){
  string label = "cloud " + id;
  return label;
}

/* Store color in r, g, b variables of point cloud with corresponding label.
 * Created since functionality is missing from PCL libraries. */
void getPointCloudRenderingProperties_Color (
    boost::shared_ptr<pcl::visualization::PCLVisualizer>& viewer,
    string label, double& r, double& g, double &b){
  vtkLODActor* actor = vtkLODActor::SafeDownCast( viewer->getCloudActorMap()->find(label)->second.actor );
  actor->GetProperty()->GetColor(r, g, b);
}

void keyboardEventOccurred (const pcl::visualization::KeyboardEvent &event,
                            void* viewer_void)
{
  boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer = *static_cast<boost::shared_ptr<pcl::visualization::PCLVisualizer> *> (viewer_void);

  //SPLIT OPERATION
  if (event.getKeySym () == "a" && event.keyDown ())
  {
    std::cout << "A was pressed" << std::endl;

    if (merge_op == merge_order.begin())
      return; //exit if splitting is at base (beginning)
    merge_op--;
    pair<treeNode*, treeNode*> mp = *merge_op; //pair of split clouds
    pcl::PointCloud<PointT>::Ptr cloud_a (new pcl::PointCloud<PointT>);
    pcl::PointCloud<PointT>::Ptr cloud_b (new pcl::PointCloud<PointT>);

    pcl::ExtractIndices<PointT> extract;
    extract.setInputCloud (cloud);
    //Type conversion necessary for matching ExtractIndices::setIndices prototype
    pcl::IndicesPtr indices_a = boost::make_shared<std::vector<int> >(*getCloudIndices(mp.first));
    pcl::IndicesPtr indices_b = boost::make_shared<std::vector<int> >(*getCloudIndices(mp.second));

    extract.setIndices (indices_a);
    extract.setNegative (false);
    extract.filter (*cloud_a);

    extract.setIndices (indices_b);
    extract.setNegative (false);
    extract.filter (*cloud_b);

    cloud_list.erase(mp.first->parent->cloud_ref);
    cloud_list.push_back(cloud_a);
    cloud_list.push_back(cloud_b);
    //Already invalidated, but just in case
    //mp.first->parent->cloud_ref = NULL;
    //Parent tree node cloud creation ref
    mp.first->parent->cloud_ref = cloud_list.end();
    mp.first->cloud_ref = cloud_list.end()--; //second last
    mp.second->cloud_ref = cloud_list.end();

    string old_label = constructCloudLabel(mp.first->parent->id);

    double old_r, old_g, old_b;
    getPointCloudRenderingProperties_Color(viewer, old_label, old_r, old_g, old_b);
    viewer->removePointCloud(old_label);

    string label_a = constructCloudLabel(mp.first->id);
    string label_b = constructCloudLabel(mp.second->id);
    viewer->addPointCloud(cloud_a, label_a);
    //reuse color for consistency
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR,
                                             old_r, old_g, old_b /*R,G,B*/, label_a);
    viewer->addPointCloud(cloud_b, label_b);

    //Random colour for one of the clouds from split
    //A value of 1.0 is equivalent to 255, a value of 0.0 to 0.
    double r = (rand() % 256) / 255.0;
    double g = (rand() % 256) / 255.0;
    double b = (rand() % 256) / 255.0;
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR,
                                                 r, g, b /*R,G,B*/, label_a);

  }
  //MERGE EVENT
  else if (event.getKeySym () == "d" && event.keyDown ()){
    std::cout << "D was pressed" << std::endl;
    if (merge_op == merge_order.end())
      return; //exit if merge is at maximum
    pair<treeNode*, treeNode*> mp = *merge_op;
    pcl::PointCloud<PointT> m_cloud;


    //mp.first->cloud_ref
   // pcl::PointCloud<PointT>::Ptr cloud_a;
    //cloud_a->
    //pcl::PointCloud<PointT>::Ptr cloud_a (new pcl::PointCloud<PointT>);
    m_cloud = **(mp.first->cloud_ref);
    m_cloud += **(mp.second->cloud_ref); //concatenate clouds
    pcl::PointCloud<PointT>::Ptr cloud_c (&m_cloud);

    //hopefully, the above worked
    cloud_list.erase(mp.first->cloud_ref);
    cloud_list.erase(mp.second->cloud_ref);

    cloud_list.push_back(cloud_c);

    //Already invalidated, but just to be accurate
    /*
    mp.first->cloud_ref = NULL;
    mp.second->cloud_ref = NULL;*/
    //Parent tree node cloud creation ref
    mp.first->parent->cloud_ref = cloud_list.end();

    merge_op++; //increment to next merge op


    string label_a = constructCloudLabel(mp.first->id);
    string label_b = constructCloudLabel(mp.second->id);
    double old_r, old_g, old_b;
    //reuse a color from the one of the two clouds
    getPointCloudRenderingProperties_Color(viewer, label_a, old_r, old_g, old_b);
    viewer->removePointCloud(label_a);
    viewer->removePointCloud(label_b);

    string new_label = constructCloudLabel(mp.first->parent->id);
    viewer->addPointCloud(cloud_c, new_label);
    //reuse color for consistency
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR,
                                             old_r, old_g, old_b /*R,G,B*/, label_a);
  }

  //MOVE this????
  viewer->spinOnce (100); //update the screen

}

