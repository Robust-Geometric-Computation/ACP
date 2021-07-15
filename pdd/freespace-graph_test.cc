#include "freespace.h"
#include "freespace-graph.h"

//#define TWO_ROOM_EXAMPLE
//#define LATTICE_EXAMPLE

void testCreateGraph(const char * blockspace_dir, const char * graph_dir, double theta, double clearance_unit, int num_levels) {
  int numRotations = floor(2*M_PI/theta);
  
  std::vector<PTRPolyhedron> close_blockspaces;
  char s[50];
  for (int i=0; i<=40; i++) {
    sprintf(s, "%s/close%02d-out.vtk", blockspace_dir, i);
    PTRPolyhedron p = loadPolyVTK(s);
    close_blockspaces.push_back(p);
  }

  std::vector<PTRPolyhedron> rough_blockspaces;
  for (int i=0; i<=40; i++) {
    sprintf(s, "%s/rough%02d-out.vtk", blockspace_dir, i);
    PTRPolyhedron p = loadPolyVTK(s);
    rough_blockspaces.push_back(p);
  }

  FreeSpaceGraph graph(close_blockspaces, rough_blockspaces, theta, clearance_unit, num_levels, graph_dir);
}

void testSearchGraph(const char * graph_dir, PTR<Point> start, PTR<Point> end, int startRotationIndex, int endRotationIndex) {
  FreeSpaceGraph test(graph_dir);
  std::vector<std::pair<PTR<Point>, double > > path;

  PTRPolyhedron p = loadPoly((std::string(graph_dir) + "/0-" + std::to_string(startRotationIndex) + ".tri").c_str());
  assert(!p->contains(start));
  int cell_index = p->containingCell(start);
  PTRFreeSpaceGraphNode startNode = test.graph[0][startRotationIndex]->get(cell_index);
  /*delete*/ p;

  p = loadPoly((std::string(graph_dir) + "/0-" + std::to_string(endRotationIndex) + ".tri").c_str());
  assert(!p->contains(end));
  cell_index = p->containingCell(end);
  PTRFreeSpaceGraphNode endNode = test.graph[0][endRotationIndex]->get(cell_index);
  /*delete*/ p;

  test.getPath(startNode, endNode, start, end, path);
  cout<<"done search graph"<<endl;

  cout<<endl<<endl;
  for (int i=0; i<path.size(); i++)
    cout<<path[i].first->getApprox().getX().mid()<<" "<<path[i].first->getApprox().getY().mid()<<" "<<path[i].first->getApprox().getZ().mid()<<" "<<path[i].second<<endl;
}

   
#ifdef TWO_ROOM_EXAMPLE
const char *blockspace_dir = "frustum-tworooms-output";
#endif

#ifdef LATTICE_EXAMPLE
const char *blockspace_dir = "skinny-lattice-output";
#endif		
																		    
int main (int argc, char *argv[]) {
  //time_t start_t,end_t;
  //time(&start_t);
  double t0 = getTime();
  //srand(7);
  int numRotations = 40;
  PTR<Object<Scalar>> tan_half_angle = new InputParameter(tanHalfAngle(numRotations));
  
  acp::enable();

  if (argc < 2) {
    cout<<"not enough arguments"<<endl;
    return 1;
  }

  char * filename =  argv[1];

  if (argc == 4) { 
    unsigned seed = atoi(argv[3]);
    cout << "seed " << seed << endl;
    srandom(seed);
  }

  PTR<Polyhedron> poly = loadPoly_VTK(filename, true);
  if (poly == NULL)
    return 1;

  PTR<Polyhedron> obstacle;
  if (argc >= 3) 
    obstacle = loadPoly_VTK(argv[2], true); 
  else {
    double bounds[6] = { 5, 7, 5, 7, 5, 7};
    obstacle = box(bounds);
  }

  FreeSpace * fs  = new FreeSpace(poly, obstacle, tan_half_angle, numRotations);

  acp::enable();

  if (argc == 2) { 
    unsigned seed = atoi(argv[1]);
    srandom(seed);
  }

#ifdef TWO_ROOM_EXAMPLE
   // --- Two Room Example ---
   const char * graph_dir = "frustum-tworooms-graph";
   blockspace_dir = "frustum-tworooms-output";
   double theta = M_PI / 20;
   double clearance_unit = 0.05;
   int num_levels = 5;
  
   testCreateGraph(blockspace_dir, graph_dir, theta, clearance_unit, num_levels);
   double t1 = getTime();
   cerr << "elapsed time: " << t1 - t0 << endl;
   //time(&end_t);
   //cout << "Elapsed time: "<< difftime (end_t,start_t)/60.0 << " minutes" << endl;

   acp::primitiveReport();
   acp::resetReport();

   PTR<Point> start = new Point(2.5,-1.0, 0.0);
   PTR<Point> end = new Point(2.5,8.0,0.0);
   int startRotationIndex = 0;
   int endRotationIndex = 0;
  
   //time(&start_t);
   testSearchGraph(graph_dir, start, end, startRotationIndex, endRotationIndex);
   //time(&end_t);
   double t2 = getTime();
   cout << "elapsed time: " << t2 - t1 << endl;
   //cout << "Elapsed time: "<< difftime (end_t,start_t) << " seconds" << endl;
   //--------------------------------
#endif




#ifdef LATTICE_EXAMPLE
  // --- Lattice Example ---
  const char * graph_dir = "skinny-lattice-graph";
  blockspace_dir = "skinny-lattice-output";
  double theta = M_PI / 20;
  double clearance_unit = 0.05;
  int num_levels = 1;
  
  //time(&start_t);

  testCreateGraph(blockspace_dir, graph_dir, theta, clearance_unit, num_levels);
  double t1 = getTime();
  cout << "elapsed time: "<< t1 - t0 << endl;
  //cout << "Elapsed time: "<< difftime (end_t,start_t)/60.0 << " minutes" << endl;

  acp::primitiveReport();
  acp::resetReport();

  PTR<Point> start = new Point(-2.0,0.0, -0.5);
  PTR<Point> end = new Point(5.0,0.0,0.0);
  int startRotationIndex = 0;
  int endRotationIndex = 0;
  
  //time(&start_t);
  testSearchGraph(graph_dir, start, end, startRotationIndex, endRotationIndex);
  double t2 = getTime();
  cout << "elapsed time: " << t2 - t1 << endl;
  //cout << "Elapsed time: "<< difftime (end_t,start_t) << " seconds" << endl;
  //--------------------------------
#endif

  acp::primitiveReport();
}
