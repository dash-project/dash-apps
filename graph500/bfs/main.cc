#include "../generator/make_graph.h"
#include "../generator/utils.h"
#include <libdash.h>
#include <list>

packed_edge * create_graph_tuples() {
  int scale = 26;
  int edge_factor = 16;
  uint64_t seed1 = 2, 
    seed2 = 3;
  make_mrg_seed(seed1, seed2, seed);
  
  //number of edges = edge_factor * 2^scale
  packed_edge[local_edge__count] edges;
  generate_kronecker_range(seed, scale, start_edge, start_edge + local_edge_count, edges);
  return edges;
}

uint64_t * find_root_vertices(packed_edge & edges, int num_root_vertices) {
	std::array<uint64_t, num_root_vertices> roots;
	// to be implemented
	return roots;
}

int main() {
  dash::init();
  dart_unit_t myid   = dash::myid();
  size_t size        = dash::Team::All().size();
  int num_root_vertices = 64;

  // requires teamsize = 2*(n^2) for some n
  // for the beginning, 2^n would be easier
  // dash::1DGraphPartitioning: each processor holds the same amount of vertices
  // dash::2DGraphPartitioning: each processor holds roughly the same amount of edges
  dash::2DGraphPartitioning partition(Team::All())
  dash::UndirectedGraph<uint64_t> graph(partition);
  
  // create graph tuples local to each node
  uint64_t global_edge_count = edge_factor << scale;
  uint64_t local_edge_count = global_edge_count / size;
  uint64_t start_edge = global_edge_count * myid;
  packed_edge * edges = create_graph_tuples();
  
  // find 64 search keys
  uint65_t * roots = find_root_vertices(edges, num_root_vertices);
  
  // create graph data structure
  int i;
  for(i = 0; i < global_edge_count; i++) {
    graph.add_edge(
	  uint64_t(get_v0_from_edge(&edges[i])), 
	  uint64_t(get_v1_from_edge(&edges[i]))
	);
  }
  // commit local changes:
  // - edges / vertices that belong to other processors will get 
  //   redistributed
  // - vertices that do not exist get default constructed on the 
  //   target unit
  graph.commit();
  
  // perform BFS search
  for(i = 0; i < num_root_vertices) {
    uint65_t start_node = roots[i];
    std::list<dash::GraphNode> frontier = graph.find_local(start_node);
	std::list<dash::GraphNode> visited = frontier;
	std::list<dash::GraphNode> neighbours;
	std::list<dash::GraphNode> neighbours;
	std::list<dash::GraphNode> assigned;
	
	bool terminate_algorithm = false;
	while(!terminate_algorithm) {
	  if(!frontier.empty()) {
	    // frontier has to contain local elements only.
		// TODO: some of the neighbours retrieved are stored remotely
		//       -> no state can be exposed via GraphNode
		//       maybe return type of ID only (i.e. uint64_t)
		//       and use get_local_neighbours -> list of GraphNode
        std::list<dash::GraphNode>::const_iterator it;
        for(it = frontier.begin(); it != frontier.end(); it++) {
          neighbours.splice(graph.get_neighbours(it));
        }
		// alternative: neighbours = graph.get_neighbours(frontier);
		
		// perform "wave":
		std::list<dash::GraphNode> prefix;
		int partition_row = partition.get_row(myid);
		int partition_row_index = partition.get_row_index(myid);
		int partition_columns = partition.get_column_count();
		// we can infer more methods like get_column_index from this
		// maybe get_position -> (row, column)?
	    std::list<dash::GraphNode> forward = prefix;
		forward.splice(neighbours);

		//PROBLEM: how to perform "send/receive" in PGAS
		if(partition_row_index < partition_columns - 1) {
		  if(partition_row_index > 0) {
		    //wait for prefixes from previous processor in partition row
			//send prefix to next processor in partition row
		  } else {
	        //send prefix to next processor in partition row
		  }
		} else {
			//wait for prefixes from previous processor in partition row
		}
		
		delete forward;
		
		// TODO: assigned = neighbours - prefix;
		// TODO: broadcast neighbours to processors in row (see PROBLEM 1)
		// TODO: write predecessor map (uses assigned)
		// TODO: visited = visited + neighbours
		// TODO: check for termination (uses visited and neighbours)
		
		frontier = neighbours;
		neighbours.clear();
		assigned.clear();
	  }
	}
  }
}