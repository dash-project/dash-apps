#include <generator/make_graph.h>
#include <generator/utils.h>

#include <libdash.h>
#include <list>

/*
 * Creates tuples for graph
 */
packed_edge * create_graph_tuples() {
  int scale = 26;
  int edge_factor = 16;
  uint64_t seed1 = 2, 
    seed2 = 3;
  make_mrg_seed(seed1, seed2, seed);
  
  // number of edges = edge_factor * 2^scale
  packed_edge[local_edge__count] edges;
  generate_kronecker_range(seed, scale, start_edge, start_edge + local_edge_count, edges);
  return edges;
}

/*
 * Finds root vertices with min. 1 edge
 */
uint64_t * find_root_vertices(packed_edge & edges, int num_root_vertices) {
  std::array<uint64_t, num_root_vertices> roots;
  // TODO: to be implemented
  return roots;
}

/*
 * Removes prefix nodes from neighbour node list
 */
std::vector<dash::GraphNode> & get_assigned_nodes(std::vector<dash::GraphNode> & neighbours const, std::vector<dash::GraphNode> & prefix const) {
  std::vector<dash::GraphNode> output;
  // TODO: to be implemented
  return output;
}

/*
 * Adds neighbours to visited vector
 */
void add_to_visited(std::vector<dash::GraphNode> & visited, std::vector<dash::GraphNode> & neighbours const) {
  // TODO: to be implemented
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
  std::vector<dash::GraphNode> prefix;
  int partition_row = partition.get_row(myid);
  int partition_row_index = partition.get_column(myid);
  int partition_columns = partition.get_column_count();
  dart_unit_t next_unit = partition.next_in_row();
  dart_unit_t prev_unit = partition.prev_in_row();
  dart_unit_t transposed_unit = partition.get_unit(partition_row_index, partition_row); //(i, j) -> (j, i)
  dart_team_t row_team = partition.team_in_row();
  // we can infer more methods like get_row_count from this
  // maybe get_position -> (row, column)?
  // TODO: use vector of nodes rather than vector of IDs?
  std::vector<dash::GraphNode> frontier;
  std::vector<dash::GraphNode> neighbours;
  std::vector<dash::GraphNode> visited;
  std::vector<dash::GraphNode> assigned;
  for(i = 0; i < num_root_vertices) {
    uint65_t start_node = roots[i];
    frontier = graph.find_local(start_node);
    visited = frontier;
    
    while(1) {
      if(!frontier.empty()) {
        // frontier has to contain local elements only.
        std::vector<dash::GraphNode>::const_iterator it;
        for(auto el : frontier) {
          // TODO: save parent for each neighbour node:
          // should we save this directly in the GraphNode?
          std::vector<dash::GraphNode> current = graph.get_neighbours(el);
          neighbours.insert(neighbours.end(), current.begin(), current.end());
          // TODO: remove the nodes that are already contained in visited
        }
        // alternative: neighbours = graph.get_neighbours(frontier);
        
        // TODO: check for termination:
        // collective operation for all units: is the neighbour list empty?
        // so maybe neighbour list should be a global data structure
        // we could do a simple neighbours.empty() then
        
        // perform "wave":
        if(partition_row_index < partition_columns - 1) {
          if(partition_row_index > 0) {
            // wait for prefixes from previous processor in partition row
            dart_recv(static_cast<void *>(prefix), , prev_unit);
            // send prefix to next processor in partition row
            std::vector<dash::GraphNode> forward = prefix;
            forward.insert(forward.end(), neighbours.begin(), neighbours.end());
            dart_send(static_cast<void *>(forward), forward.size() * sizeof(dash::GraphNode), next_unit);
            delete forward;
          } else {
            // send prefix to next processor in partition row
            dart_send(static_cast<void *>(neighbours), neighbours.size() * sizeof(dash::GraphNode), next_unit);
          }
          assigned = get_assigned_nodes(neighbours, prefix);
          // wait for broadcast from last unit in row
          neighbours.clear();
          dart_team_recv(static_cast<void *>(neighbours), , row_team);
        } else {
          // wait for prefixes from previous processor in partition row
          dart_recv(static_cast<void *>(prefix), , prev_unit);
          assigned = get_assigned_nodes(neighbours, prefix);
          // this processor now holds the whole frontier for the next iteration
          // -> broadcast back to all other processors in row
          neighbours.insert(neighbours.end(), prefix.begin(), prefix.end());
          dart_team_send(static_cast<void *>(neighbours), neighbours.size() * sizeof(dash::GraphNode), row_team);
        }
        
        // each unit has 2 vectors now:
        // neighbours: frontier for the next iteration in this row
        // assigned: part of neighbours that this unit is responsible for
        
        // TODO: write predecessor map (using assigned)
        // every unit writes different amount of data into memory,
        // but the total size is fixed: do we need a dynamic data structure
        // or is there a pattern that could support this for dash::Array? 
        for(auto ass : assigned) {
          pred_map[ass.get_id()] = ass.get_parent();
        }
        
        add_to_visited(neighbours);
        
        // get neighbour list from unit of transposed partition:
        // frontier(i, j) = neighbours(j, i)
        dart_sendrecv(static_cast<void *>(neighbours), neighbours.size() * sizeof(dash::GraphNode), transposed_unit);
        frontier = neighbours;
        neighbours.clear();
        assigned.clear();
      }
    }
  }
}
