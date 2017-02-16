#include <mpi.h>
#include <cstddef>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>

#include "Puzzle.h"





class Astar{
private:
	int world_rank, world_size, queue_ctr, exID = 0, flag = 0;
	bool done = false, first = true;
	
	std::vector<Puzzle> recv_buffer;
	std::vector<std::vector<Puzzle>> queues, send_buffers;
	std::vector<int> size_send_buffers, recv_sizes, sent_to_ctr, recv_from_ctr, noti_send_buffers, noti_recv_buffers;
	
	std::vector<MPI_Request>  size_send_requests,   size_recv_requests, 
                            puzzle_send_requests, puzzle_recv_requests,
                            noti_send_requests,   noti_recv_requests;

  MPI_Datatype puzzle_type;
	MPI_Status status;
	std::map<Puzzle, int, cmp> examined;
	
	Puzzle puzzle_buffer;

	void add_to_queue(Puzzle & p);
  void add_to_queue(Puzzle p, Puzzle & previous);	
  
	void handle_queue(int interrupt_ctr = std::numeric_limits<int>::max());
	
	void receive(int from);	
	void send(int to);
  
	void notify_others();
	void recv_notifications();
public:
	Astar();
	~Astar();
	
	int get_rank();
	void run(int break_after);
	void print(bool detailed = false);
  void print_all(bool detailed = false);
};
