//#pragma once
#include <mpi.h>
#include <cstddef>
#include <vector>
#include <map>

#include "Puzzle.h"

MPI_Datatype register_type(Puzzle const&) {
	size_t num_members = 3;
	int lengths[num_members] = {1,ROWS*COLUMNS,ROWS*COLUMNS};
	
	MPI_Aint offsets[num_members] = {
		offsetof(Puzzle, cost),
		offsetof(Puzzle, puzzle),
		offsetof(Puzzle, previous)};
	
	MPI_Datatype types[num_members] = { MPI_INT, MPI_INT, MPI_INT };
	
	MPI_Datatype type;
	MPI_Type_create_struct(num_members, lengths, offsets, types, &type);
	MPI_Type_commit(&type);
	return type;
}

void deregister_mpi_type(MPI_Datatype type) {
	MPI_Type_free(&type);
}

class Astar{
private:
	int world_rank, world_size, queue_ctr, exID = 0, flag = 0, localExaminedSum = 0, globalExaminedSum = 0, recvFailures = 0;
	bool done = false, first = true;
	
	std::vector<Puzzle> recv_buffer;
	std::vector<std::vector<Puzzle>> queues, send_buffers;
	std::vector<int> size_send_buffers, recv_sizes;
	
	std::vector<MPI_Request> size_send_requests, puzzle_send_requests, size_recv_requests, puzzle_recv_requests;
	
	MPI_Datatype puzzle_type;
	MPI_Status status;
	std::map<Puzzle, int, cmp> examined;
	
	Puzzle puzzle_buffer;

	void add_to_queue(Puzzle p) {
		auto it = examined.find(p);
		if (it != examined.end()) {
			if (p.cost < it->first.cost) {
			examined.erase(it);
			queues[p.get_responsible_process(world_size)].push_back(p);
			}
			return;
		}
		queues[p.get_responsible_process(world_size)].push_back(p);
	}
	
	void handle_queue(double interrupt_ctr = std::numeric_limits<double>::infinity()) {
		queue_ctr = 0;
		while (queues[world_rank].size() > 0 && queue_ctr < interrupt_ctr) {
			puzzle_buffer = queues[world_rank].back();
			queues[world_rank].pop_back();
		
			auto it = examined.find(puzzle_buffer);
			if (it == examined.end() || it->first.cost > puzzle_buffer.cost) {
				add_to_queue(puzzle_buffer.moveDown());
				add_to_queue(puzzle_buffer.moveUp());
				add_to_queue(puzzle_buffer.moveLeft());
				add_to_queue(puzzle_buffer.moveRight());
				examined[puzzle_buffer] = ++exID;
			}
			
			++queue_ctr;
		}
	}
	
	void receive(int from) {
		MPI_Test(&size_recv_requests[from], &flag, &status);
        if (flag) {
			//receive puzzle count
			MPI_Wait(&size_recv_requests[from],&status);
			if (recv_sizes[from] != 0) {
				//resize the buffer, receive and put the states in the queue
				recv_buffer.resize(recv_sizes[from]);
				MPI_Recv(recv_buffer.data(), recv_sizes[from], puzzle_type,from,1,MPI_COMM_WORLD,&status);
				queues[world_rank].insert(queues[world_rank].end(), recv_buffer.begin(), recv_buffer.end());
				recv_buffer.clear();
				recvFailures = 0;
				MPI_Irecv(&recv_sizes[from], 1, MPI_INT, from, 0, MPI_COMM_WORLD, &size_recv_requests[from]);
			} else {
				//first send/receive is only to make MPI_Test available
				recv_buffer.resize(1);
				MPI_Recv(recv_buffer.data(), 1, puzzle_type,from,1,MPI_COMM_WORLD,&status);
				recv_buffer.clear();
				MPI_Irecv(&recv_sizes[from], 1, MPI_INT, from, 0, MPI_COMM_WORLD, &size_recv_requests[from]);
			}
        }
	}
	
	void send(int to) {
		MPI_Test(&puzzle_send_requests[to], &flag, &status);
        if (!queues[to].empty() && flag) {
			size_send_buffers[to] = queues[to].size();
			send_buffers[to] = queues[to];
			queues[to].clear();
			MPI_Isend(&size_send_buffers[to],1,MPI_INT, to, 0, MPI_COMM_WORLD, &size_send_requests[to]);
			MPI_Isend(&send_buffers[to][0], size_send_buffers[to], puzzle_type, to, 1, MPI_COMM_WORLD, &puzzle_send_requests[to]);
			recvFailures = 0;
        }
	}
public:
	Astar() {
		MPI_Comm_size(MPI_COMM_WORLD, &world_size);
		MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
		
		puzzle_type = register_type(puzzle_buffer);
		
		queues.resize(world_size);
		send_buffers.resize(world_size);
		size_send_buffers.resize(world_size);
		recv_sizes.resize(world_size);
		
		size_send_requests.resize(world_size);
		puzzle_send_requests.resize(world_size);
		size_recv_requests.resize(world_size);
		puzzle_recv_requests.resize(world_size);
		
		if (0 == world_rank) {
			add_to_queue(puzzle_buffer);
			std::cout << "MPI started with " << world_size << " processes\n";
		}
		
		//initialize the receives
		for (int i=0; i<world_size; ++i) {
			MPI_Irecv(&recv_sizes[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &size_recv_requests[i]);
			MPI_Isend(&size_send_buffers[i],1,MPI_INT, i, 0, MPI_COMM_WORLD, &size_send_requests[i]);
			MPI_Isend(&puzzle_buffer, 1, puzzle_type, i, 1, MPI_COMM_WORLD, &puzzle_send_requests[i]);
		}
	}
	
	~Astar() {
		deregister_mpi_type(puzzle_type);
	}
	
	void run() {
		for (int a=0;a<1000000; ++a) {
		//while (!done) {
			handle_queue();
			for (int i=0; i<world_size; ++i) {
				if (i != world_rank) {
					receive(i);
					send(i);
				}
			}
		}
		std::cout << "process " << world_rank << " examined " << examined.size() << " states!\n";
	}
};