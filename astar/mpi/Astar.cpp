#include "Astar.h"

Astar::Astar() {
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	puzzle_type = register_type(puzzle_buffer);
	
	queues.resize(world_size);
	send_buffers.resize(world_size);
	size_send_buffers.resize(world_size);
	recv_sizes.resize(world_size);
	
	sent_to_ctr.resize(world_size);
	recv_from_ctr.resize(world_size);
	noti_send_buffers.resize(world_size);
	noti_recv_buffers.resize(world_size);
	
	size_send_requests.resize(world_size);
	size_recv_requests.resize(world_size);
  puzzle_send_requests.resize(world_size);
	puzzle_recv_requests.resize(world_size);
  noti_send_requests.resize(world_size);
  noti_recv_requests.resize(world_size);
	
	if (0 == world_rank) {
		add_to_queue(puzzle_buffer);
		std::cout << "MPI started with " << world_size << " processes\n";
	}
	
	//initialize the receives
	for (int i=0; i<world_size; ++i) {
		MPI_Irecv(&recv_sizes[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &size_recv_requests[i]);
		MPI_Irecv(&noti_recv_buffers[i], 1, MPI_INT, i, 2, MPI_COMM_WORLD, &noti_recv_requests[i]);
		MPI_Isend(&size_send_buffers[i],1,MPI_INT, i, 0, MPI_COMM_WORLD, &size_send_requests[i]);
		MPI_Isend(&puzzle_buffer, 1, puzzle_type, i, 1, MPI_COMM_WORLD, &puzzle_send_requests[i]);
    MPI_Isend(&noti_send_buffers[i], 1, MPI_INT, i, 2, MPI_COMM_WORLD, &noti_send_requests[i]);
	}
}

Astar::~Astar() {
	deregister_mpi_type(puzzle_type);
}

void Astar::add_to_queue(Puzzle & p) {
	auto it = examined.find(p);
	if (it != examined.end()) {
		if (p.cost < it->first.cost) {
			examined.erase(it);
			queues[p.get_responsible_process(world_size)].emplace_back(p);
		}
		return;
	}
	queues[p.get_responsible_process(world_size)].emplace_back(p);
}
 
void Astar::add_to_queue(Puzzle p, Puzzle & previous) {
   if (p.puzzle != previous.previous) {
     add_to_queue(p);
   }
 }

void Astar::handle_queue(double interrupt_ctr) {
	queue_ctr = 0;
	while (queues[world_rank].size() > 0 && queue_ctr < interrupt_ctr) {
		puzzle_buffer = queues[world_rank].back();
		queues[world_rank].pop_back();
	
		auto it = examined.find(puzzle_buffer);
		if (it == examined.end() || it->first.cost > puzzle_buffer.cost) {
			add_to_queue(puzzle_buffer.moveDown(), puzzle_buffer);
			add_to_queue(puzzle_buffer.moveUp(), puzzle_buffer);
			add_to_queue(puzzle_buffer.moveLeft(), puzzle_buffer);
			add_to_queue(puzzle_buffer.moveRight(), puzzle_buffer);
			examined[puzzle_buffer] = ++exID;
		}
		
		++queue_ctr;
	}
}

void Astar::receive(int from) {
	MPI_Test(&size_recv_requests[from], &flag, &status);
   if (flag) {
     //receive puzzle count
     MPI_Wait(&size_recv_requests[from],&status);
     if (recv_sizes[from] != 0) {
       //resize the buffer, receive and put the states in the queue
       recv_buffer.resize(recv_sizes[from]);
       ++(recv_from_ctr[from]);
       MPI_Recv(recv_buffer.data(), recv_sizes[from], puzzle_type,from,1,MPI_COMM_WORLD,&status);
       queues[world_rank].insert(queues[world_rank].end(), recv_buffer.begin(), recv_buffer.end());
       recv_buffer.clear();
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

void Astar::send(int to) {
	MPI_Test(&puzzle_send_requests[to], &flag, &status);
   if (!queues[to].empty() && flag) {
     size_send_buffers[to] = queues[to].size();
		++(sent_to_ctr[to]);
		send_buffers[to] = queues[to];
		queues[to].clear();
		MPI_Isend(&size_send_buffers[to],1,MPI_INT, to, 0, MPI_COMM_WORLD, &size_send_requests[to]);
		MPI_Isend(&send_buffers[to][0], size_send_buffers[to], puzzle_type, to, 1, MPI_COMM_WORLD, &puzzle_send_requests[to]);
   }
}
/**********************************************************
* keep track how many times I sent states to a process. When
* notifying another process, tell them how many times I received. 
* so the notified / sending process is able to tell whether all
* sent states are processed.
**********************************************************/
//notify every thread i received from, that i processed
void Astar::notify_others() {
	for (int i=0; i<world_size; ++i) {
    if (i != world_rank) {
      MPI_Test(&noti_send_requests[i], &flag, &status);
      if (flag && recv_from_ctr[i] > 0) {
        //copy into send_buffer
        noti_send_buffers[i] = recv_from_ctr[i];
        //clear local (nothing is received from this process)
        recv_from_ctr[i] = 0;
        //send
        MPI_Isend(&noti_send_buffers[i], 1, MPI_INT, i, 2, MPI_COMM_WORLD, &noti_send_requests[i]);
      }
    }
	}
}

//check whether all sent states are processed or not!
void Astar::recv_notifications() {
	for (int i=0; i<world_size; ++i) {
    if (i != world_rank) {
      MPI_Test(&noti_recv_requests[i], &flag, &status);
      if (flag) {
        MPI_Wait(&noti_recv_requests[i],&status);
        
        //-1 is received, when p0 tells everyone, that it is done
        if (noti_recv_buffers[i] == -1 ) {
          done = true;
          break;
        }
        
        sent_to_ctr[i] -= noti_recv_buffers[i];
        noti_recv_buffers[i] = 0;
        
        //init the next receive
        MPI_Irecv(&noti_recv_buffers[i], 1, MPI_INT, i, 2, MPI_COMM_WORLD, &noti_recv_requests[i]);
      }
    }
  }
   
   if (world_rank == 0 && (std::accumulate(sent_to_ctr.begin(), sent_to_ctr.end(), 0) == 0)) {
     //std::cout << "looks like everyone is done!\n";
    done = true;
    for (int i=1; i<world_size; ++i) {
	    noti_send_buffers[i] = -1;
	    MPI_Isend(&noti_send_buffers[i], 1, MPI_INT, i, 2, MPI_COMM_WORLD, &noti_send_requests[i]);
    }
  }
}

int Astar::get_rank() {
	return world_rank;
}

void Astar::run(int break_after) {
	while (!done) {
		handle_queue(break_after);
    
    for (int i=0; i<world_size; ++i) {
			if (i != world_rank) {
				receive(i);
				send(i);
			}
		}
		//if i'm done, notify others that sent me states
    if (queues[world_rank].empty()) {
			//if (std::accumulate(sent_to_ctr.begin(), sent_to_ctr.end(), 0) == 0) { 
        notify_others();
      //}
      
      recv_notifications();
		}
		
	}
  int ex_ctr = examined.size();
  int recv;
  MPI_Reduce( &ex_ctr, &recv, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
  if (world_rank == 0) {
    std::cout << "total examined states: " << recv << '\n';
  }
}

void Astar::print(bool detailed) {
	std::cout << "---------- PROCESS " << world_rank << " ----------\n";
	std::cout << "examined " << examined.size() << " states!\n";
	
  if (detailed) {
    for (int i=0; i<world_size; ++i) {
      if (i != world_rank) {
        std::cout << "communication with " << i << "(s/r): "
                  << sent_to_ctr[i] << ", " << recv_from_ctr[i] << '\n';
      }
    }
  }
}

void Astar::print_all(bool detailed) {
  for (int i=0; i<world_size; ++i) {
    if (i == world_rank) {
      print(detailed);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

