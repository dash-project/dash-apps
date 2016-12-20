#include <mpi.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <algorithm>
#include "Puzzle.h"
#include <unistd.h>



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
/*
void sendPuzzle(Puzzle const& p, MPI_Datatype type, int dest, int tag, MPI_Comm comm) {
  MPI_Send(&p, 1, type, dest, tag, comm);
}

void sendNPuzzles(std::vector<Puzzle> puzzles, MPI_Datatype type, int dest, int tag, MPI_Comm comm) {
  int len = puzzles.size();
  MPI_Send(&len, 1, MPI_INT, dest, tag, comm);

  if (0 != len) {
    MPI_Send(&puzzles[0], puzzles.size(), type, dest, tag, comm);
  }
}

void recvPuzzle(Puzzle & p, MPI_Datatype type, int src, int tag, MPI_Comm comm) {
  MPI_Status s;
  MPI_Recv(&p, 1, type, src, tag, comm, &s);
}

void recvNPuzzles(std::vector<Puzzle> & puzzles, MPI_Datatype type, int src, int tag, MPI_Comm comm) {
  int len;
  MPI_Status s;
  MPI_Recv(&len, 1, MPI_INT, src, tag, comm, &s);

  if (len != 0) {
    puzzles.resize(len);
    MPI_Recv(puzzles.data(), len, type, src, tag, comm, &s);
  } else {
    puzzles.clear();
  }
}*/

void addToQueue(std::map<Puzzle,int,cmp> &exP, std::vector<std::vector<Puzzle>> &queues, Puzzle p, const int &processes) {
  auto it = exP.find(p);
  if (it != exP.end()) {
    if (p.cost < it->first.cost) {
      exP.erase(it);
      queues[p.get_responsible_process(processes)].push_back(p);
    }
    return;
  }

  queues[p.get_responsible_process(processes)].push_back(p);
}

int main (int argc, char* argv[]) {

  MPI_Init(NULL, NULL);

  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  Puzzle p;
  MPI_Datatype puzzleType = register_type(p);

  std::vector<std::vector<Puzzle>> queues;
  queues.resize(world_size);

  std::vector<int> size_send_buffers;
  size_send_buffers.resize(world_size);
  std::vector<std::vector<Puzzle>> send_buffers;
  send_buffers.resize(world_size);

  std::vector<int> recv_sizes;
  recv_sizes.resize(world_size);
  std::vector<Puzzle> recv_buffer;

  std::vector<MPI_Request> size_send_requests;
  size_send_requests.resize(world_size);
  std::vector<MPI_Request> puzzle_send_requests;
  puzzle_send_requests.resize(world_size);

  std::vector<MPI_Request> size_recv_requests;
  size_recv_requests.resize(world_size);
  std::vector<MPI_Request> puzzle_recv_requests;
  puzzle_recv_requests.resize(world_size);

  MPI_Status status;

  std::map<Puzzle, int, cmp> examined;
  int exID = 0;

  bool done = false, first = true;
  int flag = 0, localExaminedSum = 0, globalExaminedSum = 0, openOps = 0;

  if (0 == world_rank) {
    addToQueue(examined, queues, p, world_size);
  }

  //initialize the receives
  for (int i=0; i<world_size; ++i) {
      MPI_Irecv(&recv_sizes[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &size_recv_requests[i]);
      MPI_Isend(&size_send_buffers[i],1,MPI_INT, i, 0, MPI_COMM_WORLD, &size_send_requests[i]);
      MPI_Isend(&p, 1, puzzleType, i, 1, MPI_COMM_WORLD, &puzzle_send_requests[i]);
  }

  int recvFailures = 0;

  while (!done) {
    //handle own queue
    while (queues[world_rank].size() > 0) {
      p = queues[world_rank].back();
      queues[world_rank].pop_back();
      auto it = examined.find(p);

      if (it == examined.end() || it->first.cost > p.cost) {
        addToQueue(examined, queues, p.moveDown(),  world_size);
        addToQueue(examined, queues, p.moveUp(),    world_size);
        addToQueue(examined, queues, p.moveLeft(),  world_size);
        addToQueue(examined, queues, p.moveRight(),  world_size);
        examined[p] = ++exID;
      }
    }

    //communicate
    for (int i=0; i<world_size; ++i) {
      openOps = 0;
      if (i != world_rank) {
        //receive
        MPI_Test(&size_recv_requests[i], &flag, &status);
        if (flag) {
          //receive puzzle count
          MPI_Wait(&size_recv_requests[i],&status);
          if (recv_sizes[i] != 0) {
            //resize the buffer, receive and put the states in the queue
            recv_buffer.resize(recv_sizes[i]);
            MPI_Recv(recv_buffer.data(), recv_sizes[i], puzzleType,i,1,MPI_COMM_WORLD,&status);
            queues[world_rank].insert(queues[world_rank].end(), recv_buffer.begin(), recv_buffer.end());
            recv_buffer.clear();
            recvFailures = 0;
            MPI_Irecv(&recv_sizes[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &size_recv_requests[i]);
          } else {
            //first send/receive is only to make MPI_Test available
            recv_buffer.resize(1);
            MPI_Recv(recv_buffer.data(), 1, puzzleType,i,1,MPI_COMM_WORLD,&status);
            recv_buffer.clear();
            MPI_Irecv(&recv_sizes[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &size_recv_requests[i]);
          }
        }
        //send
        MPI_Test(&puzzle_send_requests[i], &flag, &status);
        if (queues[i].size() > 0 && flag) {
          size_send_buffers[i] = queues[i].size();
          send_buffers[i] = queues[i];
          queues[i].clear();
          MPI_Isend(&size_send_buffers[i],1,MPI_INT, i, 0, MPI_COMM_WORLD, &size_send_requests[i]);
          MPI_Isend(&send_buffers[i][0], size_send_buffers[i], puzzleType, i, 1, MPI_COMM_WORLD, &puzzle_send_requests[i]);
          recvFailures = 0;
          first = false;
        }
      }
    }

    ++recvFailures;
    if (recvFailures > 990) {
      sleep(1);
    }

    if (recvFailures > 1000) {
      done = true;
    }
  }

  localExaminedSum = examined.size();
  std::cout << "process " << world_rank << " is done, exmained "<< localExaminedSum << " states" << std::endl;
  MPI_Reduce(&localExaminedSum, &globalExaminedSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (world_rank == 0) {
    std::cout << "total states examined: " << globalExaminedSum << std::endl;
  }

  Puzzle p2;
  p2 = p2.randomize(500);
  auto it = examined.find(p2);
  if (it != examined.end()) {
    std::cout << "searching for puzzle:" << std::endl;
    p2.print();
    std::cout << "found the puzzle, cost is: " << it->first.cost << std::endl;
  }

  deregister_mpi_type(puzzleType);
  MPI_Finalize();

  return 0;
}
