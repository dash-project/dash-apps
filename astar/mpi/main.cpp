#include <mpi.h>
#include <stdio.h>
#include <vector>
#include <algorithm>
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
}

void addToQueue(std::vector<Puzzle> &pv, Puzzle p) {
  auto it = std::find(pv.begin(), pv.end(), p);
  if (it == pv.end()) {
    pv.push_back(p);
    return;
  }

  if ((*it).cost < p.cost) {
    (*it) = p;
  }
}

int main (int argc, char* argv[]) {

  MPI_Init(NULL, NULL);
  Puzzle p;
  MPI_Datatype puzzleType = register_type(p);

  std::vector<Puzzle> np;

  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);




  if (0 == world_rank) {
    //addToQueue(np, p);
    //addToQueue(np, np.back().moveRight());
    //addToQueue(np, np.back().moveDown());
    //addToQueue(np, np.back().moveLeft());
    //addToQueue(np, np.back().moveDown());
    //addToQueue(np, np.back().moveRight());
    np.push_back(p);

    Puzzle p2 = p.moveDown();
    np.push_back(p2);
    //addToQueue(np,p2);
    np.push_back(np.back().moveDown());

    std::cout << "----- VECTOR: -----" << std::endl;
    for (int i=0; i<np.size(); ++i) {
      np[i].print();
    }
    std::cout << "--------------- move tests fertig ------------------" << std::endl;
    sendPuzzle(p,puzzleType, 1, 0, MPI_COMM_WORLD);
    sendNPuzzles(np,puzzleType, 1, 0, MPI_COMM_WORLD);
  } else if(1 == world_rank) {
    recvPuzzle(p,puzzleType, 0, 0, MPI_COMM_WORLD);
    p.print();
    recvNPuzzles(np,puzzleType, 0, 0, MPI_COMM_WORLD);
    np[4].print();
  }

  deregister_mpi_type(puzzleType);
  MPI_Finalize();

  return 0;
}
