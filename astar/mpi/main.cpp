#include <mpi.h>
#include <stdio.h>
#include <vector>

class Puzzle{
  public:
  Puzzle(){
    columns = 3;
    rows = 3;
    for (int i=0; i< rows*columns; ++i) {
      //puzzle.push_back(i);
      //previous.push_back(i);
      puzzle[i] = i;
      previous[i] = i;
    }
  }

  int columns, rows, puzzle[9], previous[9];
  //std::vector<int> puzzle, previous;
  
  void print() {
    std::cout << "--- Puzzle ---" << std::endl;
    std::cout << puzzle[0] << ", " << puzzle[1] << ", " << puzzle[2] << std::endl;
    std::cout << puzzle[3] << ", " << puzzle[4] << ", " << puzzle[5] << std::endl;
    std::cout << puzzle[6] << ", " << puzzle[7] << ", " << puzzle[8] << std::endl;
  }
};

MPI_Datatype register_type(Puzzle const&) {
  size_t num_members = 4;
  int lengths[num_members] = {1,1,9,9};

  MPI_Aint offsets[num_members] = {
    offsetof(Puzzle, columns),
    offsetof(Puzzle, rows),
    offsetof(Puzzle, puzzle),
    offsetof(Puzzle, previous)};

  MPI_Datatype types[num_members] = { MPI_INT, MPI_INT, MPI_INT, MPI_INT };
  
  MPI_Datatype type;
  MPI_Type_create_struct(num_members, lengths, offsets, types, &type);
  MPI_Type_commit(&type);
  return type;
}

void deregister_mpi_type(MPI_Datatype type) {
  MPI_Type_free(&type);
}

void sendPuzzle(Puzzle const& p, int dest, int tag, MPI_Comm comm) {
  MPI_Datatype type = register_type(p);
  MPI_Send(&p, 1, type, dest, tag, comm);
  deregister_mpi_type(type);
}

void sendNPuzzles(std::vector<Puzzle> puzzles, int dest, int tag, MPI_Comm comm) {
  int len = puzzles.size();
  MPI_Send(&len, 1, MPI_INT, dest, tag, comm);

  if (0 != len) {
    MPI_Datatype type = register_type(puzzles[0]);
    MPI_Send(&puzzles[0], puzzles.size(), type, dest, tag, comm);
    deregister_mpi_type(type);
  }
}

void recvPuzzle(Puzzle & p, int src, int tag, MPI_Comm comm) {
  MPI_Status s;
  MPI_Datatype type = register_type(p);
  MPI_Recv(&p, 1, type, src, tag, comm, &s);
  deregister_mpi_type(type);
}

void recvNPuzzles(std::vector<Puzzle> & puzzles, int src, int tag, MPI_Comm comm) {
  int len;
  MPI_Status s;
  MPI_Recv(&len, 1, MPI_INT, src, tag, comm, &s);

  if (len != 0) {
    puzzles.resize(len);
    MPI_Datatype type = register_type(puzzles[0]);
    MPI_Recv(puzzles.data(), len, type, src, tag, comm, &s);
    deregister_mpi_type(type);
  } else {
    puzzles.clear();
  }
}

int main (int argc, char* argv[]) {
  std::vector<int> p1 = {0,1,2,3,4,5,6,7,8};
  std::vector<int> p2 = {8,7,6,5,4,3,2,1,0};

  MPI_Init(NULL, NULL);

  std::vector<Puzzle> np;
  Puzzle p;
  np.push_back(p);
  np.push_back(Puzzle());
  np.push_back(Puzzle());
  np.push_back(Puzzle());
  np.push_back(Puzzle());
  np.push_back(Puzzle());

  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  if (0 == world_rank) {
    p.puzzle[3] = 21;
    np[4].puzzle[6] = 139;
    sendPuzzle(p, 1, 0, MPI_COMM_WORLD);
    sendNPuzzles(np, 1, 0, MPI_COMM_WORLD);
    //MPI_Send(&p, 1, ptype, 1, 0, MPI_COMM_WORLD);
  } else if(1 == world_rank) {
    //MPI_Recv(&p, 1, ptype, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    recvPuzzle(p, 0, 0, MPI_COMM_WORLD);
    p.print();
    recvNPuzzles(np, 0, 0, MPI_COMM_WORLD);
    np[4].print();
  }

  MPI_Finalize();

  return 0;
}
