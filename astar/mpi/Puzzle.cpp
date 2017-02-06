#include "Puzzle.h"
#include <iostream>

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

Puzzle::Puzzle() {
  for (int i=0; i<ROWS * COLUMNS; ++i){
    puzzle[i] = i;
    previous[i] = i;
    cost = 0;
  }
}

Puzzle Puzzle::randomize(int n) {
  Puzzle p;
  for (int i=0; i<n; ++i) {
    int r = rand() % 4;
    switch (r) {
      case 0:
        p = p.moveRight();
        break;
      case 1:
        p = p.moveLeft();
        break;
      case 2:
        p = p.moveUp();
        break;
      case 3:
        p = p.moveDown();
        break;
      default:
        break;
    }
  }
  return p;
}

void Puzzle::print() {
  std::cout << "----- PUZZLE STATE: -----" << std::endl;
  std::cout << "cost: " << cost << std::endl;
  for (int i=0; i<ROWS; ++i) {
     for (int j=0; j<COLUMNS; ++j) {
        std::cout << puzzle[i*COLUMNS + j] << " ";
     }
     std::cout << std::endl;
  }

  std::cout << "  --- PREVIOUS STATE: ---  " << std::endl;
  for (int i=0;i<ROWS;++i) {
    for (int j=0;j<COLUMNS; ++j) {
      std::cout << previous[i*COLUMNS + j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "-------------------------" << std::endl;
}



Puzzle Puzzle::moveUp() {
  int empty = findEmpty();
  if (empty < COLUMNS) {
    return *this;
  }

  Puzzle p_new(*this);
  p_new.puzzle_to_previous();
  p_new.puzzle[empty] = p_new.puzzle[empty-COLUMNS];
  p_new.puzzle[empty - COLUMNS] = 0;
  p_new.cost++;

  return p_new;
}
Puzzle Puzzle::moveDown() {
  int empty = findEmpty();
  if ((empty / COLUMNS) == ROWS - 1) {
    return *this;
  }

  Puzzle p_new = Puzzle(*this);
  p_new.puzzle_to_previous();
  p_new.puzzle[empty] = p_new.puzzle[empty + COLUMNS];
  p_new.puzzle[empty + COLUMNS] = 0;
  p_new.cost++;

  return p_new;
}
Puzzle Puzzle::moveLeft() {
  int empty = findEmpty();
  if ((empty % COLUMNS)==0) {
    return *this;
  }

  Puzzle p_new = Puzzle(*this);
  p_new.puzzle_to_previous();
  p_new.puzzle[empty] = p_new.puzzle[empty-1];
  p_new.puzzle[empty-1] = 0;
  p_new.cost++;

  return p_new;
}
Puzzle Puzzle::moveRight() {
  int empty = findEmpty();
  if ((empty % COLUMNS) == COLUMNS-1) {
    return *this;
  }

  Puzzle p_new = Puzzle(*this);
  p_new.puzzle_to_previous();
  p_new.puzzle[empty] = p_new.puzzle[empty+1];
  p_new.puzzle[empty+1] = 0;
  p_new.cost++;

  return p_new;
}

int Puzzle::findEmpty() {
  for (size_t i = 0; i<ROWS*COLUMNS; ++i) {
    if (puzzle[i] == 0) {
      return i;
    }
  }
  return 0;
}

void Puzzle::puzzle_to_previous() {
  for (int i=0; i<ROWS*COLUMNS; ++i) {
    previous[i]=puzzle[i];
  }
}

int Puzzle::get_responsible_process(int no_of_processes) {
  //int area_size = SIZE / no_of_processes;
  //return findEmpty() / area_size ;
  return findEmpty() % no_of_processes;
}

