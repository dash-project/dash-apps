#pragma once
#include <iostream>
#include <string.h>
#include <cstddef>
#include <mpi.h>

const int COLUMNS = 4, ROWS = 4;
//const int SIZE = COLUMNS * ROWS;

class Puzzle{
public:
  int cost;
  int puzzle[ROWS * COLUMNS], previous[ROWS * COLUMNS];
  int findEmpty();
  void puzzle_to_previous();

  Puzzle();

  inline bool operator==(const Puzzle& other) {
    return (memcmp(puzzle, other.puzzle, sizeof(puzzle)) == 0);
  }

  Puzzle randomize(int n);

  void print();
  Puzzle moveUp();
  Puzzle moveDown();
  Puzzle moveLeft();
  Puzzle moveRight();

  int get_responsible_process(int no_of_processes);
  
};

struct cmp{
  bool operator() (const Puzzle& lhs, const Puzzle& rhs) const {
    for (int i=0; i<ROWS*COLUMNS; ++i) {
      if (lhs.puzzle[i] < rhs.puzzle[i]) {
        return true;
      } else if (lhs.puzzle[i] > rhs.puzzle[i]) {
        return false;
      }
    }
    return false;
  }
};
