#pragma once
#include <iostream>

const int COLUMNS = 2, ROWS = 3;

class Puzzle{
public:
  int cost;
  int puzzle[ROWS * COLUMNS], previous[ROWS * COLUMNS];
  int findEmpty();
  void puzzle_to_previous();

  Puzzle();

  inline bool operator==(const Puzzle& rhs) {
    return (puzzle == rhs.puzzle);
  }



  void print();
  Puzzle moveUp();
  Puzzle moveDown();
  Puzzle moveLeft();
  Puzzle moveRight();

  int get_responsible_process(int no_of_processes);
  
};
