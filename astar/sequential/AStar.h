#pragma once
#include "Puzzle.h"
#include "iostream"
#include "vector"
#include "string"
#include "deque"
#include "map"
#include "thread"
#include "chrono"
#include "JUtils.h"

class AStar {
private:
	//std::vector<Puzzle> queue;
  std::deque<Puzzle> queue;
	std::map<Puzzle, int, cmp> examined;

  Puzzle puzzle_buffer;
	unsigned long id, ow_ctr;

	void add_to_queue(Puzzle & p);
  void add_to_queue(Puzzle p, Puzzle & previous);	

public:
	AStar();

	void run();
  void print();
	//void reconstructPath();
	
};
