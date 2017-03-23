#include "AStar.h"

AStar::AStar() {
	queue.push_back(Puzzle());
}

void AStar::add_to_queue(Puzzle & p) {
	auto it = examined.find(p);
	if (it != examined.end()) {
		if (p.cost < it->first.cost) {
			examined.erase(it);
			queue.emplace_back(p);
		}
		return;
	}
	queue.emplace_back(p);
}

void AStar::add_to_queue(Puzzle p, Puzzle & previous) {
  if (p.puzzle != previous.previous) {
    add_to_queue(p);
  }
}

void AStar::run() {
  JIpsManager im;
	while (queue.size() > 0) {
    im.update(true);
    //std::cout << "qs: " << queue.size() << '\n';
		puzzle_buffer = queue.front();
		queue.pop_front();
	
		auto it = examined.find(puzzle_buffer);
		if (it == examined.end() || it->first.cost > puzzle_buffer.cost) {
			add_to_queue(puzzle_buffer.moveDown(), puzzle_buffer);
			add_to_queue(puzzle_buffer.moveUp(), puzzle_buffer);
			add_to_queue(puzzle_buffer.moveLeft(), puzzle_buffer);
			add_to_queue(puzzle_buffer.moveRight(), puzzle_buffer);
			examined[puzzle_buffer] = ++id;
		}
	}
}

void AStar::print() {
	std::cout << "examined " << examined.size() << " states!\n";
}
