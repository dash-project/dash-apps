#pragma once
#include "Puzzle.h"
#include "iostream"
#include "vector"
#include "string"
#include "deque"
#include "map"
#include "thread"
#include "chrono"


class State
{
private:
public:
	State();
	State(std::string s, std::string p, int c);

	~State();

	bool operator==(const State & b);
	bool operator<(const State & b);
	bool operator()(const State & b);

	std::string state, prev;
	int cost;
};

struct cmp
{
	bool operator() (const State& lhs, const State& rhs) const
	{
		for (int i = 0; i < lhs.state.size(); ++i) {
			if (lhs.state[i] < rhs.state[i]) {
				return true;
			}
			else if (lhs.state[i] > rhs.state[i]) {
				return false;
			}
		}
		return false;
	}
};



class AStar
{
private:
	//std::deque<State> examinedStates, queuedStates, solution;
	std::deque<State> queuedStates, solution;
	std::map<State, int, cmp> examinedStates; //, queuedStates, solution;
	std::map<State, int, cmp>::iterator solveIterator, addIterator;
	State goal;
	unsigned long queueID, exID;

	void addToQueue(State s);
	bool found;

	//muss dann durch ein template ersetzt werden
	Puzzle p;
	
public:
	AStar();
	AStar(int rowCount, int columnCount);
	~AStar();

	void solve();
	void reconstructPath();
	
};
