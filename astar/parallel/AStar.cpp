//#include "stdafx.h"
#include "AStar.h"

AStar::AStar() : 
	queueID(0),
  exID(0),
  found(false),
  p(3, 3)
{
}

AStar::AStar(int rowCount, int columnCount):
	queueID(0),
  exID(0),
  found(false),
  p(rowCount,columnCount)
{
}


AStar::~AStar()
{
}

void AStar::solve()
{
	goal = State(p.getGoalHash(), "-1", 0);
	std::cout << "goal: " << goal.state << std::endl;
	queuedStates.push_back(State(p.getStateHash(), "start", 0));

	p.printCurrentState();
	std::cout << "hash: " << p.getStateHash() << std::endl;

	while (!queuedStates.empty() && !found) {
		solveIterator = examinedStates.find(queuedStates.front());
		if (solveIterator == examinedStates.end()) {

			p.setState(queuedStates.front().state);
			p.moveUp();
			addToQueue(State(p.getStateHash(), queuedStates.front().state, queuedStates.front().cost + 1));

			p.setState(queuedStates.front().state);
			p.moveDown();
			addToQueue(State(p.getStateHash(), queuedStates.front().state, queuedStates.front().cost + 1));

			p.setState(queuedStates.front().state);
			p.moveLeft();
			addToQueue(State(p.getStateHash(), queuedStates.front().state, queuedStates.front().cost + 1));

			p.setState(queuedStates.front().state);
			p.moveRight();
			addToQueue(State(p.getStateHash(), queuedStates.front().state, queuedStates.front().cost + 1));

			examinedStates[queuedStates.front()] = ++exID;			
		}
		queuedStates.pop_front();
	}
	reconstructPath();
}

void AStar::reconstructPath()
{
	solveIterator = examinedStates.find(goal);
	if (solveIterator == examinedStates.end()) {
		std::cout << "there is no solution" << std::endl;
		return;
	}
	while (solveIterator != examinedStates.end()) {
		solution.push_front(solveIterator->first);
		solveIterator = examinedStates.find(State(solveIterator->first.prev, "egal", 0));
	}

	for (int i = 0; i < solution.size(); ++i) {
		p.setState(solution[i].state);
		p.printCurrentState();
	}
	std::cout << "The minimal solution takes " << solution.size()-1 << " steps!" << std::endl;
}

void AStar::addToQueue(State s)
{
	if (found) { return; }

	addIterator = examinedStates.find(s);
	if (addIterator != examinedStates.end()) {
		if (s.cost < addIterator->first.cost) {
			examinedStates.erase(addIterator);
			queuedStates.push_back(s);
		}
		return;
	}

	queuedStates.push_back(s);
}

State::State()
{
}

State::State(std::string s, std::string p, int c)
{
	state = s;
	cost = c;
	prev = p;
}

State::~State()
{
}

bool State::operator==(const State & b)
{
	return (0 == state.compare(b.state));
}

bool State::operator<(const State & b)
{
	return false;
}

bool State::operator()(const State & b)
{
	return false;
}
