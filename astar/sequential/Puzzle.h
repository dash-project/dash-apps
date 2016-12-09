#pragma once
#include "string"
#include "vector"
#include "iostream"

class Coords {
public:
	Coords();
	Coords(int r, int c);
	~Coords();

	int column, row;
};

class Puzzle
{
private:
	int rowCount, columnCount;
	Coords empty;
	std::vector<std::vector<int>> puzzle;
	
	void createPuzzle();
	void randomizePuzzle();
	void swap(Coords a, Coords b);

public:
	Puzzle(int rC, int cC);
	~Puzzle();


	void moveUp();
	void moveDown();
	void moveLeft();
	void moveRight();

	std::string getGoalHash();
	std::string getStateHash();
	void setState(std::string hash);
	void printCurrentState();
};

