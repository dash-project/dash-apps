//#include "stdafx.h"
#include "Puzzle.h"

Puzzle::Puzzle(int rC, int cC)
{
  empty = Coords(0,0);
	rowCount = rC;
	columnCount = cC;
	createPuzzle();
	randomizePuzzle();
}


Puzzle::~Puzzle()
{
}

void Puzzle::createPuzzle()
{
	for (int i = 0; i < rowCount; i++) {
		std::vector<int> row;
		for (int j = 0; j < columnCount; j++) {
			row.push_back(i*columnCount + j);
		}
		puzzle.push_back(row);
	}
}

void Puzzle::randomizePuzzle() {
	for (int i = 0; i < 500; i++) {
		int r = rand() % 4;
		switch (r) {
		case 0:
			moveRight();
			break;
		case 1:
			moveLeft();
			break;
		case 2:
			moveUp();
			break;
		case 3:
			moveDown();
			break;
		default:
			break;
		}
	}
}

void Puzzle::swap(Coords a, Coords b)
{
  //std::cout << "a: " << a.row << ", " << a.column << std::endl;
  //std::cout << "b: " << b.row << ", " << b.column << std::endl;
  //std::cout << "dimensions: " << rowCount << ", " << columnCount << std::endl;
	int buf = puzzle[a.row][a.column];
	puzzle[a.row][a.column] = puzzle[b.row][b.column];
	puzzle[b.row][b.column] = buf;
}



std::string Puzzle::getGoalHash()
{
	std::string goal;
	for (int i = 0; i < rowCount*columnCount; i++) {
		goal += (char)(i+97);
	}
	return goal;
}

std::string Puzzle::getStateHash()
{
	std::string hash;
	for (int i = 0; i < puzzle.size(); i++) {
		for (int j = 0; j < puzzle[i].size(); j++) {
			hash += (char)(puzzle[i][j]+97);
		}
	}
	
	return hash;
}

void Puzzle::setState(std::string hash)
{
	for (int i = 0; i < hash.length(); i++) {
		puzzle[i / columnCount][i % columnCount] = (int)(hash[i] - 97);
		if ((int)(hash[i] - 97) == 0) {
			empty.column = i % columnCount;
			empty.row = i / columnCount;
		}
	}
}

void Puzzle::moveLeft()
{
	if (empty.column == 0) {
		return;
	}
	Coords newEmpty = Coords(empty.row, empty.column-1);
	swap(empty, newEmpty);
	empty = newEmpty;
}

void Puzzle::moveRight()
{
	if (empty.column == columnCount-1) {
		return;
	}
	Coords newEmpty = Coords(empty.row, empty.column+1);
	swap(empty, newEmpty);
	empty = newEmpty;
}

void Puzzle::moveUp()
{
	if (empty.row == 0) {
		return;
	}
	Coords newEmpty = Coords(empty.row-1, empty.column);
	swap(empty, newEmpty);
	empty = newEmpty;
}

void Puzzle::moveDown()
{
	if (empty.row == rowCount-1) {
		return;
	}
	Coords newEmpty = Coords(empty.row+1, empty.column);
	swap(empty, newEmpty);
	empty = newEmpty;
}



void Puzzle::printCurrentState()
{
	std::cout << "------ state -----" << std::endl;
	for (int i = 0; i < rowCount; i++) {
		for (int j = 0; j < columnCount; j++) {
			std::cout << puzzle[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "------------------" << std::endl;
}

Coords::Coords()
{
}

Coords::Coords(int _row, int _column)
{
	row		= _row;
	column	= _column;
}

Coords::~Coords()
{
}
