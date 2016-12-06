// ConsoleApplication1.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
#include "AStar.h"


int main(int argc, char *argv[])
{
  int rows = 3;
  int columns = 3;
  if (argc == 3) {
    std::cout << "rows: " << argv[1] << std::endl;
    std::cout << "columns: " << argv[2] << std::endl;
    rows = strtol(argv[1], NULL, 10);
    columns = strtol(argv[2], NULL, 10);
  }
  std::cout << "test" << std::endl;
  AStar a = AStar(rows,columns);
	
	a.solve();

	
	std::this_thread::sleep_for(std::chrono::seconds(5));

    return 0;
}

